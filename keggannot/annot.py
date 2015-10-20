import os, sys
from decimal import Decimal
import gzip
from collections import OrderedDict
import logging

def blast_result_iterator(blast_file):
    _ = blast_file.readline()  # burn header
    for line in blast_file:
        if not line.startswith("#"):
            yield BlastHit(line)

class BlastHit(object):
    """
    Parse ncbi blast tabular output, comments will be ignored.

    blastp produces the following columns:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    """

    col_sep = "\t"
    query_id_col = 0
    subject_id_col = 1
    evalue_col = 10
    bitscore_col = 11
    
    def __init__(self, line=None):
        self.line = line
        self.query_id = None
        self.subject_id = None
        self.evalue = None
        self.bitscore = None
        if self.line:
            cols = line.split(self.col_sep)
            self.query_id = cols[self.query_id_col].strip()
            self.subject_id = cols[self.subject_id_col].strip()
            self.evalue = Decimal(cols[self.evalue_col].strip())
            self.bitscore = Decimal(cols[self.bitscore_col].strip())
    
    def __repr__(self):
        vals = (self.__class__.__name__, self.query_id, self.subject_id,
                self.evalue, self.bitscore)
        return "%s(query_id=%s, subject_id=%s, evalue=%s, bitscore=%s)" % vals

class Gene(object):
    def __init__(self, gene_id):
        self.id = gene_id
        self.best_hit = None
        # the best hit with KO was not used because it lacked a module
        # -m flag changed this result
        self.module_skip = False
    
    def add_blast_hit(self, hit):
        logging.debug("In %s, adding blast hit %s" % (self, hit))
        if (not self.best_hit) or (hit.evalue < self.best_hit.evalue):
            logging.debug("New best hit")
            self.best_hit = hit
    
    def __repr__(self):
        vals = (self.__class__.__name__, self.id, self.best_hit, self.module_skip)
        return "%s(id=%s, best_hit=%s, module_skip=%s)" % vals

class KEGGAnnotator(object):
    def __init__(self, ko_genes_list_file, ko_module_list_file,
                 ko_pathway_list_file, ko_enzyme_list_file, ko_file,
                 module_file, pathway_file, blast_file=None,
                 require_module=False):

        # To be output, a hit must have a module definition.  The first (best)
        # hit with a module definition is used for each query.
        if require_module:
            self.add_blast_hit = self.add_blast_hit_module
        else:
            self.add_blast_hit = self.add_blast_hit_ko

        # Load KO lookup files found in genes/ko/ko_*.list
        logging.info("Creating new %s object" % self.__class__.__name__)
        logging.info("Loading ko_genes.list file")
        self.ko2genes = self.dict_from_tabular_file(ko_genes_list_file, 
                                                    remove_key_prefix=True, remove_val_prefix=False)
        logging.info("Done loading ko_genes.list file, found %i KOs" % len(self.ko2genes))
        logging.info("Loading ko_module.list file")
        self.ko2modules = self.dict_from_tabular_file(ko_module_list_file)
        logging.info("Done loading ko_module.list file, found %i KOs" % len(self.ko2modules))

        logging.info("Loading ko_pathway.list file")
        self.ko2pathways = self.dict_from_tabular_file(ko_pathway_list_file)
        logging.info("Done loading ko_pathway.list file, found %i KOs" % len(self.ko2pathways))

        logging.info("Loading ko_enzyme.list file")
        self.ko2enzymes = self.dict_from_tabular_file(ko_enzyme_list_file)
        logging.info("Done loading ko_enzyme.list file, found %i KOs" % len(self.ko2enzymes))

        # Make genes to KO lookup
        logging.info("Creating gene to KOs lookup table")
        self.gene2kos = dict()
        for ko, gene_list in self.ko2genes.items():
            for g in gene_list:
                klist = self.gene2kos.setdefault(g, [])
                klist.append(ko)
        logging.info("Done creating gene to KOs lookup table")

        # Load KO, module, pathway definition file
        logging.info("Loading ko file")
        self.kos = self.parse_ko_file(ko_file)
        logging.info("Done loading ko file")
        logging.info("Loading module file")
        self.modules = self.parse_module_file(module_file)
        logging.info("Done loading module file")
        logging.info("Loading pathway file")
        self.pathways = self.parse_pathway_file(pathway_file)
        logging.info("Done loading pathway file")

        # keep track of genes which had ko, useful to mark cases where -m changed
        # results.
        self.genes_with_ko = {}

        self.genes = {}
        if blast_file:
            logging.info("Adding new blast results file")
            self.add_blast_results(blast_file)
            logging.info("Done adding new blast results")

        
    
    def dict_from_tabular_file(self, f, key_column=0, value_column=1, 
                               remove_key_prefix=True, remove_val_prefix=True, 
                               comment_symbol=None):
        d = dict()
        max_column = max(key_column, value_column)
        line_num = 0
        for line in f:
            line_num += 1
            if comment_symbol and line[0] == comment_symbol:
                continue
            fields = line.split()
            if max_column > len(fields):
                logging.error("Highest column is greater than columns in line %i" % line_num)
                sys.exit(1)
            # Remove ko:, md:, etc. if required
            if remove_key_prefix:
                fields[key_column] = self.remove_prefix(fields[key_column])
            if remove_val_prefix:
                fields[value_column] = self.remove_prefix(fields[value_column])
            
            if fields[key_column] in d:
                #logging.debug("Duplicate key %s encountered" % fields[key_column])
                d[fields[key_column]].append(fields[value_column])
            else:
                d[fields[key_column]] = [fields[value_column]]
        return d
    
    def remove_prefix(self, kegg_id):
        """Remove ko:, md:, path:, ec: from ID"""
        return kegg_id.split(":")[1]
    
    def parse_ko_file(self, ko_file):
        """Parse genes/ko/ko file"""
        k = dict()
        cur_entry = ""
        cur_name = ""
        cur_def = ""
        line_no = 0
        for line in ko_file:
            line_no += 1
            if line.startswith("ENTRY"):
                cur_entry = line.split()[1].rstrip()
            elif line.startswith("NAME"):
                cur_name = line.split(None, 1)[1].rstrip()
            elif line.startswith("DEFINITION"):
                cur_def = line.split(None, 1)[1].rstrip()
            elif line.startswith("///"):
                if (not cur_entry) or (not cur_name):
                    sys.stderr.write("Error parsing %s at %i\n" % (ko_file.name, line_no))
                    sys.exit(1)
                k[cur_entry] = {"name": cur_name, "def": cur_def}
                cur_entry = ""
                cur_name = ""
                cur_def = ""
        return k
    
    def parse_module_file(self, module_file):
        """Parse module/module file"""
        m = dict()
        cur_entry = ""
        cur_name = ""
        cur_class = ""
        line_no = 0
        for line in module_file:
            line_no += 1
            if line.startswith("ENTRY"):
                cur_entry = line.split()[1].rstrip()
            elif line.startswith("NAME"):
                cur_name = line.split(None, 1)[1].rstrip()
            elif line.startswith("CLASS"):
                cur_class = line.split(None, 1)[1].rstrip()
            elif line.startswith("///"):
                if (not cur_entry) or (not cur_name):
                    sys.stderr.write("Error parsing %s at %i\n" % (module_file.name, line_no))
                    sys.exit(1)
                m[cur_entry] = {"name": cur_name, "class": cur_class}
                cur_entry = ""
                cur_name = ""
                cur_class = ""
        return m
    
    def parse_pathway_file(self, pathway_file):
        """Parse pathway/pathway file"""
        p = dict()
        cur_name = ""
        cur_desc = ""
        cur_class = ""
        cur_entry = ""
        line_no = 0
        for line in pathway_file:
            line_no += 1
            if line.startswith("ENTRY"):
                fields = line.split()
                if fields[1].startswith("k"):
                    cur_entry = fields[1]
            elif cur_entry and line.startswith("NAME"):
                cur_name = line.split(None, 1)[1].rstrip()
            elif cur_entry and line.startswith("DESCRIPTION"):
                cur_desc = line.split(None, 1)[1].rstrip()
            elif cur_entry and line.startswith("CLASS"):
                cur_class = line.split(None, 1)[1].rstrip()
            elif line.startswith("///"):
                if cur_entry:
                    if (not cur_name) and (not cur_desc) and (not cur_class):
                        sys.stderr.write("Error parsing %s at %i\n" % (pathway_file.name, line_no))
                        sys.exit(1)
                    p[cur_entry] = {"name": cur_name, "desc": cur_desc, "class": cur_class}
                cur_pathway = ""
                cur_name = ""
                cur_desc = ""
                cur_class = ""
                cur_entry = ""
        return p

    def add_blast_results(self, blast_file):
        cnt = 0
        for hit in blast_result_iterator(blast_file):
            self.add_blast_hit(hit)
            cnt += 1
        logging.info("Found %i blast results, %i passed filter" % (cnt, len(self.genes)))
    
    # Only add hit if it has KO annotation
    def add_blast_hit_ko(self, hit):
        logging.debug("Attempting to add hit: %s" % hit)
        qid = hit.query_id
        if (not qid in self.genes) and self.get_hit_ko_list(hit):
            self.genes[qid] = Gene(qid)
            self.genes[qid].add_blast_hit(hit)
    
    # Only add hit if it has module annotation
    def add_blast_hit_module(self, hit):
        logging.debug("Attempting to add hit: %s" % hit)
        qid = hit.query_id
        if (not qid in self.genes):
            kos = self.get_hit_ko_list(hit)
            if kos:
                if self.get_hit_module_list(kos):
                    # This hit has module annotations so add it
                    g = Gene(qid)
                    self.genes[qid] = g
                    if qid in self.genes_with_ko:
                        g.module_skip = True   # note that this results differs because of -m
                    g.add_blast_hit(hit)       # add the hit result to the Gene object
                else:
                    # This hit didn't have module annotations, but had ko so add to 
                    # genes_with_ko dict
                    self.genes_with_ko[qid] = True

    def get_hit_ko_list(self, hit):
        if not hit:
            return None
        try:
            return sorted(self.gene2kos[hit.subject_id])
        except KeyError:
            logging.debug("KOs not found for gene %s" % hit.subject_id)
            return None
    
    def get_hit_module_list(self, kos):
        answer = set()
        if kos:
            for ko in kos:
                if ko in self.ko2modules:
                    answer.update(self.ko2modules[ko])
        if len(answer):
            return sorted(answer)
        else:
            if kos:
                logging.debug("Modules not found for KOs %s" % kos)
            return None
    
    def get_hit_pathway_list(self, kos):
        answer = set()
        if kos:
            for ko in kos:
                try:
                    answer.update(self.ko2pathways[ko])
                except KeyError:
                    pass
        if len(answer):
            return sorted(answer)
        else:
            if kos:
                logging.debug("Pathways not found for KOs %s" % kos)
            return None
    
    def get_hit_enzyme_list(self, kos):
        answer = set()
        if kos:
            for ko in kos:
                try:
                    answer.update(self.ko2enzymes[ko])
                except KeyError:
                    pass
        if len(answer):
            return sorted(answer)
        else:
            if kos:
                logging.debug("Enzymes not found for KOs %s" % kos)
            return None
    
    def make_basic_report_text(self):
        yield "\t".join(["query", "gene", "KO", "KO_names", "KO_descriptions",
                         "modules", "module_names", "module_classes", "pathways", "pathway_names",
                         "pathway_classes", "EC", "evalue", "score", "module_skip"])
        annotations = self.get_gene_annotations()
        for gene_id, annot in annotations.iteritems():
            yield self.tabify_annotations(gene_id, annot)
    
    def get_gene_annotations(self):
        out = dict()
        for gene_id, gene in self.genes.iteritems():
            out[gene_id] = self.aggregate_hit_data(gene)
        return out
    
    def aggregate_hit_data(self, gene):
        hit = gene.best_hit
        out = {"gene": "",
               "kos": [],
               "ko_names": [],
               "ko_defs": [],
               "modules": [],
               "module_names": [],
               "module_classes": [],
               "pathways": [],
               "pathway_names": [],
               "pathway_classes": [],
               "enzymes": [],
               "evalue": "",
               "score": "",
               "module_skip": False
        }
        
        # Add hit gene
        if hit:
            out["gene"] = hit.subject_id
            out["evalue"] = hit.evalue
            out["score"] = hit.bitscore
            out["module_skip"] = gene.module_skip
        
        kos = self.get_hit_ko_list(hit)
        if len(kos) > 1:
            logging.info("More than one KO for %s" % gene.id)
        out.update(self.ko_annots(kos))
        return out
    
    def ko_annots(self, kos):
        """Make a dictionary of KEGG annotations for a list of KOs"""
        out = {}
        if kos:
            out["kos"] = kos
            out["ko_names"] = [self.kos[k]["name"] for k in kos]
            out["ko_defs"] = [self.kos[k]["def"] for k in kos]
        
        # Build module list
        modules = self.get_hit_module_list(kos)
        if modules:
            out["modules"] = modules
            out["module_names"] = [self.modules[m]["name"] for m in modules]
            out["module_classes"] = [self.modules[m]["class"] for m in modules]
        
        # Build pathway list
        pathways = self.get_hit_pathway_list(kos)
        if pathways:
            out["pathways"] = pathways
            out["pathway_names"] = [self.pathways[p]["name"] for p in pathways]
            out["pathway_classes"] = [self.pathways[p]["class"] for p in pathways]
        
        # Build enzyme list
        enzymes = self.get_hit_enzyme_list(kos)
        if enzymes:
            out["enzymes"] = enzymes
        
        return out
    
    def tabify_annotations(self, gene_id, annot):
        """Create tab delimited string for hit data created by aggregate_hit_data"""
        out_text = gene_id + "\t"
        out_text += "\t".join([annot["gene"],
                               self.tabify_ko_annotations(annot),
                               str(annot["evalue"]), str(annot["score"]),
                               str(annot["module_skip"])])
        return out_text
    
    def tabify_ko_annotations(self, annot):
        """Tabify output from ko_annots"""
        return "\t".join(["; ".join(annot["kos"]),
                          "; ".join(annot["ko_names"]),
                          "; ".join(annot["ko_defs"]),
                          "; ".join(annot["modules"]),
                          "; ".join(annot["module_names"]),
                          " ::: ".join(annot["module_classes"]),
                          "; ".join(annot["pathways"]),
                          "; ".join(annot["pathway_names"]),
                          " ::: ".join(annot["pathway_classes"]),
                          "; ".join(annot["enzymes"])])
