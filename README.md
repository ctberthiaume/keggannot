Package to annotate proteins based on blastp against the KEGG database.


## Install

    python setup.py install


## Usage

    usage: keggannot_genes2ko [-h] [-v] [-V] [-m] [kegg_dir] [blast_file]

    Annotate proteins with KEGG BLAST results

    positional arguments:
      kegg_dir       KEGG data directory (default: None)
      blast_file     Tabular blast results file (default: None)

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  Print version to STDOUT (default: False)
      -V, --verbose  Verbose logging output (default: False)
      -m, --module   Only use hits with a defined module (default: False)

    The following files should be available in the KEGG data directory:
    genes/ko/ko_genes.list genes/ko/ko_module.list genes/ko/ko_pathway.list
    genes/ko/ko_enzyme.list genes/ko/ko module/module pathway/pathway


## KEGG Files needed

The kegg_dir argument should point to a copy of the KEGG data directory that
contains the following files

* genes/ko/ko_genes.list
* genes/ko/ko_module.list
* genes/ko/ko_pathway.list
* genes/ko/ko_enzyme.list
* genes/ko/ko
* module/module
* pathway/pathway


## Blast tabular input

The blast_file argument should be the result of a blastp agains the KEGG
protein database.  Format of the file should be either blast+ -outfmt 6 or 
-outfmt 7 (tabular without or with comments).


## Output

By default, output shows KO annotations for query sequences.  If the best hit
did not have a KO annotation, hits farther down the list are checked until a
hit is found with KO.  When -m is used, the same rules apply except the first
hit with a module annotation is used.  This is useful when summarizing results
by module class level 3.  Often the first KO annotation for a given query does
not have a module defined, but KOs for later hits may have modules.  While the
specific annotations for these later hits may not be as accurate as the best
hit with a KO, the level 3 module class annotation is often good enough.
