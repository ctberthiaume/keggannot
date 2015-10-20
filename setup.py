from distutils.core import setup
import re

version_file = "keggannot/_version.py"
try:
    version_text = open(version_file).read()
except IOError:
    raise RuntimeError("Error: unable to open version file %s" % version_file)
version_re = re.compile(r"__version_info__\s*=\s*\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^,]+)\s*\)")
m = version_re.search(version_text)
if m:
    version_str = ".".join(m.groups())
else:
    raise RuntimeError("Error: unable to parse version string from file %s" % version_file)

setup(name="Keggannot",
    packages=["keggannot"],
    version=version_str,
    description="KEGG annotation based on protein BLAST search",
    author="Chris Berthiaume",
    author_email="chrisbee@uw.edu",
    scripts=["bin/keggannot_genes2ko"])
