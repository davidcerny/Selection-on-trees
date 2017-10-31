import os, uuid
from Bio.Phylo.PAML import codeml

results = codeml.read("/Users/David/Documents/Extracurricular/2016-17/CRG_Summer_Internship/Research/amniote_c-mos_fixedbrlen.mlc")

print results["dS tree"]
