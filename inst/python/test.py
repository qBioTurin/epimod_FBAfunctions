# inst/python/test.py

import cobra, sys, json

if len(sys.argv) != 3:
    print("Usage: python test.py input.mat output.json")
    sys.exit(1)

inp, out = sys.argv[1], sys.argv[2]
model = cobra.io.load_matlab_model(inp)

# 1) Build index maps
mets = list(model.metabolites)
rxns = list(model.reactions)
genes = list(model.genes)
met_index  = {m.id: i for i, m in enumerate(mets)}
rxn_index  = {r.id: j for j, r in enumerate(rxns)}
gene_index = {g.id: k for k, g in enumerate(genes)}

# 2) Stoichiometry entries
S_i, S_j, S_x = [], [], []
for r in rxns:
    j = rxn_index[r.id]
    for met, coef in r.metabolites.items():
        i = met_index[met.id]
        S_i.append(i)
        S_j.append(j)
        S_x.append(coef)

# 3) Gene–reaction association matrix
#    rows = reactions, cols = genes
rxnGeneMat = [[0]*len(genes) for _ in rxns]
for r in rxns:
    j = rxn_index[r.id]
    for g in r.genes:            # g is a Gene object
        k = gene_index[g.id]
        rxnGeneMat[j][k] = 1

# 4) Pack into dict
d = {
    "S_i"       : S_i,
    "S_j"       : S_j,
    "S_x"       : S_x,
    "dims"      : [len(mets), len(rxns)],
    "lb"        : [r.lower_bound          for r in rxns],
    "ub"        : [r.upper_bound          for r in rxns],
    "c"         : [r.objective_coefficient for r in rxns],
    "rxns"      : [r.id                   for r in rxns],
    "rxnNames"  : [r.name                 for r in rxns],
    "mets"      : [m.id                   for m in mets],
    "metNames"  : [m.name                 for m in mets],
    "subSystems": [r.subsystem or ""      for r in rxns],
    "rev"       : [r.reversibility        for r in rxns],
    "grRules"   : [r.gene_reaction_rule   for r in rxns],
    "genes"     : [g.id                   for g in genes],
    "rxnGeneMat": rxnGeneMat,
    "rxnKEGGID" : [r.annotation.get("kegg.reaction","") for r in rxns],
    "metKEGGID" : [m.annotation.get("kegg.compound","") for m in mets]
}

with open(out, "w") as f:
    json.dump(d, f, indent=2)

