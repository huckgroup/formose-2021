from pathlib import Path

# get the path of the project folder
# (one directory above that of this script)
path = Path(__file__)
script_dir = path.parent
parent = script_dir.parents[1]

# to be implemented
# import pandas as pd
# compound_info = pd.read_csv(parent/'COMPOUND_INFO/compound_properties.csv', index_col = 0)

# load
info_container = []
with open(parent/'COMPOUND_INFO/compound_properties.csv', 'r') as f:
    for c, line in enumerate(f):
        if c == 0:
            header = line.strip('\n').split(',')
        else:
            spl = line.strip('\n').split(',')
            info_container.append(spl)

info_container = [list(i) for i in zip(*info_container)]

props_dict = {}
for n,i in zip(header, info_container):
    props_dict[n] = i

colour_assignments = {k:v for k,v in
                              zip(props_dict['@ SMILES'], props_dict['colour'])}

for c,p in enumerate(props_dict['@@ SMILES']):
    if p in colour_assignments:
        pass
    else:
        colour_assignments[p] = props_dict['colour'][c]


for c,p in enumerate(props_dict['Other_names']):
    for s in p.split(';'):
        if s == '':
            pass
        elif s in colour_assignments:
            pass
        else:
            colour_assignments[s] = props_dict['colour'][c]

molecular_masses = {k:float(v) for k,v in
                           zip(props_dict['@ SMILES'], props_dict['Mr_gmol-1'])}
canonical_SMILES = {k:v for k,v in
                       zip(props_dict['compound_name'], props_dict['@ SMILES'])}

for a,b in zip(props_dict['Other_names'], props_dict['@ SMILES']):
    for s in a.split(';'):
        if s != '':
            canonical_SMILES[s] = b

smiles_to_names = {}
for c,v in enumerate(props_dict['compound_name']):
    spl_name = v.split(' ')[0]
    smiles_to_names[props_dict['@ SMILES'][c]] = spl_name

class_assignments =  {k:v for k,v in
                               zip(props_dict['@ SMILES'], props_dict['Class'])}

for sm,cls in zip(props_dict['@@ SMILES'], props_dict['Class']):
    class_assignments[sm] = cls

reaction_SMARTS = {}
reaction_class_colours = {}
reaction_class_short_names = {}
reaction_class_names = {}
with open(parent/'REACTION_INFO/reaction_SMARTS_templates.tsv', 'r') as f:
    for c,line in enumerate(f):
        if c==0:
            pass
        else:
            ins = line.strip('\n').split('\t')
            reaction_SMARTS[ins[0]] = ins[3]
            reaction_class_colours[ins[0]] = ins[4]
            reaction_class_short_names[ins[0]] = ins[5]
            reaction_class_names[ins[0]] = ins[6]

with open(parent/'COMPOUND_INFO/compound_numbering.txt', 'r') as f:
    lines = f.readlines()

lines = [l.strip('\n') for l in lines]
compound_numbering = {}
for l in lines:
    SMILES_num = l.split(',')
    compound_numbering[SMILES_num[0]]  = SMILES_num[1]


reaction_colours = {}
with open(parent/"REACTION_INFO/reaction_colour_assignments.csv", 'r') as f:
    lines = f.readlines()

lines = [l.strip('\n').split(',') for l in lines][1:]

reaction_colours = {l[0]:l[1] for l in lines}
