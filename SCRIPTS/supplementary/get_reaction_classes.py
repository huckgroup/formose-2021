import pickle
import json

fname = '/Users/williamrobinson/documents/nijmegen/packages/rpso/formosepaperanalysis/analysis2/information_sources/FormoseReactionNetwork.pickle'
with open(fname,'rb') as f:
    FormoseNetwork = pickle.load(f)

full_set_reaction_classes = []
for r in FormoseNetwork.NetworkReactions:
    full_set_reaction_classes.append(FormoseNetwork.get_reaction_name(r))


all_reaction_classes = list(set(full_set_reaction_classes))
all_reaction_classes.sort()

reaction_classes = {a:[] for a in all_reaction_classes}

for r in FormoseNetwork.NetworkReactions:
    key = FormoseNetwork.get_reaction_name(r)
    reaction_classes[key].append(r)

for r in reaction_classes:
    reaction_classes[r] = list(set(reaction_classes[r]))
    reaction_classes[r].sort()

json_text = json.dumps(reaction_classes)
with open('information_sources/reaction_class_assignments.json', 'w') as f:
    f.write(json_text)
