#############################################################
# This makefile runs all of the scripts required to generate
# analysis results from the primary data. Having generated 
# these files, the plots required for the main text are then 
# generated.
#############################################################

PRIMARY_DATA=$(shell find ../DATA/DATA_REPORTS -type f)
REACTION_SMARTS=../REACTION_INFO/reaction_SMARTS_templates.tsv
FORMOSE_REACTION_LIST=../FORMOSE_REACTION/FullFormoseReaction.txt

AMPLITUDES=../DATA/DERIVED_PARAMETERS/AmplitudeData.csv
AVERAGES=../DATA/DERIVED_PARAMETERS/AverageData.csv

CLUSTERS=../RESOURCES/clusters.txt
LEAVES=../RESOURCES/leaf_list.txt
DENDROGRAM=../PLOTS/dendrogram.svg

REACTION_LISTS=$(shell find ../REACTION_LISTS -type f)
REACTION_TYPE_COUNTS=../RESOURCES/reaction_expression_normalised.csv

all: $(AMPLITUDES) $(AVERAGES) $(CLUSTERS) $(LEAVES) $(DENDROGRAM) $(REACTION_LISTS) $(REACTION_TYPE_COUNTS)

FORMOSE_REACTION_LIST: $(REACTION_SMARTS) DATA_PREP/FormoseReactionGenerationAlgorithm.py
	python DATA_PREP/FormoseReactionGenerationAlgorithm.py

AMPLITUDES: $(PRIMARY_DATA) DATA_PREP/process_all_data_sets.py
	python DATA_PREP/process_all_data_sets.py

AVERAGES: $(PRIMARY_DATA) DATA_PREP/process_all_data_sets.py
	python DATA_PREP/process_all_data_sets.py

CLUSTERS: $(AMPLITUDES) $(AVERAGES) DATA_PREP/Clustering.py
	python DATA_PREP/Clustering.py

LEAVES: $(AMPLITUDES) $(AVERAGES) DATA_PREP/Clustering.py
	python DATA_PREP/Clustering.py

DENDROGRAM: $(AMPLITUDES) $(AVERAGES) DATA_PREP/Dendrogram.py
	python DATA_PREP/Dendrogram.py

REACTION_LISTS: $(AMPLITUDES) $(AVERAGES) $(FORMOSE_REACTION_LIST) DATA_PREP/FindReactionPathways.py
	DATA_PREP/FindReactionPathways.py

REACTION_TYPE_COUNTS: $(REACTION_LISTS) $(CLUSTERS) $(LEAVES) DATA_PREP/ReactionClassAnalysis.py 
	python DATA_PREP/ReactionClassAnalysis.py

