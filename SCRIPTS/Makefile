#############################################################
# This makefile runs all of the scripts required to generate
# analysis results from the primary data.
##############################################################

#############################################################
# Data Files
#############################################################

PRIMARY_DATA=$(shell find ../DATA/DATA_REPORTS -type f)
REACTION_SMARTS=../REACTION_INFO/reaction_SMARTS_templates.tsv
FORMOSE_REACTION_LIST=../FORMOSE_REACTION/FullFormoseReaction.txt

AMPLITUDES=../DATA/DERIVED_PARAMETERS/AmplitudeData.csv
AVERAGES=../DATA/DERIVED_PARAMETERS/AverageData.csv

CLUSTERS=../RESOURCES/clusters.txt
LEAVES=../RESOURCES/leaf_list.txt
DENDROGRAM_COORDINATES=../RESOURCES/dendrogram_coordinates.csv
DENDROGRAM=../PLOTS/dendrogram.svg

REACTION_LISTS=$(shell find ../REACTION_LISTS -type f)
REACTION_TYPE_COUNTS=../RESOURCES/reaction_expression_normalised.csv

DATA_FILES=$(AMPLITUDES) $(AVERAGES) $(CLUSTERS) $(LEAVES) \
					$(DENDROGRAM) $(REACTION_LISTS) $(REACTION_TYPE_COUNTS)

all: $(DATA_FILES)

#############################################################
# Data Plotting
#############################################################

$(FORMOSE_REACTION_LIST): $(REACTION_SMARTS) DATA_PREP/FormoseReactionGenerationAlgorithm.py
	python DATA_PREP/FormoseReactionGenerationAlgorithm.py
$(AMPLITUDES): $(PRIMARY_DATA) DATA_PREP/process_all_data_sets.py
	python DATA_PREP/process_all_data_sets.py
$(AVERAGES): $(PRIMARY_DATA) DATA_PREP/process_all_data_sets.py
	python DATA_PREP/process_all_data_sets.py

#############################################################
# These outputs are 'non-deterministic', so they will run every
# time make is called, thereby triggering all the files which
# depend on then downstream.
$(CLUSTERS): $(AMPLITUDES) $(AVERAGES) DATA_PREP/Clustering.py
	python DATA_PREP/Clustering.py
$(LEAVES): $(AMPLITUDES) $(AVERAGES) DATA_PREP/Clustering.py
	python DATA_PREP/Clustering.py
$(DENDROGRAM_COORDINATES): $(AMPLITUDES) $(AVERAGES) DATA_PREP/Clustering.py
	python DATA_PREP/Clustering.py
#############################################################

$(DENDROGRAM): $(AMPLITUDES) $(AVERAGES) DATA_PREP/Dendrogram.py
	python DATA_PREP/Dendrogram.py
$(REACTION_LISTS): $(AMPLITUDES) $(AVERAGES) $(FORMOSE_REACTION_LIST) DATA_PREP/FindReactionPathways.py
	python DATA_PREP/FindReactionPathways.py
$(REACTION_TYPE_COUNTS): $(REACTION_LISTS) DATA_PREP/ReactionClassAnalysis.py
	python DATA_PREP/ReactionClassAnalysis.py

