#############################################################
# This makefile runs all of the scripts required to generate
# analysis results from the primary data.
# NOTE: only tested on MacOS.
##############################################################

##########################
# Data sources definitions 
##########################
## Data
PRIMARY_DATA=$(shell find ../DATA/DATA_REPORTS -type f)
## Reaction rules
REACTION_SMARTS=../REACTION_INFO/reaction_SMARTS_templates.tsv

#################################
# Generated data files definitions
#################################
# Compound data
AMPLITUDES=../DATA/DERIVED_PARAMETERS/AmplitudeData.csv
AVERAGES=../DATA/DERIVED_PARAMETERS/AverageData.csv
CLUSTERS=../RESOURCES/clusters.txt
LEAVES=../RESOURCES/leaf_list.txt
DENDROGRAM_COORDINATES=../RESOURCES/dendrogram_coordinates.csv
DENDROGRAM=../PLOTS/dendrogram.svg
# Reaction data
FORMOSE_REACTION_LIST=../FORMOSE_REACTION/FullFormoseReaction.txt
# The line below a bit clumsy, but it refers to the last reaction list 
# generated. If that file is kept up to date, the other reaction 
# lists should also be kept up to date
REACTION_LISTS=../REACTION_LISTS/FRN093A_reaction_list.txt
#REACTION_LISTS=$(shell find ../REACTION_LISTS -type f)
REACTION_TYPE_COUNTS=../RESOURCES/reaction_expression_normalised.csv

DATA_FILES=$(AMPLITUDES) $(AVERAGES) $(CLUSTERS) $(LEAVES) \
			$(DENDROGRAM) $(REACTION_LISTS) $(REACTION_TYPE_COUNTS)

all: $(DATA_FILES)

#############################################################
# Data processing rules
#############################################################

$(FORMOSE_REACTION_LIST): $(REACTION_SMARTS)
	python DATA_PREP/FormoseReactionGenerationAlgorithm.py
$(AMPLITUDES) $(AVERAGES): $(PRIMARY_DATA)
	python DATA_PREP/process_all_data_sets.py

$(CLUSTERS) $(LEAVES) $(DENDROGRAM_COORDINATES): $(AMPLITUDES) $(AVERAGES)
	python DATA_PREP/Clustering.py
$(DENDROGRAM): $(AMPLITUDES) $(AVERAGES)
	python DATA_PREP/Dendrogram.py

$(REACTION_LISTS): $(AMPLITUDES) $(AVERAGES) $(FORMOSE_REACTION_LIST)
	python DATA_PREP/FindReactionPathways.py
$(REACTION_TYPE_COUNTS): $(REACTION_LISTS)
	python DATA_PREP/ReactionClassAnalysis.py

