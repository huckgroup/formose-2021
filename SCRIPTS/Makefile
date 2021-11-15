#############################################################
# This makefile runs all of the scripts required to generate
# analysis results from the primary data. Having generated
# these files, the plots required for the main text are then
# generated.
#############################################################

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

#############################################################
# Figures
#############################################################
# not included in the make ('outside repository')
FIGURE1_PANELS=../PLOTS/Figure1C.svg
FIGURE1_SCRIPTS=../PLOTS/Figure1C.py

FIGURE2_PANEL_A=../PLOTS/Figure2A.svg
FIGURE2_PANEL_B=../PLOTS/Figure2B_Main.svg
FIGURE2_PANEL_C=../PLOTS/Figure2C.svg
FIGURE2_PANEL_D=../PLOTS/Figure2D.svg
FIGURE2_PANEL_E=../PLOTS/Figure2E.svg
FIGURE2_PANEL_F=../PLOTS/Figure2F.svg
FIGURE2_PANELS=$(FIGURE2_PANEL_A) $(FIGURE2_PANEL_B) $(FIGURE2_PANEL_C) \
		$(FIGURE2_PANEL_D) $(FIGURE2_PANEL_E) $(FIGURE2_PANEL_F)

FIGURE3_PANEL_A=../PLOTS/Figure3A.svg
FIGURE3_PANEL_B=../PLOTS/Figure3B_BarChart.svg
FIGURE3_PANEL_C=../PLOTS/Figure3C.svg
FIGURE3_PANEL_D=../PLOTS/Figure3D.svg
FIGURE3_PANELS=$(FIGURE3_PANEL_A) $(FIGURE3_PANEL_B) $(FIGURE3_PANEL_C) $(FIGURE3_PANEL_D)

CONDITIONS_DENDROGRAMS=../PLOTS/[C=O]_Figure2B_inset.svg \
		../PLOTS/[CaCl2]_Figure2B_inset.svg \
		../PLOTS/[NaOH]_Figure2B_inset.svg \
		../PLOTS/[O=C(CO)CO]_Figure2B_inset.svg \
		../PLOTS/residence_time_Figure2B_inset.svg \
		../PLOTS/temperature_Figure2B_inset.svg

RESIDENCE_TIME_SERIES=../PLOTS/Residence_time_series.svg

# not included in the make ('requires chromprocess')
GCMS_HPLC_COMPARISON=../PLOTS/GCMS_HPLC_Comparison.svg
#$(GCMS_HPLC_COMPARISON): DATA_PRESENT/GCMS_HPLC_Comparison.py
#	python DATA_PRESENT/GCMS_HPLC_Comparison.py
#$(GCMS_HPLC_COMPARISON)

SI_FIGURES=$(CONDITIONS_DENDROGRAMS) $(RESIDENCE_TIME_SERIES) $(RESIDENCE_TIME_SERIES)

FIGURES=$(FIGURE2_PANELS) $(FIGURE3_PANELS) $(SI_FIGURES)

all: $(DATA_FILES) $(FIGURES)

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

#############################################################
# Figure Plotting
#############################################################

# Figure 2
$(FIGURE2_PANEL_A): $(AVERAGES) DATA_PRESENT/Figure2A.py
	python DATA_PRESENT/Figure2A.py
$(FIGURE2_PANEL_B): $(AVERAGES) $(AMPLITUDES) $(DENDROGRAM_COORDINATES) DATA_PRESENT/Figure2B_MainPanel.py
	python DATA_PRESENT/Figure2B_MainPanel.py
$(FIGURE2_PANEL_C): $(AVERAGES) DATA_PRESENT/Figure2C.py
	python DATA_PRESENT/Figure2C.py
$(FIGURE2_PANEL_D): $(AVERAGES) $(DENDROGRAM_COORDINATES) DATA_PRESENT/Figure2D.py
	python DATA_PRESENT/Figure2D.py
$(FIGURE2_PANEL_E): $(AVERAGES) DATA_PRESENT/Figure2E.py
	python DATA_PRESENT/Figure2E.py
$(FIGURE2_PANEL_F): $(AVERAGES) $(DENDROGRAM_COORDINATES) DATA_PRESENT/Figure2F.py
	python DATA_PRESENT/Figure2F.py

# Figure 3
$(FIGURE3_PANEL_A): $(PRIMARY_DATA) DATA_PRESENT/Figure3A.py
	python DATA_PRESENT/Figure3A.py
$(FIGURE3_PANEL_B): $(AMPLITUDES) DATA_PRESENT/Figure3_bar_amplitudes.py
	python DATA_PRESENT/Figure3_bar_amplitudes.py
$(FIGURE3_PANEL_C): $(AMPLITUDES) DATA_PRESENT/Figure3C.py
	python DATA_PRESENT/Figure3C.py
$(FIGURE3_PANEL_D): $(REACTION_TYPE_COUNTS) DATA_PRESENT/Figure3D.py
	python DATA_PRESENT/Figure3D.py

# SI
$(CONDITIONS_DENDROGRAMS): $(DENDROGRAM_COORDINATES) DATA_PRESENT/Figure2B_Insets.py
	python DATA_PRESENT/Figure2B_Insets.py
$(RESIDENCE_TIME_SERIES): $(AVERAGES) DATA_PRESENT/SI_Residence_Time_Figure.py
	python DATA_PRESENT/SI_Residence_Time_Figure.py
