# formose-2021
Data and analysis for the "Environmental conditions drive self-organization of reaction pathways in a prebiotic reaction network" (Robinson, Daines, van Duppen, de Jong, Huck, 2021).

Requires [NorthNet](https://github.com/Will-Robin/NorthNet.git), numpy, scipy, networkx and pandas.

## Abstract

The origin of life from an abiotic environment required a gradual process of chemical evolution towards greater molecular complexity. Although elaborate prebiotically-relevant synthetic routes to many of the building blocks of life have been established, it is unclear how functional systems evolved with only the inherent chemical reactivity of molecules as driving forces for organization. Using the formose reaction as a model system, we demonstrate that complex systems of chemical reactions self-organize in response to changes in the environment, strongly narrowing the potential compositional complexity of the reaction outcome. We observe how Breslow’s cycle contributes to the reaction composition by feeding C2 building blocks into the network, alongside reaction pathways dominated by formaldehyde-driven chain growth. The emergence of organized systems of chemical reactions in response to changes in the environment offers a potential mechanism for a chemical evolution process that bridges the gap between prebiotic chemical building blocks and the origin of life.

## Contents
### COMPOUND_INFO

Contains a compound numbering scheme for the manuscript (compound_numbering.txt) and general compound information (compound_properties.csv).

### DATA

##### DATA_REPORTS

Contains folders named according to experiment codes containing data reports of compound concentrations determined from GC-MS and HPLC data.

#### DERIVED_PARAMETERS

Contains the average concentrations (AmplitudeData.csv) and concentration amplitudes (AmplitudeData.csv) determined for each experiment.

### EXPERIMENT_INFO

Contains experimental conditions (Experiment_parameters.csv) for each experiment, and selected experiment code series corresponding to systematic variations in conditions (Series_info.csv).

### FORMOSE_REACTION

Contains a full list of rule-generated formose reactions (FullFormoseReaction.txt) and the reaction classes assigned to each reaction (reaction_class_assignments.json).

### NOTES  
Contains notes on the data analysis.

### PLOTS

A container for plotting data (note: .png files will not be pushed to GitHub). This should be empty.

### REACTION_INFO

Contains a set of reaction SMARTS strings used to generate the formose reaction (reaction_SMARTS_templates.tsv). These reaction SMARTS strings define reaction classes.

### REACTION_LISTS

Lists of reactions for each modulated experiment in as reaction SMILES strings. Files are named according to experiment code.

### SCRIPTS

Python scripts used to analyse and plot the data.
