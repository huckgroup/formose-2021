# formose-2021
Data and analysis for the formose paper (Robinson, Daines, van Duppen, de Jong, Huck, 2021).

Requires installation of NorthNet (ask W.E.R.), numpy, scipy

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

A container for plotting data (note: .png files will not be pushed to GitHub).

### REACTION_INFO

Contains a set of reaction SMARTS strings used to generate the formose reaction (reaction_SMARTS_templates.tsv). These reaction SMARTS strings define reaction classes.

### REACTION_LISTS

Lists of reactions for each modulated experiment in as reaction SMILES strings. Files are named according to experiment code.

### SCRIPTS

Python scripts used to analyse and plot the data.
