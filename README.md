# formose-2021
Data and analysis for the "Environmental conditions drive self-organization of reaction pathways in a prebiotic reaction network" (Robinson, Daines, van Duppen, de Jong, Huck, 2021).

Requires [NorthNet](https://github.com/Will-Robin/NorthNet.git), numpy, scipy, networkx and pandas.

## Abstract

The origin of life from an abiotic environment required a gradual process of chemical evolution towards greater molecular complexity. Although elaborate prebiotically-relevant synthetic routes to many of the building blocks of life have been established, it is unclear how functional systems evolved with only the inherent chemical reactivity of molecules as driving forces for organisation. Using the formose reaction as a model system, we demonstrate that complex systems of chemical reactions self-organise in response to changes in the environment, strongly narrowing the potential compositional complexity of the reaction outcome. We observe how Breslow’s cycle contributes to the reaction composition by feeding C2 building blocks into the network, alongside reaction pathways dominated by formaldehyde-driven chain growth. The emergence of organised systems of chemical reactions in response to changes in the environment offers a potential mechanism for a chemical evolution process that bridges the gap between prebiotic chemical building blocks and the origin of life.

## Contents
### COMPOUND_INFO

Contains a compound numbering scheme for the manuscript (compound_numbering.txt) and general compound information (compound_properties.csv).

### DATA

##### DATA_REPORTS

Contains folders named according to experiment codes containing data reports of compound concentrations determined from GC-MS and HPLC data.

#### DERIVED_PARAMETERS

Contains the average concentrations (AverageData.csv) and concentration amplitudes (AmplitudeData.csv) determined for each experiment.

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

## Installation

The programs contained in this repository require Python 3.9.2 and the rdkit, numpy, scipy, networkx and pandas matplotlib and sklearn libraries (see `environment.yaml`)

This software should work on all systems capable of installing the dependencies described above. It has been run successfully on MacOS (10.15) and Windows (Windows 10) machines.

A typical installation time should take less than 1 hour. The software can be installed as follows:

### 1. Clone

Clone the repository from GitHub. If you're having trouble with cloning, just download the zip file and store it on your computer wherever you would like.

### 2. Create a virtual environment

[Anaconda/Miniconda](https://www.anaconda.com/products/individual-b#Downloads, 'Anaconda') (`conda`) or `pip` can be used to install the software dependencies.

First, create a virtual environment:

#### Using conda:

Create a virtual environment with conda:

`conda create --name name-env`

Activate the virtual environment:

Mac:

`conda activate name-env`

Windows:

`activate name_of_my_env`

**Go to Install dependencies.**

#### Using pip:

Create a virtual environment:

`pip install virtualenv`

`virtualenv name-env`

Activate the virtual environment:

Mac:

`source name-env/bin/activate`

Windows:

`name-env\Scripts\activate`

**Go to Install dependencies.**

### Install dependencies

Use `pip` or `conda` to install the following dependencies.

e.g. using conda:
- `conda install -c anaconda scipy`
- `pip install networkx`
- `conda install matplotlib`
- `conda install pandas`
- `conda install -c rdkit`

### Install NorthNet
In command line/terminal, navigate to the folder containing the NorthNet code, then type:

conda:
  - `conda develop NorthNet`

  (you may have to run `conda install conda-build` first)

pip:
  - Install from the repository root using pip: `pip install .`,
  - Or in editable mode (so edits are immediately reflected): `pip install -e .`
