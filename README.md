# Environmental conditions drive self-organization of reaction pathways in a prebiotic reaction network

Data and analysis for the "Environmental conditions drive self-organization of reaction pathways in a prebiotic reaction network" (Robinson, Daines, van Duppen, de Jong, Huck, 2021).

## Abstract

The evolution of life from the prebiotic environment required a gradual process of chemical evolution towards greater molecular complexity. Elaborate prebiotically-relevant synthetic routes to the building blocks of life have been established. However, it is still unclear how functional chemical systems evolved with direction using only the interaction between inherent molecular chemical reactivity and the abiotic environment. Here, we demonstrate how complex systems of chemical reactions exhibit well-defined self-organisation in response to varying environmental conditions. This self-organisation allows the compositional complexity of the reaction products to be controlled as a function of factors such as feedstock and catalyst availability. We observe how Breslow's cycle contributes to the reaction composition by feeding C<sub>2</sub> building blocks into the network, alongside reaction pathways dominated by formaldehyde-driven chain growth. The emergence of organised systems of chemical reactions in response to changes in the environment offers a potential mechanism for a chemical evolution process that bridges the gap between prebiotic chemical building blocks and the origin of life.

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

A container for plotting data (note: .png files are ignored by the git repository). This should be empty.

### REACTION_INFO

Contains a set of reaction SMARTS strings used to generate the formose reaction (reaction_SMARTS_templates.tsv). These reaction SMARTS strings define reaction classes.

### REACTION_LISTS

Lists of reactions for each modulated experiment in as reaction SMILES strings. Files are named according to experiment code.

### SCRIPTS

Python scripts used to analyse and plot the data. It is possible to run each script individually, priovided that the correct data files are present in the repository. Reading and writing of files is designed to work within the directory structure of this repository. The primary source of data with respect to the analysis in this repository are the data files in `DATA/DATA_REPORTS` and the reaction rules outlined in `REACTION_INFO/reaction_smarts_templates.csv`. All other data files are derived from these sources using the scripts in `SCRIPTS/DATA_PREP`. For more details on the dependencies of the data, please see the `Makefile`.

The scripts in `SCRIPTS/DATA_PRESENT` are for data visualisation purposes (for example, for inclusion as figures in the main text and supplementary information). They convert the various data files into graphical formats.

## Installation

The programs contained in this repository require Python 3.9.2 and the rdkit, numpy, scipy, networkx, pandas matplotlib and sklearn libraries (see `environment.yml`).

This software should work on all systems capable of installing the dependencies described above. It has been run successfully on MacOS (10.15) and Windows (Windows 10) machines.

A typical installation time should take less than 1 hour. The software can be installed as follows:

### 1. Clone

Clone the repository from GitHub. If you're having trouble with cloning, just download the zip file and store it on your computer wherever you would like.

### 2. Create a virtual environment

[Anaconda/Miniconda](https://www.anaconda.com/products/individual-b#Downloads, 'Anaconda') (`conda`) or `pip` can be used to install the software dependencies.

First, create a virtual environment:

#### Using conda:

Create a virtual environment with conda:

`conda create --name name-env python=3.9.2`

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
- `conda install -c conda-forge rdkit rdkit`
- `conda install -c anaconda scikit-learn `
- `pip install graphviz`
- `conda install --channel conda-forge pygraphviz`

### Install NorthNet

First, download or clone a copy of the NorthNet code (please use [the v0.1 release](https://github.com/Will-Robin/NorthNet/releases/tag/v0.1)). In the command line/terminal, navigate to the repository directory and run:

conda:

  - `conda develop NorthNet`

  (you may have to run `conda install conda-build` first)

pip:

  - Install from the repository root using pip: `pip install .`,
  - Or in editable mode (so edits are immediately reflected): `pip install -e .`

## Expected run times

Each of the programs should not take more than a few minutes to run on a modern laptop.
