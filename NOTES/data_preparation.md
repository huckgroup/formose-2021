Data preparation
---

---
In DATA_PREP
### 1. Copy data reports into project directory
- copy_GCMS_files.py
- copy_HPLC_files.py

### 2. Clean data sets

- clean_up_data_reports.py

Data sets are copied, cleaned and deposited in a new directory.

### 3. Extract parameters from cleaned data sets
- process_all_data_sets.py

The averages, amplitudes, time lags and errors are extracted from each data set and deposited in their own .csv files (see information_sources).

### 4. Generate reaction networks
- 4_GenerationAlgorithm.py
Uses reaction rules from the information folder in NorthNet

(deposited in information_sources)

### 5. Find reaction pathways in all relevant data sets
- Find_reaction_pathways_2.py

Uses the network generated in step 4 to find reaction pathways within data from which amplitude info is present.

(see also: reaction_pathway_searching.md)

---

---
In DATA_PRESENT
### 6. Data presentation
Use scripts to generate figures from the data.
- bar_chart_panel.py
- dendrogram_panel.py
- plot_series_2.py
- plot_series_grid.py

---
