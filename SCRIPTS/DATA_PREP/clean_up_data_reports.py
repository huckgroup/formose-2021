import os
import numpy as np
from pathlib import Path
from NorthNet import Classes
from NorthNet import info_params
from NorthNet.file_loads import data_loads, info_loads

base_directory = Path('/Users/williamrobinson/Documents/Nijmegen')
# base_directory = Path(r'C:\Users\willi\Documents')
report_directory = base_directory/'safestore_DEP'
cleaned_output_directory = report_directory/'Cleaned_data_reports'
exp_info = info_loads.import_Experiment_information(report_directory/"Experiment_parameters.csv")

point_removals = {}
with open(report_directory/"Data_point_removals.csv", 'r') as f:
    for c,line in enumerate(f):
        if c == 0:
            pass
        else:
            ins = [x for x in line.strip('\n').split(',')]

            point_removals[ins[0]] = {'HPLC':[int(x) for x in ins[1].split(';')if x!=''],
                                      'GCMS':[int(x) for x in ins[2].split(';')if x!='']}

data_reports = []
for e in exp_info:
    for file in os.listdir(exp_info[e].path):
        fname = '{}/{}'.format(exp_info[e].path,file)
        if file.endswith('.csv'):
            data_reports.append(Classes.DataReport(file = fname))

# clean up input concentration keys
for d in data_reports:
    condition_replacements = {}
    for c in d.conditions:
        if ' M' in c:
            standardised_key = c.split('_')[0].split('/')[0]
            clef = info_params.canonical_SMILES[standardised_key]
            condition_replacements[c] = clef + '/ M'

    for c in condition_replacements:
        d.conditions[condition_replacements[c]] = d.conditions[c]
        del d.conditions[c]

# clean up flow profile keys
for d in data_reports:
    condition_replacements = {}
    for c in d.conditions:
        if 'flow' in c and not 'time' in c:
            standardised_key = c.split(' ')
            clef = info_params.canonical_SMILES[standardised_key[0].split('_')[0]]
            condition_replacements[c] = clef + '_flow_rate/ ' + "µl/h"
    for c in condition_replacements:
        d.conditions[condition_replacements[c]] = d.conditions[c]
        del d.conditions[c]

def check_zero_complimentarity(arr1, arr2):
    '''
    Assumes arr1, arr2 are 1D arrays of equal length
    '''
    zeros_1 = arr1 > 0.0
    zeros_2 = arr2 > 0.0
    summed = np.sum(zeros_1 + zeros_2)

    if summed == len(arr1):
        return True
    else:
        return False

# clean up repeated compound entries
# clean up columns containing zeroes
# clean up column headers
get_rt= lambda x: float(x[x.find("(")+1:x.find(")")])
for d in data_reports:
    data_header_replacements = {}
    compound_list = []
    for dep in d.data:
        standardised_key = dep.split('/')[0].split(' ')[0]
        compound_list.append(standardised_key)
        data_header_replacements[dep] = standardised_key + '/ M'

    # check for repeats
    for c in compound_list:
        temp_key_list = []
        for dep in d.data:
            if c in dep:
                temp_key_list.append(dep)

        temp_arr = np.zeros((len(temp_key_list),len(d.series_values)))
        for c,k in enumerate(temp_key_list):
            temp_arr[c] = d.data[k]

        # check entries pairwise for complimentary and
        # if they are not complimentary in terms of zero entries,
        # delete the entry with a lower summed concentration
        del_list = []
        combine_list = []
        for x in range(0,len(temp_arr)):
            for y in range(0,len(temp_arr)):
                if x == y:
                    pass
                elif check_zero_complimentarity(temp_arr[x], temp_arr[y]):
                    combine_list.append((temp_key_list[x],temp_key_list[y]))
                elif np.sum(temp_arr[x]) > np.sum(temp_arr[y]):
                    del_list.append(temp_key_list[y])
                elif np.amax(temp_arr[y]) == 0.0:
                    del_list.append(temp_key_list[y])
                else:
                    pass

        # removing entries containing only zeroes
        for dep in d.data:
            if np.sum(d.data[dep]) == 0.0:
                del_list.append(dep)

        for c in combine_list:
            if c[0] in del_list:
                pass
            else:
                d.data[c[0]] = d.data[c[0]] + d.data[c[1]]
                d.errors[c[0]] = d.errors[c[0]] + d.errors[c[1]]
                del_list.append(c[1])

        if len(del_list) == 0:
            pass
        else:
            del_list = list(set(del_list))
            for x in del_list:
                del d.data[x]
                del d.errors[x]
                del data_header_replacements[x]

    for dep in data_header_replacements:
        d.data[data_header_replacements[dep]] = d.data[dep]
        d.errors[data_header_replacements[dep]] = d.errors[dep]
        del d.data[dep]
        del d.errors[dep]

    # correct xylitol SMILES
    xylitol_key = False
    glyc_dimer_key = False
    for dep in d.data:
        if 'OC[C@H](O)C(O)[C@H](O)CO' in dep:
            xylitol_key = dep
            new_key = dep.replace('OC[C@H](O)C(O)[C@H](O)CO', 'OC[C@H](O)[C@@H](O)[C@H](O)CO')
        if 'CC1COC(C)CO1' in dep:
            glyc_dimer_key = dep
            new_key2 = dep.replace('CC1COC(C)CO1', 'OC1COC(O)CO1')

    if xylitol_key:
        d.data[new_key] = d.data[xylitol_key]
        d.errors[new_key] = d.errors[xylitol_key]

        del d.data[xylitol_key]
        del d.errors[xylitol_key]

    if glyc_dimer_key:
        d.data[new_key2] = d.data[glyc_dimer_key]
        d.errors[new_key2] = d.errors[glyc_dimer_key]

        del d.data[glyc_dimer_key]
        del d.errors[glyc_dimer_key]

# remove internal standard entries
for d in data_reports:
    del_list = []
    for dep in d.data:
        if 'CCCCCCCCCCCCCC' in dep:
            del_list.append(dep)
    for dep in del_list:
        del d.data[dep]
        del d.errors[dep]


# remove point removals from the data
for d in data_reports:
    exp_code = d.experiment_code
    print(exp_code)
    if exp_code in point_removals:

        if 'Chromatography_method' in d.analysis_details:
            c_type = d.analysis_details['Chromatography_method'][0]
        elif 'GCMS' in d.filename:
            c_type = 'GCMS'
        elif 'HPLC' in d.filename:
            c_type = 'HPLC'

        remove_idx = point_removals[exp_code][c_type]
        print(len(d.series_values))
        d.series_values = np.delete(d.series_values, remove_idx)
        for dep in d.data:
            d.data[dep] = np.delete(d.data[dep],remove_idx)
            d.errors[dep] = np.delete(d.errors[dep],remove_idx)
        print(len(d.series_values))
        print()

for d in data_reports:
    store_path = cleaned_output_directory/d.experiment_code
    os.makedirs(store_path, exist_ok = True)
    exp_code = d.experiment_code
    if 'Chromatography_method' in d.analysis_details:
        data_type = d.analysis_details['Chromatography_method'][0]
    elif 'HPLC' in d.filename:
        data_type = 'HPLC'
    elif 'GCMS' in d.filename:
        data_type = 'GCMS'
    else:
        data_type = 'not_specified'

    file_name = '{}_{}_concentration_report'.format(exp_code, data_type)
    d.write_to_file(filename = file_name, path = store_path)
