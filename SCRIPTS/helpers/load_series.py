def load_series_sequences(fname):

    series_dict = {}
    with open(fname, 'r') as f:
        for line in f:
            ins = line.strip('\n').split(',')
            series_dict[ins[0]] = [x for x in ins[1:] if x != '']

    return series_dict
