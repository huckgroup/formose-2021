import numpy as np
from scipy import fftpack
from scipy import interpolate
from scipy.optimize import curve_fit
from NorthNet import info_params

def sine_wave(time, period, amplitude, phase, offset):
    '''
    For generating a sine wave flow profile.

    Parameters
    ----------
    period: float
        The period of the sine wave in seconds.

    amplitude: float
        The amplitude of the sine wave (any unit: think of output).

    phase:
        Phase shift

    offset: float
        The centre of the sine wave.

    Returns
    -------
    wave: 1Darray
        An array of flow rate values.
    '''

    wave = amplitude*np.sin(2*np.pi*time/period + phase) + offset

    return wave

def fit_sine_wave(time,signal, estimate):
    '''
    For fitting a sine wave to a signal trace

    Paramters
    ---------
    time: 1D numpy array
        Time array
    signal: 1D numpy array
        Intensity of signal over time.

    Returns
    -------
    popt: list
        list of fitting parameters: [period, amplitude, phase, offset]

    '''

    if np.average(signal) == 0.0:
        return {"period":1, "amplitude":0, "phase":0, "offset":0}
    else:
        # period, amplitude, phase, offset
        initial_guess = estimate
        lowerbounds   = [0,
                        (np.average(signal)- 2*np.amax(signal)),
                        0,
                        0.0]
        upperbounds   = [20*60,
                         (2*np.amax(signal)+np.average(signal)),
                         np.pi,
                         2*(np.average(signal))]

        boundarr  = [lowerbounds,upperbounds]

        popt, pcov = curve_fit(sine_wave, time, signal, p0=initial_guess,
                                bounds = boundarr)

        if np.average(signal) - popt[1] < 0:
            popt[1] = 0

        return {"period":popt[0], "amplitude":popt[1], "phase":popt[2], "offset":popt[3]}

def autocorrelation(X,Y):
    roll_vals = np.arange(-(len(X))+1,len(X))
    t_step = X[1]-X[0]

    y_i = Y - Y.mean()
    con = np.correlate(y_i,y_i, mode = "full")
    ccor = con/(len(X)*y_i.std()**2)

    return roll_vals, ccor

def fourier_transform(X,Y):
    # interpolate axes to get even time steps
    f = interpolate.interp1d(X, Y, kind = "linear")

    T = np.diff(X).mean()
    new_x = np.arange(X.min(), X.max(), T)
    idx = np.where(new_x < X.max())
    new_x = new_x[idx]
    
    new_y = f(new_x)

    N = len(new_x)

    yf = fftpack.fft(new_y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    amplitudes = 2.0/N * np.abs(yf[0:N//2])

    return xf, amplitudes

def fourier_amplitude(X,Y, cutoffs = [0.0005,0.1]):
    '''
    Get amplitudes from X,Y data via fourier transform.

    Parameters
    ----------
    X: 1D numpy array
        array of x values.
    Y: 1D numpy array
        array of y values.
    cutoff: list of floats
        lower and upper frequencies (in Hz) in which to narrow the area to be
        searched for maximum amplitudes.

    Returns
    -------
    np.nan_to_num(max_amp): 1D numpy array
        An array of maximum amplitude values within the set cutoff frequency
        bounds.
    '''

    xf, amplitudes = fourier_transform(X,Y)
    max_amp = np.amax(amplitudes)

    max_amp = 0.0

    return np.nan_to_num(max_amp)

def t_lag_corr_to_drive(dataset, modulated_input = 'O=C(CO)CO'):
    '''
    parameters
    -----------
    dataset: NorthNet DataReport object
        data
    modulated_input: str
        SMILES string corresponding to the modulated input

    returns
    ------
    t_lag_corr_mat: numpy nd array
        results
    '''
    t_lag_correlations = {}

    for dr in dataset.data_reports:

        report = dataset.data_reports[dr]
        if len(report.series_values) == 0:
            continue

        # some general variables
        time = report.series_values

        modulated_input_key = '{}/ M'.format(modulated_input)
        flow_profile_key = '{}_flow_rate/ Âµl/h'.format(modulated_input)
        flow_time = report.conditions['flow_profile_time/ s']


        input_conc_0 = report.conditions[modulated_input_key]
        modulated_input_flow = report.conditions[flow_profile_key]

        net_flow = np.zeros(len(flow_time))
        for c in report.conditions:
            if 'flow' and not 'time' in c:
                net_flow += report.conditions[c]

        print(report.experiment_code)
        modulated_concentration = (input_conc_0*modulated_input_flow)/net_flow

        roll_vals = np.arange(-(len(time))+1,len(time))
        t_step = time[1]-time[0]

        t_lags = np.zeros(len(roll_vals))

        for r in range(0,len(roll_vals)):
            t_lags[r] = roll_vals[r]*t_step

        # interpolate flow profile to same number of data points
        # as data

        f = interpolate.interp1d(flow_time, modulated_concentration, kind = "linear")
        downsampled_flow = f(time)

        x_i = downsampled_flow - downsampled_flow.mean()

        for e,d in enumerate(report.data):
            x_j = report.data[d] - report.data[d].mean()

            con = np.correlate(x_i,x_j, mode = "full")
            ccor = con/(len(time)*x_i.std()*x_j.std())

            t_lag_correlations[d] = [t_lags,np.nan_to_num(ccor)]


    return t_lag_correlations

def time_lags_to_drive(t_lag_correlations):

    time_lag_dict = {}
    for d in t_lag_correlations:
        import matplotlib.pyplot as plt

        abs_line = np.abs(t_lag_correlations[d][1])
        idx = np.where((abs_line == np.amax(abs_line)))[0]
        time_lag_dict[d]  = t_lag_correlations[d][0][idx]

    return time_lag_dict

def get_amplitudes(dataset):
    '''
    Get amplitudes of the variables in a Dataset objects via Fourier transform.

    Parameters
    ----------
    dataset: Dataset object
        Object containig data.

    Returns
    -------
    amps_dict: dict
        Dictionary of amplitudes of species indexed by their names from Dataset.
    '''

    amps_dict = {}

    for ds in dataset.compounds:
        time, signal = dataset.get_entry(ds)
        amps = fourier_amplitude(time,signal)
        amps_dict[ds] = amps

    return amps_dict
