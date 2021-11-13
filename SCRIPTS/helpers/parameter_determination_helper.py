import numpy as np
from scipy import fftpack
from scipy import interpolate

def fourier_transform(X,Y):
    '''
    Function for determining concentration amplitudes from data.
    '''
    # interpolate axes to get even time steps
    f = interpolate.interp1d(X, Y, kind = "linear")

    # get timestep from average of timepoints
    T = np.diff(X).mean()
    # create new time axis with time step
    new_x = np.arange(X.min(), X.max(), T)
    # cut the new x axis down so it does not
    # exceed the range of the experimental 
    # x axis.
    idx = np.where(new_x < X.max())
    new_x = new_x[idx]

    # create interpolated signal
    new_y = f(new_x)

    # Perform fourier transform
    N = len(new_x)

    yf = fftpack.fft(new_y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    amplitudes = 2.0/N * np.abs(yf[0:N//2])

    return xf, amplitudes
