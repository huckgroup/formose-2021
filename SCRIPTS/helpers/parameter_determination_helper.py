import numpy as np
from scipy import fftpack
from scipy import interpolate
from scipy.optimize import curve_fit
from NorthNet import info_params

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
