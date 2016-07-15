#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq, fftshift

if __name__ == "__main__":
    """
    """
    print("\n ### Testing FFT with a sinusoidal curve ###\n")
    N= 2000 # Sample points
    dt = 0.01 # t-step
   #T = N*dt # Period
    print(" Sample points N:", N)
    print(" Sample spacing dt:", dt)
   #print(" Period T = N*dt:", T)
   #print(" Frequency F=1/T:", 1/T)
   #df = 1/T # Frequency step
   #dw = 2*np.pi/T # Pulse step
    t = np.linspace(-N*dt, N*dt, N)
    #print(t)
    # Defining a functin in time
    f = lambda x: np.sin(5.0 * 2.0*np.pi*x) + 0.5*np.sin(8.0 * 2.0*np.pi*x)
    yt = f(t)
    f = np.random.rand(len(t))*np.sin(2*np.pi*t*10)
    yt = f
    # Printing the function in time domain
    plt.figure(); plt.grid()
   #plt.xlim(0,N*dt)
    plt.plot(t,yt)
    # Performing the Fourier transform
    fft = fft(yt)
    xf = fftfreq(N,dt)
    #w = np.fft.fftfreq(N)*N*dw
    # Shifting the frequencies for visualization
    sxf = fftshift(xf)
    sfft = fftshift(fft)
   #sw = np.fft.fftshift(w)
    #plt.figure()
    #plt.plot(f,fft)
    #plt.figure()
    #plt.plot(w,fft)
   #plt.figure()
   #plt.plot(t,xf,t,sxf)
   #plt.plot(xf,abs(fft))
    # Plotting the frequency spectrum
    plt.figure(); plt.grid()
    plt.xlim(0,sxf[-1])
    plt.plot(sxf,abs(sfft)/N)
    plt.show()
