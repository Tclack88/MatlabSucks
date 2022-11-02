from scipy import signal
import matplotlib.pyplot as plt
import sympy as sym
import cmath
import numpy as np
import mpmath as mp
#import wave
#import pyaudio
from scipy.fftpack import fft, fftfreq
from scipy.io import wavfile # get the api
from pylab import *


s = sym.Symbol('s')
t = sym.Symbol('t')


def highlight(word):
    space = ' '*round(len(word)+1)
    center = '#'+space+word+space+'#'
    end = '#'*len(center)
    print(f'\n{end}\n{center}\n{end}')

showall = 0

if showall:
    highlight('Part 1')
    highlight('Q1-2')
    
    num = [7.2, 720, 3888]
    denom = [1, 124.3, 2935.2, 20173.8, 67194, 151308, 0]
    
    sys = signal.TransferFunction(num,denom)
    w,mag,phase = signal.bode(sys)
    plt.figure()
    plt.semilogx(w,mag)
    plt.grid()
    plt.title('mag plot')
    
    plt.figure()
    plt.semilogx(w,phase)
    plt.title('phase plot')
    plt.grid()
    plt.show()
    
    print('low gain at high frequencies')
    print('low (negative) phase at higher frequencies')
    
    highlight('Q3')
    print('should be -106')
    
    Gol = (7.2*s**2 + 720*s + 3888)/(s**6 + 124.3*s**5 + 2935.2*s**4 + 20173.8*s**3 + 67194*s**2 + 151308*s)
    
    lap_sin = sym.laplace_transform(sym.sin(t), t, s, noconds=True)
    #print(lap_sin)
    A = 500
    def lap_f(s):
        return sym.simplify(A*lap_sin*Gol,rational=True)
    
    #f = lap_f(s)
    #print(f)
    
    # Failed to do inverse laplace transform symbolically
    # F = sym.inverse_laplace_transform(f,s,t)
    # print(F)
    # #sym.plot(f(s).
    # import sys; sys.exit(0)
    def f(s):
        # printed from above
        #(because there's an error if I don't write it explicitly)
        return 36000*(s**2 + 100*s + 540)/(s*(s**2 + 1)*(10*s**5 + 1243*s**4 + 29352*s**3 + 201738*s**2 + 671940*s + 1513080))
    t_range = np.linspace(0.001,50,100)
    G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t_range)
    glist = list(G)
    c = np.array(glist, dtype=float)
    plt.plot(t_range,c)
    plt.title('Q3')
    plt.show()
    print('steady state value', c[-1],'\nalso looks like it\'s off by (Q3) -90')
    
    
    highlight('Q4')
    A = 10
    print('with A=500, it\'s 12.885, so with A=10, it\s .258?')
    def lap_f(s):
        return A*lap_sin*Gol
    #print(lap_f(s))
    
    def f(s):
        # printed from above
        #(because there's an error if I don't write it explicitly)
        return 10*(7.2*s**2 + 720*s + 3888)/((s**2 + 1)*(s**6 + 124.3*s**5 + 2935.2*s**4 + 20173.8*s**3 + 67194*s**2 + 151308*s))
    
    
    t_range = np.linspace(0.001,50,100)
    G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t_range)
    glist = list(G)
    c = np.array(glist, dtype=float)
    plt.plot(t_range,c)
    plt.title('Q4')
    plt.show()
    print('steady state value (Q4)', c[-1],'\nalso looks like it\'s off by -90')
    
    
    highlight('Q5')
    A = 1000
    w = .5
    lap_sin = sym.laplace_transform(A*sym.sin(w*t), t, s, noconds=True)
    def lap_f(s):
        return lap_sin*Gol
    #print(lap_f(s))
    
    def f(s):
        # printed from above
        #(because there's an error if I don't write it explicitly)
        return 2000.0*(7.2*s**2 + 720*s + 3888)/((4.0*s**2 + 1)*(s**6 + 124.3*s**5 + 2935.2*s**4 + 20173.8*s**3 + 67194*s**2 + 151308*s))
    
    t_range = np.linspace(0.001,100,500)
    G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t_range)
    glist = list(G)
    c = np.array(glist, dtype=float)
    plt.plot(t_range,c)
    plt.title('Q5')
    plt.show()
    print('steady state value', c[-1],'\nalso looks like it\'s off by -90')
    
    
    highlight('Q6')
    w = 50
    lap_sin = sym.laplace_transform(A*sym.sin(w*t), t, s, noconds=True)
    def lap_f(s):
        return lap_sin*Gol
    #print(lap_f(s))
    
    def f(s):
        # printed from above
        #(because there's an error if I don't write it explicitly)
        return 50000*(7.2*s**2 + 720*s + 3888)/((s**2 + 2500)*(s**6 + 124.3*s**5 + 2935.2*s**4 + 20173.8*s**3 + 67194*s**2 + 151308*s))
    
    t_range = np.linspace(0.001,50,100)
    G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t_range)
    glist = list(G)
    c = np.array(glist, dtype=float)
    plt.plot(t_range,c)
    plt.show()
    plt.title('Q6')
    max_val = max(c)
    print(f'max deviation (max val - ss): {max_val} - {c[-1]} = {abs(max_val - c[-1])}')
    
    highlight('Q7')
    
    highlight('Part 2')
    Gol = (.35*(s+5.7))/(s*(s+5.8)*(s+1.4+3.4j)*(1+1.4-3.4j))
    w = 50
    K = 15

    highlight('Part 3')
    highlight('Q38')
    print('low signal is stil high while the others are attenuated. So this is a LOW PASS FILTER')
    
    highlight('Q39')
    
    
    filename = 'tones1.wav'
    #sound_file = wave.open(filename, 'r')
    
    
    def extract_peak_frequency(data, sampling_rate):
        fft_data = np.fft.fft(data)
        freqs = np.fft.fftfreq(len(data))
    
        peak_coefficient = np.argmax(np.abs(fft_data))
        top_coefficients = np.argsort(np.abs(fft_data))
        peak1 = top_coefficients[-1]
        peak2 = top_coefficients[-2]
        peak3 = top_coefficients[-3]
        #peak_freq = freqs[peak_coefficient]
        peak1_freq = freqs[peak1]
        peak2_freq = freqs[peak2]
        peak3_freq = freqs[peak3]
        p = np.array([peak1_freq,peak2_freq,peak3_freq])
    
        return abs(p * sampling_rate)
    
    sample_rate, data = wavfile.read(filename)
    peak_freqs = extract_peak_frequency(data, sample_rate)
    print(peak_freqs)
    
    # another attempt
    fig, axs = plt.subplots(5,2)
    outputs = ['outputa.wav',  'outputb.wav', 'outputc.wav', 'outputd.wav',  'outpute.wav']
    tones = ['tones1.wav','tones2.wav', 'tones3.wav', 'tones4.wav','tones5.wav']
    
    for i,t in enumerate(tones):
        sample_rate, data = wavfile.read(t)
        peak_freqs = extract_peak_frequency(data, sample_rate)
        #X = fft(data)
        X = abs(fft(data))
        freqs = fftfreq(len(data))*sample_rate
        axs[i,0].set_xlim([-2000,2000])
        axs[i,0].plot(freqs,X)
        axs[i,0].set_title(t)
        p = extract_peak_frequency(data,sample_rate)
        print(t,'\n\t',p)
    
    print()
    for i,o in enumerate(outputs):
        sample_rate, data = wavfile.read(o)
        peak_freqs = extract_peak_frequency(data, sample_rate)
        X = abs(fft(data))
        freqs = fftfreq(len(data))*sample_rate
        axs[i,1].set_xlim([-2000,2000])
        axs[i,1].plot(freqs,X)
        axs[i,1].set_title(o)
        p = extract_peak_frequency(data,sample_rate)
        print(o,'\n\t',p)
    
    
    plt.show()
    

print('1 e \n 2 c\n3 a\n4 b\n5 d')


highlight('Q40')

import scipy.io
import control.matlab as cnt
freq = scipy.io.loadmat('freq.mat')
tf_data = scipy.io.loadmat('tf_data-1.mat')
freq = freq['f']
tf_data = tf_data['tf_data']
tf_data[5]
fig, axs = plt.subplots(2)
for i in [0,1]:#range(len(tf_data)):
    print(i)
    y = tf_data[i].reshape(1,25001)
    a = cnt.frd(freq, y)
    print(dir(a))
    axs[i].plot(freq,y)
#
#plt.show()
#import scipy
#import control.matlab as cnt
#freq = scipy.io.loadmat('freq.mat')
#tf_data = scipy.io.loadmat('tf_data-1.mat')
#print(tf_data.keys())
#a = cnt.frd(freq, tf_data[0])
#print(type(a))

import sys; sys.exit(0)



#plt.plot(freq,tf_data)



sys = signal.TransferFunction(num,denom)
w,mag,phase = signal.bode(sys)
plt.figure()
plt.semilogx(w,mag)
plt.grid()
plt.title('mag plot')

plt.figure()
plt.semilogx(w,phase)
plt.title('phase plot')
plt.grid()
plt.show()

#### ANOTHER attempt
# highlight('newest attempt')
# sample_rate, data = wavfile.read(filename)
# 
# def get_peak(data):
#     DatArr = data.T
#     y = DatArr
#     t = np.linspace(0, len(DatArr)-1, len(DatArr))
#     
#     y = signal.detrend(y,type='linear')
#     f1, ax1 = plt.subplots()
#     ax1.plot(t,y)
#     f1.show()
#     
#     ft = np.fft.fft(y,n=16*len(DatArr))
#     ftnorm = abs(ft)
#     ps = ftnorm**2
#     
#     xvals = np.fft.fftfreq(len(ps), 1.0/len(DatArr))
#     f2, ax2 = plt.subplots()
#     ax2.plot(xvals,ps)
#     ax2.set_xlim(0,int(np.max(xvals)))
#     plt.title("Mikayla Smells")
#     f2.show()
#     
#     peak = np.argmax(ps)
#     #print("Peak Frequency:",xvals[peak])
#     return peak
# 
# 
# 
# for i,t in enumerate(tones):
#     sample_rate, data = wavfile.read(t)
#     peak = get_peak(data)
#     print(t,'\n\t',p)
# 
# print()
# 
# for i,o in enumerate(outputs):
#     sample_rate, data = wavfile.read(o)
#     peak = get_peak(data)
#     print(o,'\n\t',p)
# 
# 
# 
# 
# import sys; sys.exit(0)
# 
# def f(filename):
#     fs, data = wavfile.read(filename) # load the data
#     a = data.T[0] # this is a two channel soundtrack, I get the first track
#     a = data.T[1] # this is a two channel soundtrack, I get the first track
#     print(data)
#     print(a)
#     print(type(a))
#     b=[(ele/2**8)*2-1 for ele in data] # this is 8-bit track, b is now normalized on [-1,1)
#     c = fft(b) # create a list of complex number
#     d = len(c)/2  # you only need half of the fft list
#     plt.plot(abs(c[:(d-1)]),'r')
#     plt.show()
#     #savefig(filename+'.png',bbox_inches='tight')
# 
# f(filename)

# N = SAMPLE_RATE * DURATION
# 
# yf = fft(normalized_tone)
# xf = fftfreq(N, 1 / SAMPLE_RATE)
# 
# plt.plot(xf, np.abs(yf))
# plt.show()
# 
# 
# file_length = sound_file.getnframes()
# 
# sound = np.zeros(file_length)
# mean_square = []
# sound_square = np.zeros(file_length)
# for i in range(file_length):
#     data = sound_file.readframes(1)
#     data = struct.unpack("<h", data)
#     sound[i] = int(data[0])
#     
# sound = np.divide(sound, float(2**15))  # Normalize data in range -1 to 1
# 
# 
