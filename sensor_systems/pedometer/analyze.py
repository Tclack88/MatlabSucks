import pandas as pd
import seaborn as sb
from scipy import signal
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np

## some files that exist
# background.csv

# 80_hand.csv
# 100_hand.csv
# 250_hand.csv

# 100steps.csv  (collected in pocket)
# 100steps_pause.csv (collected in pocket) <-------I think the frequency was not 6Hz
# 80_pocket.csv <-------I think the frequency was not 6Hz
# 120_pocket.csv
# 250_pocket_pause.csv


dat = 'https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/pedometer/100steps.csv'

df = pd.read_csv(dat)

df.columns = ['time','x','y','z']

df['norm'] = np.sqrt((df.x)**2 + (df.y)**2 + (df.z)**2)

df.time = df.time.apply(pd.to_datetime)

fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
axs[0].plot(df['time'], df.x)
axs[0].plot(df['time'], df.y)
axs[0].plot(df['time'], df.z)
axs[0].set_ylabel('Individual Accelerations\nm/s${^2}$')

axs[1].plot(df['time'], df.norm)
axs[1].set_ylabel('Average Magnitude Acceleration\nm/s${^2}$')
axs[1].set_xlabel('Time (s)')

axs[0].legend(['x','y','z'],loc='lower right')
axs[0].title.set_text('Acceleration in the Time Domain. Individual and Average')

######### FFT Plot ###########


def plot_freq(df):
  dat = np.array(df.z)
  # dt = (df.time.iloc[-1] - df.time.iloc[0]).total_seconds()
  dt = (df.time.iloc[1] - df.time.iloc[0]).total_seconds()
  dat_f=fft(dat, norm='forward') 

  n = dat.size   # The number of points in the data
  freq = fftfreq(n, d=dt)

  height_threshold=0.5


  peaks_index, properties = find_peaks(np.abs(dat_f), height=height_threshold)

  plt.plot(freq, np.abs(dat_f),'-', freq[peaks_index],properties['peak_heights'],'x')
  plt.xlabel("Frequency (Hz)")
  plt.ylabel("Amplitude")
  plt.title("Discrete Fourier Transform of Pedometer Frequencies\n100 steps in pocket")
  plt.show()
  srt = sorted(list(np.abs(dat_f)), reverse=True)[:7]
  corr_freqs = []
  for val in srt:
    idx = np.where(np.abs(dat_f) == val)[0][0]
    corr_freqs.append(freq[idx])
  amp_freq = dict(zip(srt,corr_freqs))

  print('Highest Value frequencies:\n\nAmplitude\t\tFrequency\n')
  for a,f in amp_freq.items():
    print(f'{a}\t{f}')

plot_freq(df)

####### Filtering Signals with Bandpass #######


# Load the data from CSV file
# data = pd.read_csv('data.csv')

# Extract the z-axis data and sampling frequency
x = df['norm']
fs = 6.0

# Define the passband frequency range
f_low = 1.5
f_high = 2.5

# Apply the bandpass filter
b, a = signal.butter(2, [2*f_low/fs, 2*f_high/fs], btype='band')
filtered = signal.filtfilt(b, a, x)

# Plot the original and filtered signals
fig, axs = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
axs[0].plot(df['time'], x)
axs[0].set_ylabel('Original signal',size=18)
axs[1].plot(df['time'], filtered)
axs[1].set_ylabel('Filtered signal',size=18)
axs[1].set_xlabel('Time (s)',size=18)
axs[0].set_title('Original Signal and Filtered Signal\n100 steps',size=20)

plt.show()


def prep_data(fn,start=None):
  if start == None:
    start = 'https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/pedometer/'
  dat = start+fn
  df = pd.read_csv(dat)
  df.columns = ['time','x','y','z']
  bias = {'x':-.1,'y':-.1,'z':.2}
  df.x = df.x - bias['x']
  df.y = df.y - bias['y']
  df.z = df.z - bias['z']
  df['norm'] = np.sqrt((df.x)**2 + (df.y)**2 + (df.z)**2) ###
  df.time = df.time.apply(pd.to_datetime)
  df.time = (df.time - df.time[0]).dt.total_seconds()
  return df

def getnum(fn):
  return int(''.join([n for n in list(fn) if n.isnumeric()]))

def percent_diff(a,b):
  return 100*abs(a - b)/a

def calculate_steps(factor,df):
  d = df['norm']
  df.norm = df.norm - df.norm.mean() # Subtract off influence of gravity. I've found this doesn't really make a difference

  # Apply the bandpass filter
  b, a = signal.butter(2, [2*f_low/fs, 2*f_high/fs], btype='band')
  filtered = signal.filtfilt(b, a, d)
  threshold = factor*max(filtered)



  # Find the times when the filtered signal crosses the threshold in the positive direction
  crossings = np.where(np.diff((filtered > threshold).astype(int)) == 1)[0]


  # Count the number of crossings, which corresponds to the number of steps
  # Num_steps is calculated by taking average of x,y,z steps.
  # ----> doubled because steps are not felt equally for each step
  num_steps = 2*len(crossings)
  return round(num_steps)

fs = 6.0
delta = .5

freq = 2.0
f_low = freq - delta
f_high = freq + delta

files = ['80_hand.csv','100_hand.csv','250_hand.csv','100steps.csv','120_pocket.csv','250_pocket_pause.csv']

# thresholds = [.3,.4,.5,.6,.7,.8,.9]
thresholds = [.42,.45,.48,.5,.52,.55,.58]



p_diffs = []
for t in thresholds:
  tot_pdiff = 0
  print(f'\n\tusing threshold {t}:')
  for fn in files:
    df = prep_data(fn)
    a = getnum(fn)
    b = calculate_steps(t, df)
    pdiff = percent_diff(a,b)
    tot_pdiff += pdiff
    print(f'{fn} \t\t: {a}, {b}\t {pdiff}')
  avg_pdiff = tot_pdiff/len(files)
  p_diffs.append(avg_pdiff)

plt.plot(thresholds,p_diffs)
plt.xlabel('threshold')
plt.ylabel('percent difference from true value')
plt.title('Average error by selected cutoff threshold value')

fn = "100steps.csv"

fs = 6.0
delta = .5
threshold = .52

freq = 2.0
f_low = freq - delta
f_high = freq + delta

df = prep_data(fn,'')
a = getnum(fn)
b = calculate_steps(threshold, df)
pdiff = percent_diff(a,b)

print(f"""
Actual number of steps taken:\t{a}
Steps calculated by pedometer:\t{b}\n
\tpercent difference:\t{pdiff}%""")

peaks, _ = find_peaks(df.norm, threshold=.12*max(df.norm))
ax, fig = plt.subplots(figsize=[20,7])
sb.lineplot(data=df,x='time',y='norm')
plt.plot(df.time[peaks], df.norm[peaks], "xr"); plt.legend(['signal','detected peaks'])
