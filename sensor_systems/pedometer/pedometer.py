import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks, butter, filtfilt




######## Functions #########
def prep_data(fn,start=None):
  # Conveniency for me to just use the filename in my repo
  if start == None:
    start = 'https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/pedometer/'
  data = start+fn
  df = pd.read_csv(data)
  df.columns = ['time','x','y','z']
  # Get actual steps (from fn title, leading numbers. eg. 100_pocket.csv -> 100)
  actual_steps = getnum(fn)
  # Bias measured on my own device (changes from device to device)
  bias = {'x':-.1,'y':-.1,'z':.2}
  df.x = df.x - bias['x']
  df.y = df.y - bias['y']
  df.z = df.z - bias['z']
  df['norm'] = np.sqrt((df.x)**2 + (df.y)**2 + (df.z)**2)
  df.time = df.time.apply(pd.to_datetime)
  df.time = (df.time - df.time[0]).dt.total_seconds()
  return actual_steps, df

def plot_raw(df):
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
  plt.show()



def plot_time_domain(df):
  # Extract the z-axis data and sampling frequency
  x = df['norm']
  fs = 6.0

  # Define the passband frequency range
  f_low = 1.5
  f_high = 2.5

  # Apply the bandpass filter
  b, a = butter(2, [2*f_low/fs, 2*f_high/fs], btype='band')
  filtered = filtfilt(b, a, x)

  # Plot the original and filtered signals
  fig, axs = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
  axs[0].plot(df['time'], x)
  axs[0].set_ylabel('Original signal',size=18)
  axs[1].plot(df['time'], filtered)
  axs[1].set_ylabel('Filtered signal',size=18)
  axs[1].set_xlabel('Time (s)',size=18)
  axs[0].set_title(f'Original Signal and Filtered Signal Time Domain \n{actual_steps} steps',size=20)

  plt.show()





def plot_freq_domain(df):
  # Extract the z-axis data and sampling frequency
  z = df['norm']
  fs = 6.0

  # Define the passband frequency range
  fw = 2
  f_low = 1.5
  f_high = 2.5

  # Apply the bandpass filter
  b, a = butter(fw, [2*f_low/fs, 2*f_high/fs], btype='band')
  filtered = filtfilt(b, a, z)

  # Perform Fourier transform on original and filtered signals
  z_fft = np.abs(np.fft.fft(z))
  filtered_fft = np.abs(np.fft.fft(filtered))

  # Create frequency axis
  freq = np.fft.fftfreq(len(z), 1/fs)

  # Plot frequency domain of original and filtered signals
  fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
  axs[0].plot(freq, z_fft)
  axs[0].set_ylabel('Original signal')
  axs[1].plot(freq, filtered_fft)
  axs[1].set_ylabel('Filtered signal')
  axs[1].set_xlabel('Frequency (Hz)')
  axs[0].set_title(f'Original Signal and Filtered Signal Frequency Domain \n{actual_steps} steps',size=20)
  plt.show()


  height_threshold=0.5
  peaks_index, properties = find_peaks(np.abs(z_fft), height=height_threshold)
  

  srt = sorted(list(np.abs(filtered_fft)), reverse=True)[:7]
  corr_freqs = []
  for val in srt:
    idx = np.where(np.abs(filtered_fft) == val)[0][0]
    corr_freqs.append(np.abs(freq[idx]))
  amp_freq = dict(zip(srt,corr_freqs))

  print('Highest Value frequencies:\n\nAmplitude\t\tFrequency\n')
  for a,f in amp_freq.items():
    print(f'{a}\t{f}')



def getnum(fn):
  return int(''.join([n for n in list(fn) if n.isnumeric()]))

def percent_diff(a,b):
  return 100*abs(a - b)/a

def calculate_steps(factor,df):
  d = df['norm']
  df.norm = df.norm - df.norm.mean() # Subtract off influence of gravity. I've found this doesn't really make a difference

  # Apply the bandpass filter
  b, a = butter(2, [2*f_low/fs, 2*f_high/fs], btype='band')
  filtered = filtfilt(b, a, d)
  threshold = factor*max(filtered)

  # Find the times when the filtered signal crosses the threshold in the positive direction
  crossings = np.where(np.diff((filtered > threshold).astype(int)) == 1)[0]

  # Count the number of crossings, which corresponds to the number of steps
  # Num_steps is calculated by taking average of x,y,z steps.
  # ----> doubled because steps are not felt equally for each step
  num_steps = 2*len(crossings)
  return round(num_steps)


# data = 'https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/pedometer/100steps.csv'
# data = '100steps.csv'
# data = '50_m1ix.csv'
# data = '80_hand.csv'
# data = '100_hand.csv'
data = '250_hand.csv'
# data = '120_pocket.csv'
# data = '250_pocket_pause.csv'

# data = '21_juwan_test.csv'
# data = '532.csv'



threshold = .52
fs = 6   # Sampling frequency
f_low = 1.5  # lower bound of walking frequency
f_high = 2.5 # upper bound of walking speed


actual_steps, df = prep_data(data) # Use for data uploaded to GitHub
# actual_steps, df = prep_data(data,'')  #local uploads only, use this


plot_raw(df)
plot_time_domain(df)
plot_freq_domain(df)


calculated_steps = calculate_steps(threshold, df)
pdiff = percent_diff(actual_steps,calculated_steps)

print(f'Actual Steps\t Calculated Steps\n {actual_steps}\t{calculated_steps}\n\tpercent difference: {pdiff}%')
