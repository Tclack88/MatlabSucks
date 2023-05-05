import pandas as pd
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


def getnum(fn):
  return int(''.join([n for n in list(fn) if n.isnumeric()]))


def percent_diff(a,b):
  return round(100*abs(a - b)/a,1)


def calculate_dist1(factor, pace, df):
  """ estimates steps by filtering input and finding peak frequencies
      Then gets distance by multiplying pace (distance/step) with steps
  """
  #fs = 6
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
  dist = num_steps*pace
  return round(dist,1)

def calculate_dist2(df):
    """estimates steps by integrating acceleration data
    """
    acc = df['norm']
    vel = np.abs(np.cumsum(acc*dt))
    # Numerically integrate the velocity to get distance
    dist = np.cumsum(vel*dt)
    total_dist = np.array(dist)[-1]
    return round(total_dist,1), dist

def linear_est(time,vel):
    dist_arr = vel*time
    return (time,dist_arr)
    

def plot_distances(dat1,dat2,dat3,true_dist):
    t1,d1 = dat1[0], dat1[1]
    t2,d2 = dat2[0], dat2[1]
    t3,d3 = dat3[0], dat3[1]
    plt.plot(t1,d1,t2,d2,t3,d3)
    plt.xlabel('Time (s)')
    plt.ylabel('Distance (m)')
    plt.title(f'True Distance traveled ({true_dist}m) and compared to \ncalculated distance through step counter and integration')
    plt.legend(['True','Step Counter','Integration'])
    plt.show()



data = '532_hand.csv' # Walked around 6th lane of track 438.33 m in 532 steps
threshold = .52 # experimentally determined. Yields best avg. results from all data
threshold = .65# experimentally determined. Yields best avg. results from all data
f_low = 1.5
f_high = 2.5
actual_steps, df = prep_data(data) # Use for data uploaded to GitHub
# actual_steps, df = prep_data(data,'')  # Use for local data only
total_dist = 438.33 # Actual
pace = total_dist/actual_steps # meters/step
tot_time = (df.time.iloc[-1]- df.time[0])
true_vel = total_dist/tot_time
time = df.time
dt = time[1] - time[0]
fs = 1/dt # get sampling freq from collected data

# Estimate the total distance from two methods:
calculated_dist1 = calculate_dist1(threshold, pace, df)
calculated_dist2, dist_arr = calculate_dist2(df)

# construct estimated linear data based on calculations above:
estimated_vel = calculated_dist1/tot_time
dat_real = linear_est(time,true_vel)
dat_est = linear_est(time,estimated_vel)

plot_distances(dat_real, dat_est, (time,dist_arr), total_dist)

# Determine the error
pdiff1 = percent_diff(total_dist,calculated_dist1)
pdiff2 = percent_diff(total_dist,calculated_dist2)

print(f"Actual Distance\t {total_dist}\n\t Calculated (step counter)\t Calculated (integration)\n\t\t\t{calculated_dist1}\t\t{calculated_dist2}\npercent diff.:\t\t{pdiff1}%\t\t{pdiff2}%")
