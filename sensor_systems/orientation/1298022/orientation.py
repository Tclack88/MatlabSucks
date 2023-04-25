import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt
from filterpy.kalman import KalmanFilter
#from step_counter_functions import * # bring in functions from pedometer


# This collected data is 90 degrees CCW, 90 degrees CW. 90 degrees CW and another 90 CW (for a total of 180).

##### Step counter functions brought in from pedometer ######
def prep_data(fn,start=None):
  # Conveniency for me to just use the filename in my repo
  if start == None:
    start = 'https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/pedometer/'
  data = start+fn
  df = pd.read_csv(data)
  df.columns = ['time','x','y','z']
  # Get actual steps (from fn title, leading numbers. eg. 100_pocket.csv -> 100)
  # Bias measured on my own device (changes from device to device)
  bias = {'x':-.1,'y':-.1,'z':.2}
  df.x = df.x - bias['x']
  df.y = df.y - bias['y']
  df.z = df.z - bias['z']
  df['norm'] = np.sqrt((df.x)**2 + (df.y)**2 + (df.z)**2)
  df.time = df.time.apply(pd.to_datetime)
  df.time = (df.time - df.time[0]).dt.total_seconds()
  return df


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



# load accelerometer and magnetometer data
acc_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/orientation/physphox/Accelerometer.csv"
gyro_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/orientation/physphox/Gyroscope.csv"
mag_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/orientation/physphox/Magnetometer.csv"

accel_df = pd.read_csv(acc_dat)
accel_df.columns = ['time','x','y','z']
accel_df['time'] = pd.to_datetime(accel_df['time'], unit='s')
accel_df = accel_df.set_index('time')

mag_df = pd.read_csv(mag_dat)
mag_df.columns = ['time','Bx','By','Bz']
mag_df['time'] = pd.to_datetime(mag_df['time'], unit='s')
mag_df = mag_df.set_index('time')

gyro_df = pd.read_csv(gyro_dat)
gyro_df.columns = ['time','wx','wy','wz']
gyro_df['time'] = pd.to_datetime(gyro_df['time'], unit='s')
gyro_df = gyro_df.set_index('time')

# Interpolate the magnetometer and gyroscope data to match the time stamps of the accelerometer data
accel_df_interp = accel_df.interpolate(method="time")
mag_df_interp = mag_df.interpolate(method="time")
gyro_df_interp = gyro_df.interpolate(method="time")

# Merge the three dataframes based on their timestamps using merge_asof()
# combined_df = pd.merge_asof(accel_df, mag_df_interp, on="time")
# combined_df = pd.merge_asof(combined_df, gyro_df_interp, on="time")
combined_df = pd.merge_asof(mag_df,accel_df_interp, on="time")
combined_df = pd.merge_asof(combined_df, gyro_df_interp, on="time")
combined_df = combined_df.dropna()



# Filter and plot headings from Magnetometer, gryoscope and Kalman Filtered fusion

df = combined_df

# Calculate the magnetic heading for each row
df['heading'] = np.degrees(np.arctan2(df['By'], df['Bx']))
df['heading'] += np.where(df['heading'] < 0, 360, 0)

# Compute the step direction using the accelerometer data
step_direction = np.degrees(np.arctan2(df['y'], df['x']))

df = df.set_index('time')
dt = df.index.to_series().diff().dt.total_seconds()[1:].mean() # get mean of the sample time from the collected data



# "Integrate" Gyroscope angular velocity (wz) data to get the change in heading over time
# dtheta = gyro_df.wz*dt
dtheta = df.wz*dt

# Compute the heading angle at each moment
df['g_heading'] = -(np.cumsum(dtheta)*180/3.1415925) + df['heading'][0] # heading as per the gyroscope (starts at 0) confert rads to degrees add magnetic field heading as starting point so they align

# Define the Kalman filter
kf = KalmanFilter(dim_x=3, dim_z=2)

# Define the state transition matrix
kf.F = np.array([[1, dt, 0], [0, 1, 0], [0, 0, 1]])

# Define the measurement function
kf.H = np.array([[1, 0, 0], [1, 0, 0]])

# Define the measurement noise covariance matrix
kf.R = np.diag([0.05, 0.05])

# Define the process noise covariance matrix
kf.Q = np.array([[0.1*dt**3/3, 0.1*dt**2/2, 0], [0.1*dt**2/2, 0.1*dt, 0], [0, 0, 0.01]])

# Set the initial state
kf.x = np.array([df.iloc[0]['heading'], 0, 0])

# Initialize the filtered data array
filtered_data = [kf.x]

# Loop over the data and apply the Kalman filter
for i in range(1, len(df)):
    # Compute the time difference since the last measurement
    dt = (df.index[i] - df.index[i-1]).total_seconds()

    # Predict the next state
    kf.predict()

    # Get the measurements
    measurement = np.array([df.iloc[i]['heading'], df.iloc[i]['g_heading']])

    # Update the state estimate with the measurements
    kf.update(measurement)

    # Add the filtered state to the filtered data array
    filtered_data.append(kf.x)

filtered_headings = [a[0] for a in filtered_data]
df['filtered_heading'] = filtered_headings

# plot
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(df.index, df.heading, label='Bfield Heading')
ax.plot(df.index, df.g_heading, label='gyro Heading')
ax.plot(df.index, filtered_headings, label='Filtered Heading')
ax.set_xlabel('Time')
ax.set_ylabel('Heading (degrees)')
ax.legend()
plt.show()


# Calculate steps from accelerometer, use this in conjunction with heading data to estimate a path taken

threshold = .52
fs = 6   # Sampling frequency
# fs = 1/(df.time[1] - df.time[0]) # technically this is right from matlab's terrible UI
f_low = 1.5  # lower bound of walking frequency
f_high = 2.5 # upper bound of walking speed


steps_df = prep_data(acc_dat,'')  #local uploads only, use this

calculated_steps = calculate_steps(threshold, steps_df)

df['steps'] = np.linspace(0,calculated_steps,df.shape[0])
df['steps_diff'] = df.steps[1] - df.steps[0]


# Plot heading from magnetometer, gyroscope and Kalman filtered data using both to normalize the data
g_coords = [(0, 0)]
b_coords = [(0, 0)]
f_coords = [(0, 0)]
for i in range(1, len(df)):
    x1, y1 = g_coords[i-1]
    x2, y2 = b_coords[i-1]
    x3, y3 = f_coords[i-1]
    b_angle = df['heading'][i]
    g_angle = df['g_heading'][i]
    f_angle = df['filtered_heading'][i]
    distance = df['steps_diff'][i]
    dx1, dy1 = distance * np.cos(np.deg2rad(g_angle)), distance * np.sin(np.deg2rad(g_angle))
    dx2, dy2 = distance * np.cos(np.deg2rad(b_angle)), distance * np.sin(np.deg2rad(b_angle))
    dx3, dy3 = distance * np.cos(np.deg2rad(f_angle)), distance * np.sin(np.deg2rad(f_angle))
    g_coords.append((x1+dx1, y1+dy1))
    b_coords.append((x2+dx2, y2+dy2))
    f_coords.append((x3+dx3, y3+dy3))

# Plot the path
fig, ax = plt.subplots()
x1, y1 = zip(*g_coords)
x2, y2 = zip(*b_coords)
x3, y3 = zip(*f_coords)
ax.plot(y1,x1)
ax.plot(y2,x2)
ax.plot(y3,x3)
plt.title('Path according to Gyroscope and Magnetometer')
plt.xlabel('Horizontal Steps taken')
plt.ylabel('Vertical Steps taken')
plt.legend(['Gyroscope','Magnetometer','Filtered'])
ax.set_aspect('equal')
plt.show()
