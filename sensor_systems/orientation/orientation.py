from filterpy.kalman import KalmanFilter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt
import plotly.graph_objs as go
from pathlib import Path


# Try one of the several collected runs below


#data = 45 # actual path is not correct. Gyroscope is pretty accurate
#data = 157 #ok
#data = 175 #bad
data = 327 
#data = 550 #kinda bad


##### Step counter functions brought in from pedometer ######

def prep_data(fn):
  """Assumes pathlib object. Accelerometer data must contain number of actual
  steps taken for accuract"""
  data = fn
  df = pd.read_csv(data)
  df.columns = ['time','x','y','z']
  # Get actual steps (from fn title, leading numbers. eg. 100_pocket.csv -> 100)
  actual_steps = getnum(fn.name)
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


def load_data(num):
    """loads data in folder (number of folder indicates actual steps taken)"""
    acc_dat = Path.cwd().joinpath('physphox',f'{num}',f'{num}_Accelerometer.csv')
    gyro_dat = Path.cwd().joinpath('physphox',f'{num}','Gyroscope.csv')
    mag_dat = Path.cwd().joinpath('physphox',f'{num}','Magnetometer.csv')
    pos_dat = Path.cwd().joinpath('physphox',f'{num}','Location.csv')

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

    pos_df = pd.read_csv(pos_dat)
    pos_df.columns = ['time','lat','long','height','vel','dir','hor_acc','vert_acc']
    pos_df['time'] = pd.to_datetime(pos_df['time'], unit='s')

    return accel_df, mag_df, gyro_df, pos_df, acc_dat


def ApplyKalmanFilter(df):
    # 3-parameter Kalman filter
    kf = KalmanFilter(dim_x=3, dim_z=3)
    kf.F = np.array([[1, dt, 0], [0, 1, 0], [0, 0, 1]]) # state transition matrix
    kf.H = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0]]) # measurement function)
    kf.R = np.diag([0.05, 0.05, 0.05]) # measurement noise covariance matrix
    kf.Q = np.array([[0.1*dt**3/3, 0.1*dt**2/2, 0], [0.1*dt**2/2, 0.1*dt, 0], [0, 0, 0.01]]) # process noise covariance matrix
    kf.x = np.array([df.iloc[0]['Bx'], df.iloc[0]['By'], df.iloc[0]['wz']]) # Set the initial state
    filtered_data = [kf.x] # Initialize the filtered data array
    
    # Loop over the data and apply the Kalman filter
    for i in range(1, len(df)):
        kf.predict()
        measurement = np.array([df.iloc[i]['Bx'], df.iloc[i]['By'], df.iloc[i]['wz']])
        kf.update(measurement)
        # Calculate the heading from the filtered Bx, By, and wz
        Bx_f, By_f, wz_f = kf.x
        heading_f = np.arctan2(By_f, Bx_f) * 180 / np.pi
        filtered_data.append(np.array([heading_f, wz_f]))

    return filtered_data

def calculate_lat(heading,steps,lat0):
  # subtract 90 because 0 is north on a compass but the positive x axis 
  # by trig relations
  # add 11.85 because in Melbourne the magnetic declination is +11.85
  lat = [lat0]
  for i, h in enumerate(heading-90 + declination):
    dlat = np.cos(pi/180*h)*steps[i]*pace*lat_convert
    lat.append(lat[-1] + dlat)
  return np.array(lat[1:])

def calculate_long(heading,steps,long0):
  # subtract 90 because 0 is north on a compass but the positive x axis 
  # by trig relations
  # add 11.85 because in Melbourne the magnetic declination is +11.85
  lon = [long0]
  for i,h in enumerate(heading-90+ declination):
    dlong = np.sin(pi/180*h)*steps[i]*pace*long_convert
    lon.append(lon[-1] + dlong)
  return np.array(lon[1:])

#####################################




accel_df, mag_df, gyro_df, pos_df, acc_dat = load_data(data)

# Interpolate the magnetometer and gyroscope data to match the time stamps of the accelerometer data
accel_df_interp = accel_df.interpolate(method="time")
mag_df_interp = mag_df.interpolate(method="time")
gyro_df_interp = gyro_df.interpolate(method="time")

# Merge the three dataframes based on their timestamps using merge_asof()
combined_df = pd.merge_asof(mag_df,accel_df_interp, on="time")
combined_df = pd.merge_asof(combined_df, gyro_df_interp, on="time")
combined_df = combined_df.dropna()


# Filter and plot headings from Magnetometer, gryoscope and Kalman Filtered fusion

df = combined_df


declination = 11.85 #11 degrees 51 minutes   https://www.magnetic-declination.com/

# Calculate the magnetic heading for each row
df['heading'] = np.degrees(np.arctan2(df['By'], df['Bx'])) #- declination # negate because CW turn causes B field to move CCW and viceversa
# df['heading'] = -np.degrees(np.arctan2(df['By'], df['Bx'])) + declination # negate because CW turn causes B field to move CCW and viceversa
df['heading'] += np.where(df['heading'] < 0, 360, 0)
# df['heading'] -= initial_heading# - df['heading'][0]#align
# df['heading'] += initial_heading# - df['heading'][0]#align
11
# Compute the step direction using the accelerometer data
step_direction = np.degrees(np.arctan2(df['y'], df['x']))

df = df.set_index('time')
dt = df.index.to_series().diff().dt.total_seconds()[1:].mean() # get mean of the sample time from the collected data

# "Integrate" Gyroscope angular velocity (wz) data to get the change in heading over time
dtheta = df.wz*dt

# Compute the gyroscope heading angle at each moment
df['g_heading'] = -(np.cumsum(dtheta)*180/3.1415925) + df['heading'][0] # heading as per the gyroscope (starts at 0) confert rads to degrees add magnetic field heading as starting point so they align





# Kalman Filter
filtered_data = ApplyKalmanFilter(df)


# Extract the filtered heading estimates from the Kalman filter output
filtered_headings = [a[0] for a in filtered_data]
df['filtered_heading'] = filtered_headings


# Calculate steps separately (as per pedometer) and add this data to a new column
threshold = .42
fs = 6   # Sampling frequency
# fs = 1/(df.time[1] - df.time[0]) # technically this is more correct
f_low = 1.5  # lower bound of walking frequency
f_high = 2.5 # upper bound of walking speed

actual_steps, steps_df = prep_data(acc_dat)  # local uploads only here

calculated_steps = calculate_steps(threshold, steps_df)

df['steps'] = np.linspace(0,calculated_steps,df.shape[0])
df['steps_diff'] = df.steps[1] - df.steps[0]



# https://www.usgs.gov/faqs/how-much-distance-does-a-degree-minute-and-second-cover-your-maps

lat_convert = 3.28084/364000  # degrees/meter (1 degree/364000 ft)(3.28084 ft/m)
long_convert = 3.28084/288200 # degrees/meter (1 degree/288200 ft)(3.28084 ft/m)
pace = 438.33/532 # meters/step (measured from previous experiment in pedometer part)
pi = 3.14159

declination = 11.85 #11 degrees 51 minutes   https://www.magnetic-declination.com/




df['calc_lat'] = calculate_lat(df.g_heading,df.steps_diff,pos_df.lat[5])
df['calc_long'] = calculate_long(df.g_heading,df.steps_diff,pos_df.long[5])


df['latitudes1'] = calculate_lat(df.heading,df.steps_diff,pos_df.lat[5])
df['longitudes1'] = calculate_long(df.heading,df.steps_diff,pos_df.long[5])

df['latitudes2'] = calculate_lat(df.g_heading,df.steps_diff,pos_df.lat[5])
df['longitudes2'] = calculate_long(df.g_heading,df.steps_diff,pos_df.long[5])

df['latitudes3'] = calculate_lat(df.filtered_heading,df.steps_diff,pos_df.lat[5])
df['longitudes3'] = calculate_long(df.filtered_heading,df.steps_diff,pos_df.long[5])



trace1 = go.Scattermapbox(
    lat=pos_df['lat'],
    lon=pos_df['long'],
    mode='markers',
    marker=dict(
        size=8,
        color='green'
    ),
    name='Actual Path'
)


trace2 = go.Scattermapbox(
    lat=df['latitudes1'],
    lon=df['longitudes1'],
    mode='markers',
    marker=dict(
        size=8,
        color='yellow'
    ),
    name='Magnetometer'
)

trace3 = go.Scattermapbox(
    lat=df['latitudes2'],
    lon=df['longitudes2'],
    mode='markers',
    marker=dict(
        size=8,
        color='blue'
    ),
    name='Gyroscope'
)

trace4 = go.Scattermapbox(
    lat=df['latitudes3'],
    lon=df['longitudes3'],
    mode='markers',
    marker=dict(
        size=8,
        color='red'
    ),
    name='Kalman Filter'
)

layout = go.Layout(
    mapbox=dict(
        style='open-street-map',
        center=dict(
            lat=pos_df['lat'].mean(),
            lon=pos_df['long'].mean()
        ),
        zoom=16
    )
)

fig = go.Figure(data=[trace1,trace2,trace3,trace4], layout=layout)
fig.show()
