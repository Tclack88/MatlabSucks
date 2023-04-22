import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# This collected data is 90 degrees CCW, 90 degrees CW. 90 degrees CW and another 90 CW (for a total of 180).
# I ONLY HAVE THE HEADINGS RIGHT NOW, I WAS UNABLE TO FIGURE OUT THE KALMAN FILTERING
# TO PUT IT ON A PHYSICAL UP DOWN LEFT RIGHT DIRECTION

# Load the accelerometer and magnetometer data into separate dataframes
accel_df = pd.read_csv("Accelerometer.csv")
accel_df.columns = ['time','x','y','z']
accel_df['time'] = pd.to_datetime(accel_df['time'], unit='s')
accel_df = accel_df.set_index('time')

mag_df = pd.read_csv("Magnetometer.csv")
mag_df.columns = ['time','Bx','By','Bz']
mag_df['time'] = pd.to_datetime(mag_df['time'], unit='s')
mag_df = mag_df.set_index('time')

# Interpolate the magnetometer data to match the time stamps of the accelerometer data
mag_df_interp = mag_df.interpolate(method="time")

# Combine the two data frames based on the time stamps
combined_df = pd.merge_asof(accel_df, mag_df_interp, on="time")

# Remove missing values (Nans) (usually the first few columns of either if the
# time stamp value isn't 
df = combined_df.dropna()


# Calculate the magnetic heading for each row
df['heading'] = np.degrees(np.arctan2(df['By'], df['Bx']))
df['heading'] += np.where(df['heading'] < 0, 360, 0)

# Compute the step direction using the accelerometer data
step_direction = np.degrees(np.arctan2(df['y'], df['x']))

# Combine the step direction and magnetic heading to get the overall orientation
df['orientation'] = df['heading'] - step_direction
df['orientation'] = np.where(df['orientation'] < 0, df['orientation'] + 360, df['orientation'])

# Plot
fig, ax = plt.subplots()
ax.plot(df['orientation'])
ax.set_xlabel('Time')
ax.set_ylabel('Orientation')
ax.set_title('Magnetic Heading and Step Direction')
plt.show()
