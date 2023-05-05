4 datasets have been provided. Data was collected by physphox, a superior sensor API to matlab.

Non-standard libraries:

filterpy
plotly

techinically non-standard, but very common or included in Anaonda or other IDEs
pandas
numpy
matplotlib


To run on a specific dataset, simply uncomment the dataset of interest at the top of the `orientation.py` file (The number represents the number of steps taken, necessary for checkin pedometer accuracy)

#data = 45 # actual path is not correct. Gyroscope is pretty accurate
#data = 157 #ok
#data = 175 #bad
data = 327
#data = 550 #kinda bad


