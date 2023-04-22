raw_dat = '130_hand_mag_gps.mat';
dat = load(raw_dat);

% convert struct to table
dataTable = struct2table(dat, 'AsArray', true);

% extract accelerometer and magnetometer data
accelData = dataTable.Acceleration;
magData = dataTable.MagneticField;

% extract separate columns for each axis
accelArray = table2array(accelData(:, {'X', 'Y', 'Z'}));
magArray = [magData.X{1:end}, magData.Y{1:end}, magData.Z{1:end}];