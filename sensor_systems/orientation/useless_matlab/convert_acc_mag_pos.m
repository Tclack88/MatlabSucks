raw_dat = '130_hand_mag_gps.mat'
dat = load(raw_dat)
% convert struct to table
dataTable = struct2table(dat,AsArray=true)

% extract accelerometer and magnetometer data
accelData = dataTable.Acceleration;
magData = dataTable.MagneticField;
posData = dataTable.Position;

% extract data from Timetables and convert to arrays
accelTimetable = accelData{1};
accelArray = table2array(accelTimetable(:, {'X', 'Y', 'Z'}));
magTimetable = magData{1};
magArray = table2array(magTimetable(:, {'X', 'Y', 'Z'}));
posTimetable = posData{1};
posArray = table2array(posTimetable(:, {'X', 'Y', 'Z'}));

% create separate columns for each axis
dataTable = array2table(accelArray, 'VariableNames', {'accel_x', 'accel_y', 'accel_z'});
dataTable.mag_x = magArray(:, 1);
dataTable.mag_y = magArray(:, 2);
dataTable.mag_z = magArray(:, 3);

% remove original columns
dataTable = removevars(dataTable, {'accelData', 'magData'});

% create separate columns for each axis
%dataTable.accel_x = accelData{:,1};
%dataTable.accel_y = accelData{:,2};
%dataTable.accel_z = accelData{:,3};
%
%dataTable.mag_x = magData{:,1};
%dataTable.mag_y = magData{:,2};
%dataTable.mag_z = magData{:,3};
%
% remove original accelerometer and magnetometer columns
%dataTable = removevars(dataTable, {'Acceleration', 'MagneticField'});