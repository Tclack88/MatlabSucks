%filename = '200steps_sensorlog_20230228_113308.mat';
filename = '100step_sensorlog_20230302_100113.mat';
filename = 'Background_sensorlog_20230303_180709.mat'
raw_dat = load(filename);
dat = struct2table(raw_dat)
T = splitvars(dat)
writetable(T ,'background.csv');
%writetable(T ,'100steps.csv');