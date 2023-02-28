filename = '200steps_sensorlog_20230228_113308.mat';
raw_dat = load(filename);
dat = struct2table(raw_dat)
T = splitvars(dat)
writetable(T ,'200steps.csv');
