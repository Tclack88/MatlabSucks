filename = '100step_sensorlog_20230302_100113.mat';
filename = '100steps_pausesensorlog_20230316_120210.mat';
%filename = 'Background_sensorlog_20230307_110521.mat'
raw_dat = load(filename);
dat = struct2table(raw_dat);
T = splitvars(dat);
writetable(T ,'100steps_pause.csv');
%writetable(T ,'100steps.csv');
%writetable(T ,'background.csv');

%% Following is an attempt to do fft using disgusting matlab. Better reserved for python
% x = T.X;
% t = T.Acceleration_Timestamp;
% 
% Y = fft(x);
% L = seconds(t(end) - t(1));
% Fs = 4;
% T = 1/Fs;
% 
% %t = (0:L-1)*T
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = Fs*(0:(L/2))/L;
% plot(f,P1)
