% filename = '100step_sensorlog_20230302_100113.mat';
% filename = '250_pocket_pause_sensorlog_20230316_134817.mat';
%filename = 'Background_sensorlog_20230307_110521.mat'

filename = '/home/tclack/Downloads/100_run_pause_sensorlog_20230321_112444.mat'
filename = '/home/tclack/Downloads/180_run_sensorlog_20230321_112548.mat'
filename = '/home/tclack/Downloads/50_hand_slow_sensorlog_20230321_113224.mat'
% filename = '/home/tclack/Downloads/50_mix_sensorlog_20230318_130007.mat'
filename = '/home/tclack/Downloads/60_pocket_slow_sensorlog_20230321_112311.mat'
ilename = '/home/tclack/Downloads/70_pause_demo_sensorlog_20230321_112801.mat'
filename = '/home/tclack/Downloads/80_run_hand_sensorlog_20230321_114218.mat'

raw_dat = load(filename);
dat = struct2table(raw_dat);
T = splitvars(dat);
writetable(T ,'80_run_hand.csv');
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
