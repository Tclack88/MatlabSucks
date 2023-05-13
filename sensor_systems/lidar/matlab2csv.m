filename = 'pointcloudA.mat';
raw_dat = load(filename);
matrix = raw_dat.pointcloudA;
writematrix(matrix, 'pointcloud.csv');

