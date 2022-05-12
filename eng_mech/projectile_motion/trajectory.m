file_path=input('Enter the file path: ');
obj = VideoReader(file_path);
obj 
disp('Use image(read(obj,frame number) to find a frame for calibration and to indicate the range of frames for trajectory.')