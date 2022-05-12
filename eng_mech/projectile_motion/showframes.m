first_frame=input('Enter the number of first frame of motion: ');
last_frame=input('Enter the number of last frame of motion: ');
n=0;
for i=first_frame:last_frame
    image(read(obj,i))
    pause(0.2)
    
end