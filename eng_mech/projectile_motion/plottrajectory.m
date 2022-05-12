calib_frame=input('Enter the frame number for calibration: ');
image(read(obj,calib_frame))
disp('Click on up-left and down-right corners on the calibration paper and press Enter.'); 
[x_pix,y_pix]=getpts();
phys_X=input('Enter the horizontal physical distance between two cornenrs in cm:  ');
phys_Y=input('Enter the verticalal physical distance between two cornenrs in cm:  ');

first_frame=input('Enter the number of first frame of motion: ');
last_frame=input('Enter the number of last frame of motion: ');

n=0;
for i=first_frame:last_frame
    n=n+1;
    image(read(obj,i))
    disp('Click on ball and press Enter.')
    [x(n),y(n)]=getpts();
end

%using 10cmX10cm calibration paper

x_p=x.*phys_X/(x_pix(2)-x_pix(1));
y_p=y.*phys_Y/(y_pix(2)-y_pix(1));

X_Ball=max(x_p)-x_p;
Y_Ball=max(y_p)-y_p;

disp('The X and Y positions of the ball are saved as X_Ball and Y_Ball.')

figure(1)
plot(X_Ball,Y_Ball,'o')
xlabel('x(cm)')
ylabel('y(cm)')
        
        
        
        
        
        
        
        