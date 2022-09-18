vidObj = VideoReader("pendulum_vid3.mp4");
[t, red, blue, dt]=ColourTracking(vidObj);

plot(t,blue)