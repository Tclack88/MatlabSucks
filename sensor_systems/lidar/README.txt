
Tested in a terminal environment on Ubuntu. May behave differently on other 
machines. (Such as how long an image is displayed for). During testing, using
<q> key will close the image and continue. It also has been tested in a
notebook environment (such as google colab). showing CV2 images requires using
cv2_imshow (commented out in the code) instead of cv2.imshow

`from google.colab.patches import cv2_imshow` may be necessary

some nonstandard libraries needed:

pip install opencv-python
pip install opencv-contrib-python # used for SIFT feature extraction

