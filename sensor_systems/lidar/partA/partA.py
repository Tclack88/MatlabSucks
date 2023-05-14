import cv2
import pandas as pd
import numpy as np
import random

img_src = 'file1.png'




def extract_and_display_sift_keypoints(img_src,dat_src,sigma,nOctaveLayers=3):
  img = cv2.imread('file1.png')
  pointcloud = pd.read_csv(dat_src)

  # Extract SIFT feature keypoints
  sift = cv2.xfeatures2d.SIFT_create(nOctaveLayers=nOctaveLayers, sigma=sigma)
  kp, des = sift.detectAndCompute(img, None)

  # Sort keypoints by SIFT score
  kp_sorted = sorted(kp, key=lambda x: x.response, reverse=True)

  # Display the ten strongest keypoints on the image
  # Create a new image for the keypoints
  img_with_kp = np.zeros_like(img)

  # Draw the ten strongest keypoints on the image
  img_with_kp = cv2.drawKeypoints(img, kp_sorted[:10], img_with_kp, color=(255, 255, 0), flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)



  # cv2.imshow('Keypoints', img_with_kp)
  # img = cv2.imread("img.png")
  cv2.imshow("Imgage with kp", img_with_kp)
  #cv2_imshow(img_with_kp)
  cv2.waitKey(0)
  cv2.destroyAllWindows()

extract_and_display_sift_keypoints(img_src, dat_src, 1.6,1)

extract_and_display_sift_keypoints(img_src, dat_src, 5,1)


########### FEATURE MATCHING #################


img1 = cv2.imread('file1.png')
img2 = cv2.imread('file2.png')

# Extract SIFT feature descriptors from both images
sift = cv2.xfeatures2d.SIFT_create(nOctaveLayers=3, sigma=1.6)
kp1, des1 = sift.detectAndCompute(img1,None) # 2nd argument passes in a mask, but we are interested in the whole image, so no mask layer
kp2, des2 = sift.detectAndCompute(img2,None)


### "By Hand" ####
def match_features(des1, des2):
    """
    Match SIFT feature descriptors from two images using sum of squared differences as distance metric.
    Return a list of matched features (as DMatch objects for convenient display later.
    This takes about 2 minutes to run
    """
    matches = []
    for i in range(len(des1)):
        best_match_idx = -1 # no initial best match
        best_match_dist = float('inf')
        for j in range(len(des2)):
            dist = np.sum((des1[i] - des2[j]) ** 2)
            if dist < best_match_dist:
                best_match_idx = j
                best_match_dist = dist
        second_best_match_dist = float('inf')
        for j in range(len(des2)):
            if j != best_match_idx:
                dist = np.sum((des1[i] - des2[j]) ** 2)
                if dist < second_best_match_dist:
                    second_best_match_dist = dist
        if best_match_dist / second_best_match_dist < 0.8: # distance ratio test
            matches.append(cv2.DMatch(i, best_match_idx, best_match_dist))
    return matches


matches = match_features(des1, des2)


def show_matches(imgs, kps, matches, rand_pts=0, n=None):
  img1 = imgs[0]
  img2 = imgs[1]
  kp1 = kps[0]
  kp2 = kps[1]
  matchesN = matches
  if rand_pts:
    print('random')
    if n is None:
      print('n is none')
      n=10 # set a default value here in case it's forgotten and rondom is chosen
    matchesN = random.sample(matches,n)
  else:
    print('not random, chooseing first n ponts')
  if n:
    print('n exists, no need to change')
    matchesN = matches[:n]

  match_img = cv2.drawMatches(img1, kp1, img2, kp2, matchesN, None)
  cv2.imshow("Matches", match_img)
  #cv2_imshow(match_img)
  cv2.waitKey(0)
  cv2.destroyAllWindows()

show_matches((img1,img2), (kp1,kp2), matches,0,10)


### Using openCV built-in module ###
def match_features_built_in(des1, des2):
  # Normalize the feature descriptors
  des1_norm = cv2.normalize(des1, None, norm_type=cv2.NORM_L2)
  des2_norm = cv2.normalize(des2, None, norm_type=cv2.NORM_L2)

  # Create a brute-force matcher object
  bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)

  # Match the feature descriptors
  matches = bf.match(des1_norm, des2_norm)

  # Sort the matches by distance
  matches = sorted(matches, key=lambda x: x.distance)
  return matches


built_in_matches = match_features_built_in(des1,des2)


show_matches((img1,img2), (kp1,kp2), built_in_matches,0,10)
