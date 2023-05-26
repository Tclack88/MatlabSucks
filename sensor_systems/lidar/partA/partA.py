import cv2
import pandas as pd
import numpy as np
import random

img_src = 'file1.png'



def extract_and_display_sift_keypoints(img_src,sigma,nOctaveLayers=3):
  img = cv2.imread(img_src)
  # Extract SIFT feature keypoints
  sift = cv2.xfeatures2d.SIFT_create(nOctaveLayers=nOctaveLayers, sigma=sigma)
  kp, des = sift.detectAndCompute(img, None)

  # Sort keypoints by SIFT score
  kp_sorted = sorted(kp, key=lambda x: x.response, reverse=True)

  # Display the ten strongest keypoints on the image
  img_with_kp = np.zeros_like(img) # new image for the keypoints
  img_with_kp = cv2.drawKeypoints(img, kp_sorted[:10], img_with_kp, color=(255, 255, 0), flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
  print('\n\nshowing image, press <q> to close and continue\n')
  cv2.imshow('Keypoints', img_with_kp) # local
  #cv2_imshow(img_with_kp) # colab patch
  cv2.waitKey(0)
  cv2.destroyAllWindows()


extract_and_display_sift_keypoints(img_src, 1.6,1)

extract_and_display_sift_keypoints(img_src, 5,1)


########### FEATURE MATCHING #################


img1 = cv2.imread('file1.png')
img2 = cv2.imread('file2.png')

# Extract SIFT feature descriptors from both images
sift = cv2.xfeatures2d.SIFT_create(nOctaveLayers=3, sigma=1.6)
kp1, des1 = sift.detectAndCompute(img1,None) # 2nd argument passes in a mask, but we are interested in the whole image, so no mask layer
kp2, des2 = sift.detectAndCompute(img2,None)


### "By Hand" ####

def match_features(des1, des2, threshold=0.8):
    """
    My match features (not using openCV built-in)
    Match SIFT feature descriptors from two images using sum of squared differences as distance metric.
    Return a list of matched features (as DMatch objects for convenient display later.
    Summary: -for each des1, compare to to each des2 to find the best sum-squared matches
             -go through des2 again to find 2nd best match (2nd best will be larger of course)
             -if ratio of best/second_best below some threshold (default .8), don't count it as a good match
             - Match indices are filtered at the end to check for unique one-to-one matches
             - only return "good" matches
    This takes about 2 minutes to run because I can't find a smarter way to avoid nested for loops
    """
    print('Matching featues manually, this will take some minutes....')
    matches = []
    for i in range(len(des1)):
        best_match_idx = -1 # no initial best match
        best_match_dist = float('inf')
        second_best_match_idx = -1 # no initial best match
        second_best_match_dist = float('inf')

        for j in range(len(des2)):
            dist = np.sum((des1[i] - des2[j])**2)

            if dist < best_match_dist:
                second_best_match_idx = best_match_idx
                second_best_match_dist = best_match_dist
                best_match_idx = j
                best_match_dist = dist
            elif dist < second_best_match_dist:
                second_best_match_idx = j
                second_best_match_dist = dist

        if best_match_dist/second_best_match_dist < threshold:
            matches.append(cv2.DMatch(i, best_match_idx, best_match_dist))

    # Filter for one-to-one matching
    filtered_matches = []
    matched_indices = set()
    for match in matches:
        if match.queryIdx not in matched_indices and match.trainIdx not in matched_indices:
            filtered_matches.append(match)
            matched_indices.add(match.queryIdx)
            matched_indices.add(match.trainIdx)

    return filtered_matches


matches = match_features(des1, des2)


def show_matches(imgs, kps, matches, rand_pts=0, n=None):
  """
  TODO: change to nicer display (overlap)
  show matches between images using openCV's built in for display
    imgs = tuple (img1,img2)
    matches = list of matches DMatch objects
  options:
    rand_pts = 0,1  (flag for choosing random points)
    n = int         (determine how many matches to show)
  """
  img1 = imgs[0]
  img2 = imgs[1]
  kp1 = kps[0]
  kp2 = kps[1]
  matchesN = matches
  if rand_pts:
    if n is None:
      print('n is none')
      n=10 # set a default value here in case it's forgotten and rondom is chosen
    print(f'Showing {n} random points')
    matchesN = random.sample(matches,n)
  else:
    print(f'Showing first {n} ponts')
  if n:
    matchesN = matches[:n]

  match_img = cv2.drawMatches(img1, kp1, img2, kp2, matchesN, None)
  #cv2_imshow(match_img)
  print('\n\nShowing image, press <q> to close and continue\n')
  cv2.imshow('Matches', match_img)
  cv2.waitKey(0)
  cv2.destroyAllWindows()


show_matches((img1,img2), (kp1,kp2), matches,0,10)


### Using openCV built-in module ###


def match_features_built_in(des1, des2, threshold=0.8):
    """
    Match SIFT feature descriptors from two images using OpenCV's built-in matcher.
    Return a list of "good" matches based on the  ratio of best/second_best as
    """
    # Normalize the feature descriptors
    des1_norm = cv2.normalize(des1, None, norm_type=cv2.NORM_L2)
    des2_norm = cv2.normalize(des2, None, norm_type=cv2.NORM_L2)

    # Brute Force (BF) matcher object
    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)

    matches = bf.match(des1_norm, des2_norm)
    # matches = sorted(matches, key=lambda x: x.distance)  # Sort the matches by distance

    good_matches = []
    for i in range(len(matches)-1):
        best_match = matches[i]
        second_best_match = matches[i+1]

        # Apply the distance ratio test
        if best_match.distance / second_best_match.distance < threshold:
            good_matches.append(best_match)

    good_matches = sorted(good_matches, key=lambda x: x.distance)  # Sort the matches by distance
    return good_matches

built_in_matches = match_features_built_in(des1,des2)


show_matches((img1,img2), (kp1,kp2), built_in_matches,0,10)
