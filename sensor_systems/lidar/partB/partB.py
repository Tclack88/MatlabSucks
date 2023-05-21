import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from scipy.linalg import svd


# Load the pointcloud data from CSV
dat_src = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/sensor_systems/lidar/pointcloud.csv"
# dat_src = "pointcloud.csv" #local
pointcloud = pd.read_csv(dat_src)

# Extract the x, y, z, and intensity values from the DataFrame
x = pointcloud['x']
y = pointcloud['y']
z = pointcloud['z']
intensity = pointcloud['intensity']


def plot_pointcloud(x,y,z,intensity):
  """ 3D scatter plot using intensity as color"""
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(x, y, z, c=intensity, cmap='gray')
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  plt.title(f'Pointcloud plot using {len(x)} points')
  plt.show()



def downsample_pointcloud(dat,voxel_size=.25):
  pointcloud_df = pd.read_csv(dat)
  pointcloud = df[['x', 'y', 'z']].values
  voxel_size = voxel_size

  # Create NearestNeighbors object (with n=1 to get unique matches)
  nn = NearestNeighbors(n_neighbors=1)
  nn.fit(pointcloud)

  # Compute voxel centers
  min_coords = np.min(pointcloud, axis=0)
  max_coords = np.max(pointcloud, axis=0)
  voxel_centers = np.mgrid[min_coords[0]:max_coords[0]:voxel_size,
                          min_coords[1]:max_coords[1]:voxel_size,
                          min_coords[2]:max_coords[2]:voxel_size].reshape(3, -1).T

  # Find the nearest point in the point cloud for each voxel center
  distances, indices = nn.kneighbors(voxel_centers)

  # Filter any voxel centers matched to multiple nearest neighbors
  unique_indices, unique_counts = np.unique(indices, return_counts=True)
  filtered_indices = unique_indices[unique_counts == 1]
  filtered_voxel_centers = voxel_centers[filtered_indices]

  # Extract the downsampled point cloud
  downsampled_cloud = pointcloud[indices] # not strictly downsampled. a data point from the original set can still be matched to multiple voxel centers
  filtered_downsampled_cloud = pointcloud[filtered_indices] # remove multiple points assigned to the same voxel center
  filtered_intensities = pointcloud_df.intensity[filtered_indices]  # assuming you have an array named 'intensity'

  return filtered_downsampled_cloud, filtered_intensities

################################################################

plot_pointcloud(x,y,z,intensity) #unfiltered (original)


filtered_downsampled_cloud, filtered_intensities = downsample_pointcloud(dat_src,voxel_size=.25)

x = filtered_downsampled_cloud[:, 0]
y = filtered_downsampled_cloud[:, 1]
z = filtered_downsampled_cloud[:, 2]
plot_pointcloud(x,y,z,filtered_intensities)

# Define pointcloud A and pointcloud B from rotation (used in next part to verify ICP)
# Convert to numpy array (for easy manipulation)
pointcloudA = pointcloud[['x', 'y', 'z']].values

# Rotate and translate points
rotation_angle = np.radians(30)
rotation_matrix = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle), 0],
                            [np.sin(rotation_angle), np.cos(rotation_angle), 0],
                            [0, 0, 1]])
translation = np.array([2, 4, 0])


pointcloudB = np.dot(rotation_matrix, pointcloudA.T).T + translation

xb = pointcloudB[:, 0]
yb = pointcloudB[:, 1]
zb = pointcloudB[:, 2]
intensity = pointcloud.intensity # likely not relevant anymore, my plotting function just expects an "intensity entry"

plot_pointcloud(xb,yb,zb,intensity)


###### ICP procedure (using nearest neighbors)  ########

# convergence criteria
max_iterations = 100
threshold = 1e-6



def icp(pointcloudA, pointcloudB, num_iterations=max_iterations, threshold=threshold):
    # Initialize guesses
    # R = np.array([[np.cos(.5), -np.sin(.5), 0],
    #               [np.sin(.5), np.cos(.5), 0],
    #               [0, 0, 1]])
    # t = np.array([1.8, 3.8, 0])

    R = np.eye(3)
    # R = np.array([[.8,-.4,0],[.4,.8,0],[0,0,.9]])
    t = np.array([1,1,1]).T
    # t = np.array([1.8,3.8,.2]).T


    for i in range(num_iterations):
        # apply current rotation and translation
        transformedA = R.dot(pointcloudA.T).T + t.T

        # Find nearest neighbors
        # distances = np.zeros((transformedA.shape[0], pointcloudB.shape[0]))
        # for i in range(transformedA.shape[0]):
        #   for j in range(pointcloudB.shape[0]):
        #       diff = transformedA[i] - pointcloudB[j]
        #       distance = np.linalg.norm(diff)
        #       distances[i, j] = distance
        # indices = np.argmin(distances, axis=1)
        indices = np.argmin(np.linalg.norm(transformedA[:, None] - pointcloudB[None, :], axis=2), axis=1) #Accomplishes the same thing as the nested for loop above just in a smarter way

        print(len(indices))
        matched_pointsA = transformedA
        matched_pointsB = pointcloudB[indices]

        # Compute centroids
        mu_p = np.mean(matched_pointsA, axis=0)
        mu_q = np.mean(matched_pointsB, axis=0)

        # Find optimal R and t
        X = (matched_pointsA - mu_p).T # gives 3x434
        Y = (matched_pointsB - mu_q).T # gives 3x434

        U, _, Vt = svd(np.dot(X,Y.T)) # Singular Value Decomposition
        R = np.dot(Vt.T,U
        t = mu_q.reshape(3, 1) - R.dot(mu_p.reshape(3, 1)) #was (3,), correct errors with reshape

        # Check convergence to cancel early
        delta = np.linalg.norm(t, ord=None)
        if delta < threshold:
            print(f'ended after {i} iterations')
            break

    return R, t

R, t = icp(pointcloudA, pointcloudB)

# Print the recovered transformations
print(f"Recovered Rotation:\n{R}")
print(f"Recovered Translation:\n{t}")



# compare originalA, originalB and A trying to be rotated to fit B:
# TODO better to make under one plot with several subplots

# original A
xa1 = pointcloudA[:, 0]
ya1 = pointcloudA[:, 1]
za1 = pointcloudA[:, 2]
plot_pointcloud(xa1,ya1,za1,intensity)

# original B
xb1 = pointcloudB[:, 0]
yb1 = pointcloudB[:, 1]
zb1 = pointcloudB[:, 2]
plot_pointcloud(xb1,yb1,zb1,intensity)

# rotatd B (using ICP procedure)
pointcloudB2 = (np.dot(R, pointcloudA.T) + t).T
xb2 = pointcloudB2[:, 0]
yb2 = pointcloudB2[:, 1]
zb2 = pointcloudB2[:, 2]
plot_pointcloud(xb2,yb2,zb2,intensity)
