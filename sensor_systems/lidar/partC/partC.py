# part C
import scipy.io as sio
import pandas as pd
from pathlib import Path
import numpy as np
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import accuracy_score
from sklearn import svm

def extract_features(df):   
    # Extract features from intensity column
    intensity_values = df.iloc[:, -1]  # Access the last column (intensity)

    # Compute some features from intensity
    mean_intensity = intensity_values.mean()
    std_intensity = intensity_values.std()
    max_intensity = intensity_values.max()



    # PCA for shape ceature
    coordinates = df.iloc[:, :-1].values  # Extract the 3D coordinates (x, y, z)
    pca = PCA(n_components=3)
    pca.fit(coordinates)
    eigenvalues = pca.explained_variance_


    # Assign the features to the matrix
    features = np.zeros(6)  # Initialize
    features[0] = mean_intensity
    features[1] = std_intensity
    features[2] = max_intensity
    features[3:6] = eigenvalues

    return features



# Extract class names into a list
label_dat = sio.loadmat('lidarLabel.mat')
labels = [label[0] for label in label_dat['lidarLabel'][0]]


data_list = [] # list of dictionaries turn to pd dataframe


for i, label in enumerate(labels):
    filepath = Path.cwd().joinpath('lidarData',f'lidarData_{i+1}.csv')
    label_df = pd.read_csv(filepath)
    # Extract features for the current lidar scan
    features = extract_features(label_df)
    # Create a dictionary with the data
    data_dict = {
        'file_path': filepath,
        # 'features': features,
        'mean_intensity': features[0],
        'std_intensity' : features[1],
        'max_intensity' : features[2],
        'PCA1' : features[3],
        'PCA2' : features[4],
        'PCA3' : features[5],
        'label': label
    }
    data_list.append(data_dict)

df = pd.DataFrame(data_list)


#### Now use SVM to split ####

random_state=1298022 #student ID


features = df[['mean_intensity', 'std_intensity', 'max_intensity', 'PCA1', 'PCA2', 'PCA3']]
labels = df['label']
X_train, X_test, y_train, y_test = train_test_split(features,labels,random_state=random_state)

## An attempt to stratified the data to account for the different size of the data 
## (eg. over 100 pedestrians but only about 20 buildings in the dataset
## But there was no difference and in fact it was slightly lower. Uncomment to try

#sss = StratifiedShuffleSplit(n_splits=1, test_size=0.3, random_state=random_state)
#train_indices, test_indices = next(sss.split(features, labels))
#X_train, y_train = features.iloc[train_indices], labels.iloc[train_indices]
#X_test, y_test = features.iloc[test_indices], labels.iloc[test_indices]

model = svm.SVC()

model.fit(X_train,y_train)
predictions = model.predict(X_test)
print('Accuracy Score (using intinsity + PCA1 + PCA2 + PCA3)')
print(accuracy_score(predictions,y_test))



# Try again using only intensity data
features1 = df[['mean_intensity', 'std_intensity', 'max_intensity']]
X_train1, X_test1, y_train1, y_test1 = train_test_split(features1,labels,random_state=random_state) # tried using a different means of splitting, but it wasn't much more accurate,

model1 = svm.SVC()

model1.fit(X_train1,y_train1)
predictions1 = model1.predict(X_test1)
print('Accuracy score (using intensity data only)')
print(accuracy_score(predictions1,y_test1))




# Try again using intensity data + PCA1
features2 = df[['mean_intensity', 'std_intensity', 'max_intensity','PCA1']]
X_train2, X_test2, y_train2, y_test2 = train_test_split(features2,labels,random_state=random_state) # tried using a different means of splitting, but it wasn't much more accurate,

model2 = svm.SVC()

model2.fit(X_train2,y_train2)
predictions2 = model2.predict(X_test2)
print('Accuracy Score (using intinsity + PCA1)')
print(accuracy_score(predictions2,y_test2))




# Try again using intensity data + PCA1 + PCA2
features3 = df[['mean_intensity', 'std_intensity', 'max_intensity','PCA1','PCA2']]
X_train3, X_test3, y_train3, y_test3 = train_test_split(features3,labels,random_state=random_state) # tried using a different means of splitting, but it wasn't much more accurate,

model3 = svm.SVC()

model3.fit(X_train3,y_train3)
predictions3 = model3.predict(X_test3)
print('Accuracy Score (using intinsity + PCA1 + PCA2)')
print(accuracy_score(predictions3,y_test3))


# When using all PCAs and adding sucessive ones from the intensity data usually
# brought it down. With the choice made, it actually rises very slightly
# But intensity alone does the bulk of the work
