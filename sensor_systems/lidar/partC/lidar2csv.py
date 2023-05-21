""" A python script to move that pesky matlab table of tables into individual csv files. I couldn't be bothered to figure out how to do this in matlab (unlike other simple load and export, it was more involved) so this code is brought to you almost exclusively by chatGPT"""

import scipy.io as sio
import numpy as np
import pandas as pd

## Load lidarData.mat
#data = sio.loadmat('lidarData.mat')
#
## Iterate over each element and save as CSV
#for i, element in enumerate(data['lidarData'][0]):
#    # Convert element to DataFrame
#    df = pd.DataFrame(element, columns=['x', 'y', 'z', 'intensity'])
#
#    # Save as CSV
#    df.to_csv(f'lidarData_{i+1}.csv', index=False)




# Load lidarLabel.mat
labels = sio.loadmat('lidarLabel.mat')

# Extract class names as a list
class_names = [label[0] for label in labels['lidarLabel'][0]]

### IMPORTANT #### If this is run again, delete the old 'lidarLabels.txt', otherwise thie existing one will just be appended to and more than 534 lines will exist. It's commented out now to avoid that
# with open('lidarLabels.txt','w+') as f:
#    for c in class_names:
#        f.write(f'{c}\n')
    
