import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Patch
import scipy.io as sio
from scipy.sparse import tril, triu

# get file from command line
import sys
file = sys.argv[1]

# custom colormap
cmap = colors.ListedColormap(['black', 'white'])

# load the Matrix Market file
matrix = sio.mmread(file)

# only reflect if the original matrix is lower triangular (csfmatrix stores lower triangular for symmetric matrices)
if matrix.shape[0] == matrix.shape[1] and (triu(matrix, 1).nnz == 0):
    # reflect the matrix
    matrix = triu(matrix, 1) + tril(matrix, -1).T

# create a binary mask of the matrix
mask = matrix.astype(bool).toarray()

# display the mask using matplotlib
plt.imshow(mask, cmap=cmap, interpolation='nearest')
plt.title(f'Nonzero Pattern of {file}')
plt.xlabel('Column')
plt.ylabel('Row')
legend_elements = [Patch(facecolor='black', edgecolor='black', label='Zero'),
                     Patch(facecolor='white', edgecolor='white', label='Nonzero')]
plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.show()