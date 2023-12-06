import matplotlib.pyplot as plt
from scipy.io import mmread

def visualize_matrix_market(file_path):
    matrix = mmread(file_path)

    non_zero_indices = matrix.nonzero()
    rows, cols = non_zero_indices

    plt.scatter(cols, rows, s=5, c='blue', marker='o')
    plt.title(f'Non-Zero Pattern of Matrix {file_path}')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.grid(True)
    plt.show()

def visualize_matrix_market_binary(file_path):
    matrix = mmread(file_path)

    binary_matrix = np.zeros_like(matrix.toarray())
    binary_matrix[matrix.nonzero()] = 1

    plt.imshow(binary_matrix, cmap='Blues', interpolation='none', aspect='auto')
    plt.title(f'Non-Zero Pattern of Matrix {file_path}')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.show()

matrix_market_file_path = 'your_matrix_market_file.mtx'
visualize_matrix_market(matrix_market_file_path)