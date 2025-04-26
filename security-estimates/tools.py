import math
import numpy as np


def rot_vector(arr):
    if arr.size < 2:
        return arr

    last_element = arr[-1, 0]  # 获取最后一个元素
    arr = np.insert(arr, 0, -last_element, axis=0)  # 将最后一个元素的相反数插入到第一位
    arr = np.delete(arr, -1, axis=0)  # 删除最后一个元素
    return arr


def rotate_matrix(column_vector):
    n = column_vector.shape[0]
    matrix = np.zeros((n, n))
    for i in range(n):
        if i == 0:
            matrix[:, i] = column_vector.flatten()
        else:
            rotated_vector = rot_vector(column_vector)
            matrix[:, i] = rotated_vector.flatten()
            column_vector = rotated_vector
    return matrix


def rot(column_vectors):
    processed_vectors = []

    for column_vector in column_vectors:
        processed_vector = rotate_matrix(column_vector)
        processed_vectors.append(processed_vector)

    result = np.vstack(processed_vectors)

    return result


def create_vector(d, eta, num_vectors):
    s = np.full((num_vectors, d), eta)
    return s


def create_column_vector(d, eta, num_vectors):
    s = np.full((d, num_vectors), eta)
    return [np.array([s[:, i]]).T for i in range(num_vectors)]


def max_singular_value(matrix):
    _, singular_values, _ = np.linalg.svd(matrix)
    max_singular_value = np.max(singular_values)
    return max_singular_value


def create_unit_vector(d):
    unit_vector = np.zeros((d, 1))
    unit_vector[0] = 1
    return unit_vector


def create_unit_vector_group(d, m):
    unit_vectors = []
    for i in range(m):
        unit_vector = np.zeros((d, 1))
        unit_vector[0] = 1
        unit_vectors.append(unit_vector)
    return unit_vectors
