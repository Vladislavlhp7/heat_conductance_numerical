import numpy as np
import matplotlib.pyplot as plt


def relative_error(vec):
    """
    Computes the relative error between elements of a vector.

    Args:
    - vec: A list or numpy array containing the elements of the vector.

    Returns:
    - errors: A list of the relative errors between adjacent elements of the vector.
    """
    errors = []
    for i in range(1, len(vec)):
        try:
            error = abs((vec[i] - vec[i - 1]) / vec[i])
        except ZeroDivisionError:
            error = 0.0
        errors.append(error)
    return errors



def absolute_error(v):
    """
    Computes the absolute error between adjacent elements of a vector `v`.

    Parameters:
    v (numpy array): Vector of numbers

    Returns:
    numpy array: Vector of absolute errors between adjacent elements of `v`
    """
    abs_errors = np.zeros_like(v)
    for i in range(1, len(v)):
        abs_errors[i] = np.abs(v[i] - v[i-1])
    abs_errors[0] = np.nan
    return abs_errors


def first_zero_index(lst):
    """
        Find the index of the first zero element in a list.

        Parameters:
        lst (list or ndarray): The list to search for a zero element.

        Returns:
        index (int or float): The index of the first zero element, or NaN if no zero element is found.

        Example:
        >>> first_zero_index([1, 2, 3, 0, 4, 5])
        3
    """
    try:
        return lst.index(0.0)
    except ValueError:
        return np.nan


def error_in_threshold(lst, threshold: float = 0.0001):
    """
        Find the first element in a list that is below a specified threshold.

        Parameters:
        lst (list or ndarray): The list to search for an element below the threshold.
        threshold (float): The threshold value.

        Returns:
        index (int or float): The index of the first element below the threshold, or NaN if no such element is found.
        element (float): The first element below the threshold, or NaN if no such element is found.

        Example:
        >>> error_in_threshold([0.1, 0.05, 0.001, 0.0005], threshold=0.01)
        (2, 0.001)
        """
    index = np.where(np.array(lst) <= threshold)[0]
    if index.size == 0:
        return np.nan, np.nan
    else:
        return index[0], lst[index[0]]


def plot_line(y, x=None, title='Relative Errors', xlabel='Element', ylabel='Error', color='tab:blue', grid=True):
    if x is None:
        x = np.arange(len(y))
    fig, ax = plt.subplots()
    ax.plot(x, y, color=color, linewidth=2)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis='both', which='major', labelsize=12)
    if grid:
        ax.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.show()


def plot_histogram(data, num_bins=10, color='blue', edgecolor='black', alpha=1.0, title='Histogram', xlabel='Value', ylabel='Frequency', grid=True):
    fig, ax = plt.subplots()
    ax.hist(data, bins=num_bins, color=color, edgecolor=edgecolor, alpha=alpha)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis='both', which='major', labelsize=12)
    if grid:
        ax.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.show()


def initial_temperature_distribution(x, g: float = 0.1, c: float = 400.0, L: int = 100):
    """
    Calculates the initial temperature distribution of a rod using the given parameters.

    Args:
        x (float or numpy.ndarray): The position(s) along the rod to calculate the temperature distribution at.
        g (float): The value of the function g in the initial boundary conditions. Defaults to 0.1.
        c (float): The value of the constant c in the initial boundary conditions. Defaults to 400.0.
        L (int): The length of the rod. Defaults to 100.

    Returns:
        float or numpy.ndarray: The initial temperature distribution at the given position(s).
    """
    return -x * g * (x - L) + c
