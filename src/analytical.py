import os
from math import sin, pi, exp, pow
from tqdm import tqdm
from helper import *


def analytical_method_update_scheme(x, t, N=10_000_001, L=100.0, k=0.835):
    """
        Computes the analytical solution of a one-dimensional heat equation for a homogeneous and isotropic medium.

        Args:
            x (float): The position at which to evaluate the solution, in centimeters.
            t (float): The time at which to evaluate the solution, in seconds.
            N (int, optional): The number of terms in the summation. Defaults to 101.
            L (float, optional): The length of the medium, in meters. Defaults to 100.0.
            k (float, optional): The thermal diffusivity of the medium, in square meters per second. Defaults to 0.835.

        Returns:
            float: The temperature at position `x` and time `t`, in degrees Celsius.

        Raises:
            AssertionError: If `N` is even.
    """
    assert N % 2 == 1, '`N` must be odd in order to compute the analytical solution'
    summation = 0.0

    px_over_L = pi * x / L
    ktpi2_over_L2 = k * t * pow(pi / L, 2)

    for n in range(1, N + 1, 2):
        summation += (1.0 / n) * sin(n * px_over_L) * exp(- n ** 2 * ktpi2_over_L2)
    res = (2000 / pi) * summation
    return res


def tune_analytical_solution():
    # Generate random x and t lists to test the analytical solution
    x_list = np.random.randint(0, 100, size=100)
    t_list = np.random.randint(0, 100, size=100)
    # Generate a list of odd N values to use in the analytical solution
    N_list = np.arange(1, 102, 2)
    # Initialize a list to store the indices of the first zero errors
    ind_list_error_zero = []

    # Loop over each x and t value
    for x in tqdm(x_list):
        for t in t_list:
            # Evaluate the analytical solution for each N value
            res = []
            for N in N_list:
                res.append(analytical_method_update_scheme(x=x, t=t, N=N))
            # Calculate the relative errors for each N value
            errors = relative_error(res)
            # Find the index of the first zero error
            ind_list_error_zero.append(first_zero_index(errors))

    # Convert the list of indices to a NumPy array
    ind_list_error_zero = np.array(ind_list_error_zero)
    # Plot a histogram of the indices of the first zero errors
    plot_histogram(data=ind_list_error_zero, num_bins=20, title='Minimum `N` term for Analytical Convergence')
    # plot_line(y=errors, title='Relative Errors: Tuning the analytical solution of the Heat equation')


def tune_analytical_solution_input_fixed(x: float, t: float) -> None:
    """
    Calculate the temperature distribution in a heat equation at a specific position and time using the analytical method
    with an updated scheme. This function is designed to help tune the parameters of the analytical method for optimal
    accuracy and efficiency.

    Parameters:
    x (float): The position at which to calculate the temperature distribution
    t (float): The time at which to calculate the temperature distribution

    Returns:
    None

    Example:
    >>> tune_analytical_solution_input_fixed(50, 300)
    """
    # Generate a list of odd N values to use in the analytical solution
    N_list = np.array([101, 1001, 5007, 10001, 50007, 100001, 500_007, 1_000_001, 5_000_007, 10_000_001])

    # Evaluate the analytical solution for each N value
    res = []
    for N in tqdm(N_list):
        res.append(analytical_method_update_scheme(x=x, t=t, N=N))

    # Calculate the absolute error between each temperature distribution and the exact solution
    errors = absolute_error(res)

    # Print the resulting temperature distributions and errors
    print(f"Temperature distributions: {res}")
    print(f"Absolute errors: {errors}")

    # Find the first N value that produces an error below a specified threshold
    threshold = 0.0001
    min_error_index, min_error_N = error_in_threshold(errors, threshold=threshold)
    if np.isnan(min_error_index):
        print(f"No N values produced an error below the threshold of {threshold}")
    else:
        print(f"The first N value that produces an error below the threshold of {threshold} is {min_error_N}")


def analytical_method(step_space=20, step_time=100, L=100, time_limit=600, k=0.835, N=10_000_001):
    """
        Calculate the temperature distribution in a heat equation using the analytical method with a fixed scheme.

        Parameters:
        step_space (int): The spatial step size.
        step_time (int): The time step size.
        L (float): The length of the rod.
        time_limit (int): The final time.
        k (float): The thermal diffusivity coefficient.
        N (int): The number of terms to use in the series expansion.

        Returns:
        comp_matrix (ndarray): The temperature distribution matrix.

        Example:
        >>> analytical_method(step_space=10, step_time=50, L=200, time_limit=800, k=0.5, N=10_000_001)
    """
    # Define the data folder path
    data_folder = "data"
    # Create the data folder if it does not exist
    os.makedirs(data_folder, exist_ok=True)
    filename = f"results_{step_space}_{step_time}_{L}_{time_limit}_{k}_{N}.out"
    # Join the data folder and filename paths
    filepath = os.path.join(data_folder, filename)

    # Check if the file exists
    if os.path.exists(filepath):
        comp_matrix = np.loadtxt(filepath)
    else:
        # Set up the time and space ranges.
        time_range = range(0, time_limit + step_time, step_time)
        space_range = range(0, L + step_space, step_space)
        # Initialize the computational matrix with NaN values.
        comp_matrix = np.array([[np.nan for _ in space_range] for _ in time_range])
        for i, t in tqdm(enumerate(time_range), desc='Calculating the Analytical solution'):
            for j, x in enumerate(space_range):
                comp_matrix[i][j] = analytical_method_update_scheme(x, t, L=L, k=k, N=N)

        # Save the NumPy array to file
        np.savetxt(filepath, comp_matrix)

    return comp_matrix
