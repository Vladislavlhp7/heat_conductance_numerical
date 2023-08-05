import os

from tqdm import tqdm

from helper import *


def implicit_method(step_space=20, step_time=100, HEAT_MIN=0, HEAT_MAX=500, L=100, time_limit=600, k=0.835, initial_cond: bool = False):
    """
    Simulates the heat equation using an implicit finite difference method.

    Args:
    step_space (int): Space step size for discretization. Default is 20.
    step_time (int): Time step size for discretization. Default is 100.
    HEAT_MIN (int): Minimum heat value for boundary conditions. Default is 0.
    HEAT_MAX (int): Maximum heat value for boundary conditions. Default is 500.
    L (int): Length of the 1D space domain. Default is 100.
    time_limit (int): Time limit for simulation. Default is 600.
    k (float): Thermal conductivity of the material. Default is 0.835.

    Returns:
    comp_matrix (np.array): 2D numpy array with computed heat values at each spatial and time step.
    """
    lambda_ = k * step_time / step_space ** 2
    print(f'Δt={step_time}, Δx={step_space}, λ = {lambda_}')

    # Define the data folder path
    data_folder = "data"
    # Create the data folder if it does not exist
    os.makedirs(data_folder, exist_ok=True)
    filename = f"results_implicit_{step_space}_{step_time}_{L}_{time_limit}_{k}_{initial_cond}.out"
    # Join the data folder and filename paths
    filepath = os.path.join(data_folder, filename)

    # Check if the file exists
    if os.path.exists(filepath):
        comp_matrix = np.loadtxt(filepath)
    else:
        # Set up the time and space ranges.
        time_range = np.arange(0, time_limit + step_time, step_time)
        space_range = np.arange(0, L + step_space, step_space)

        # Initialize the computational matrix with NaN values.
        comp_matrix = np.array([[np.nan for _ in space_range] for _ in time_range])

        # Set the initial boundary conditions for the computational matrix.
        if initial_cond:
            comp_matrix[0, :] = initial_temperature_distribution(x=np.array(space_range))
        else:
            comp_matrix[0, :] = HEAT_MAX
        comp_matrix[:, 0] = HEAT_MIN
        comp_matrix[:, -1] = HEAT_MIN

        for j, t in tqdm(enumerate(time_range[:-1]), desc='Calculating Crank-Nicholson Matrices'):  # No need to compute for the last time step
            A = np.zeros(shape=(len(space_range), len(space_range)))
            B = np.zeros(shape=(len(space_range), len(space_range)))

            # for i, x in tqdm(enumerate(space_range), desc=f'Calculating Crank-Nicholson Matrix for t={t}'):
            for i, x in enumerate(space_range):
                if i == 0 or i == len(space_range) - 1:
                    A[i, i] = 1
                    B[i, i] = 1
                else:
                    A[i, i - 1] = -lambda_
                    A[i, i] = 2 + 2 * lambda_
                    A[i, i + 1] = -lambda_
                    B[i, i - 1] = lambda_
                    B[i, i] = 2 - 2 * lambda_
                    B[i, i + 1] = lambda_

            # Solve the system of linear equations: A * U_{j+1} = B * U_j
            comp_matrix[j + 1] = np.linalg.solve(A, B @ comp_matrix[j])

        # Save the NumPy array to file
        np.savetxt(filepath, comp_matrix)

    return comp_matrix
