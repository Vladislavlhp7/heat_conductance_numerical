from helper import *


def explicit_method(step_space=20, step_time=100, HEAT_MIN=0, HEAT_MAX=500, L=100, time_limit=600, k=0.835, initial_cond: bool = False):
    """
    Compute the solution to the heat equation using the explicit method.

    The heat equation describes the diffusion of heat over time in a medium with constant thermal conductivity, and is given
    by the partial differential equation:

        dT/dt = k d^2T/dx^2

    where T is the temperature, t is time, x is position, and k is the thermal conductivity.

    The explicit method is a finite difference method that approximates the solution to the heat equation by discretizing
    the time and space domains and using a forward difference scheme to calculate the time derivative and a central
    difference scheme to calculate the second spatial derivative.

    Args:
        step_space (int): The spacing between adjacent points in the spatial domain. Default is 20.
        step_time (int): The spacing between adjacent points in the time domain. Default is 100.
        HEAT_MIN (float): The initial temperature at the boundary points of the spatial domain. Default is 0.
        HEAT_MAX (float): The initial temperature at all other points in the spatial domain. Default is 500.
        L (int): The length of the spatial domain. Default is 100.
        time_limit (int): The duration of the simulation. Default is 600.
        k (float): The thermal conductivity of the medium. Default is 0.835.

    Returns:
        cm (ndarray): A 2D array containing the computed values of the heat equation at each point in the
        time-space domain.

    Raises:
        None.

    Examples:
        # Compute the solution to the heat equation using default parameters.
        explicit_method()

        # Compute the solution to the heat equation with custom parameters.
        explicit_method(step_space=10, step_time=50, HEAT_MAX=1000, L=200, time_limit=1200, k=1.0)
    """
    # Set up the time and space ranges.
    time_range = np.arange(0, time_limit + step_time, step_time)
    space_range = np.arange(0, L + step_space, step_space)

    # Initialize the computational matrix with NaN values.
    cm = np.array([[np.nan for _ in space_range] for _ in time_range])

    # Set the initial boundary conditions for the computational matrix.
    if initial_cond:
        cm[0, :] = initial_temperature_distribution(x=np.array(space_range))
    else:
        cm[0, :] = HEAT_MAX
    cm[:, 0] = HEAT_MIN
    cm[:, -1] = HEAT_MIN

    # Calculate the explicit values from the computational molecules.
    lambda_ = step_time * k / step_space**2
    print(f'Δt={step_time}, Δx={step_space}, λ = {lambda_}')
    if 0 < lambda_ <= 0.5:
        print('Explicit solution is STABLE')
    else:
        print('Explicit solution is UNSTABLE')
    for j, t in enumerate(time_range[:-1]):
        for i, x in enumerate(space_range[:-1]):
            if np.isnan(cm[j+1][i]):
                cm[j+1][i] = lambda_ * (cm[j][i-1] + cm[j][i] * (-2.0 + 1.0/lambda_) + cm[j][i+1])

    return cm
