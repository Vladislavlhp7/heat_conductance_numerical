import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from implicit import implicit_method

colours = ['blue', 'green', 'orange', 'red']
linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
markers = ['o', 'v', 's', 'p', 'P', '*', 'X', 'D', 'd', '1', '2', '3', '4', '8', 'h', 'H', '+', 'x', 'X', 'D', 'd', '|',
           '_']


def lambda_search(step_spaces, step_times, fix: str = 'x', L=100, k=0.835, period=60, initial_cond=True):
    """
    Searches for the optimal value of lambda (k * step_time / step_space^2)
    for the given step space and step time values by applying the implicit
    method for the heat equation.

    Args:
    step_spaces (numpy.ndarray): array of step space values
    step_times (numpy.ndarray): array of step time values
    L (int): length of the spatial domain
    k (float): thermal diffusivity constant

    Returns:
    pandas.DataFrame: dataframe containing the final temperature values for
    each combination of step space and step time values and their corresponding
    lambda values.
    """
    final_temps = {}
    # Iterate through the step space and step time values
    for step_space, step_time in zip(step_spaces, step_times):
        # Calculate lambda for the given step space and step time values
        lambda_ = k * step_time / step_space ** 2

        # Generate a run ID using the step space, step time, and lambda values
        run_id = f'Δt={step_time}, Δx={np.round(step_space, 2)}'#, λ={np.round(lambda_, 2)}'

        # Apply the implicit method to calculate the final temperature values
        comp_matrix = implicit_method(step_time=step_time, step_space=step_space, initial_cond=initial_cond,
                                      time_limit=period, L=L, k=k)
        ind_at20 = int(comp_matrix.shape[1] / (L / step_space)) + 1
        temp_dist = comp_matrix[:, ind_at20]
        final_temps[run_id] = temp_dist

    fig, ax = plt.subplots()
    plt.grid()
    plt.title(f'Temperature Oscillations at Δx=20 for t <= 60s')
    plt.ylabel('Temperature (Celsius)')
    plt.xlabel('Time (s)')
    for i, (k, v) in enumerate(final_temps.items()):
        xrange = np.arange(0, period + 1e-1, step_times[i])
        markevery = 10 if step_times[i] < 1 else 1
        sns.lineplot(y=v, x=xrange, label=k, ax=ax, linestyle='dashed', color=colours[i],
                     marker=markers[i], markersize=5, markevery=markevery)
    ax.legend()
    plt.savefig(f'data/temp_dist_for_task3_fix_{fix}_initial_cond_{initial_cond}.png')
    # return final_temp_df


# comp_matrix_implicit = implicit_method(initial_cond=True)
# print(comp_matrix_implicit.round())

# comp_matrix_implicit = implicit_method(initial_cond=True, step_space=10)
# print(comp_matrix_implicit.round())
#
# comp_matrix_implicit = implicit_method(initial_cond=True, step_space=1, step_time=1)
# print(comp_matrix_implicit.round())

# comp_matrix_implicit = implicit_method(initial_cond=True, step_space=0.1, step_time=100)
# print(comp_matrix_implicit.round())

step_times = np.array([3, 1, 1e-1], dtype=float)
step_spaces = np.full_like(step_times, 1e-1)
lambda_search(step_spaces=step_spaces, step_times=step_times, fix='x', L=100, k=0.835, period=60, initial_cond=False)

step_spaces = np.array([1, 1e-1, 1e-2], dtype=float)
step_times = np.full_like(step_spaces, 1)
lambda_search(step_spaces=step_spaces, step_times=step_times, fix='t', L=100, k=0.835, period=60, initial_cond=False)

# print(implicit_method(step_space=10, step_time=300, initial_cond=True).round(2))
#
# print(implicit_method(step_space=10, step_time=100, initial_cond=True).round(2))
