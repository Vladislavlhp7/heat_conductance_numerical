import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import seaborn as sns
from analytical import analytical_method
from explicit import explicit_method
from implicit import implicit_method


methods = [analytical_method, implicit_method, explicit_method]
methods_str = ['Analytical', 'Implicit', 'Explicit']
colours = ['blue', 'green', 'red']
linestyles = ['solid', 'dashed', 'dotted', 'dashdot']


def experiment1(title: str, step_times, step_spaces, method, round_mode: int = 0):
    """
    Runs the given method for different combinations of step times and step spaces,
    and prints the resulting matrix with the corresponding Δt and Δx values.

    Args:
        title (str): A string to label the type of method being used.
        step_times (list): A list of step time values to use.
        step_spaces (list): A list of step space values to use.
        method (function): The numerical method to use.
        round_mode (int, optional): The number of decimal places to round the matrix values to.
            Defaults to 0.
    """
    for step_time, step_space in zip(step_times, step_spaces):
        comp_matrix = method(step_time=step_time, step_space=step_space)
        print(f"{title} Matrix: Δt = {step_time}, Δx = {step_space}")
        print(comp_matrix.round(round_mode))
        print()


def temperature_distribution_and_evolution(step_time: int = 100, period: int = 600, round_v: int = 6, L: int = 100,
                                           step_space: int = 20, methods_str=None):
    """
        Calculate temperature distribution and its evolution along a rod using analytical, implicit, and explicit methods.
        Store the temperature distribution at intervals, final temperature, and the evolution of the temperature.
        Save the final temperature distribution and the evolution of the temperature into csv files.
        Plot the temperature distribution along the rod at 100-second intervals, the temperature distribution for all
        three methods, the evolution of temperature for all three methods, and the relative error of the temperature
        evolution between implicit and explicit methods.

        Args:
            step_time (int): Time step. Default is 100.
            period (int): Period of time. Default is 600.
            round_v (int): Number of decimal places to round values. Default is 6.
            L (int): Length of the rod. Default is 100.
            step_space (int): Space step. Default is 20.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: The final temperature distribution and the evolution of the temperature.
        """
    if methods_str is None:
        methods_str = methods_str
    temp_dist_at_intervals = {}
    final_temps = {}
    temp_evols = {}
    for method, method_str in zip(methods, methods_str):
        comp_matrix = method()
        temp_dist_at_intervals[method_str] = comp_matrix
        final_temps[method_str] = comp_matrix[-1, :]
        temp_evols[method_str] = comp_matrix[:, 1]

    final_temps_df = pd.DataFrame().from_dict(final_temps)
    final_temps_df['x'] = np.arange(0, L + 1, step_space)
    final_temps_df.round(round_v).to_csv('data/temp_dist_600t.csv')

    temp_evol_df = pd.DataFrame().from_dict(temp_evols)
    temp_evol_df['t'] = np.arange(0, period + 1, step_time)
    temp_evol_df.round(round_v).to_csv('data/temp_evolution_at_x20.csv')

    # Plot the temperature distribution along the bar at 100s intervals
    for i, temp_dist_at_interval in enumerate(temp_dist_at_intervals.values()):
        fig, ax = plt.subplots()#figsize=(8, 8))
        plt.grid()
        plt.title(f'Temperature distribution at 100s intervals for {methods_str[i]} method')
        plt.xlabel('Length (cm)')
        plt.ylabel('Temperature')
        for j, temp_dist in enumerate(temp_dist_at_interval):
            sns.lineplot(x=np.arange(0, L + 1, step_space), y=temp_dist, marker='v',
                         label=methods_str[i] + f' t={j * step_time}s', ax=ax)
        ax.legend()
        # plt.ylim(325, 510)
        plt.savefig(f'data/temp_dist_at_intervals_{methods_str[i]}.png')

    # Plot the distribution of temperature
    fig, ax = plt.subplots()
    plt.grid()
    plt.title(f'Distribution of temperature along the rod with Δt = {step_time}, Δx = {step_space}')
    plt.xlabel('Length (cm)')
    plt.ylabel('Temperature')
    for i, method_str in enumerate(methods_str):
        sns.lineplot(data=final_temps_df, y=method_str, x='x', label=method_str, ax=ax, marker="o",
                     linestyle=linestyles[i], color=colours[i])
    ax.legend()
    plt.savefig('data/temp_dist_600t.png')

    # Plot the evolution temperature
    fig, ax = plt.subplots()
    plt.grid()
    plt.title(f'Evolution of temperature along the rod with Δt = {step_time}, Δx = {step_space}')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature')
    for i, method_str in enumerate(methods_str):
        sns.lineplot(data=temp_evol_df, y=method_str, x='t', label=method_str, ax=ax, marker="o",
                     linestyle=linestyles[i], color=colours[i])
    ax.legend()
    plt.savefig('data/temp_evolution_at_x20.png')

    # Calculate error
    relative_errors = [np.abs(temp_evols['Analytical'] - temp_evols['Implicit']) / temp_evols['Analytical'],
                       np.abs(temp_evols['Analytical'] - temp_evols['Explicit']) / temp_evols['Analytical']]

    # Plot the relative error
    fig, ax = plt.subplots()
    plt.grid()
    plt.title(f'Relative error of temperature evolution with Δt = {step_time}, Δx = {step_space}')
    plt.xlabel('Time (s)')
    plt.ylabel('Relative Error')
    for i, relative_error in enumerate(relative_errors):
        sns.lineplot(data=temp_evol_df, y=relative_error, x='t', label=methods_str[i + 1], ax=ax, marker="o",
                     linestyle=linestyles[i], color=colours[i + 1])
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
    ax.legend()
    plt.savefig('data/temp_evolution_error_at_x20.png')

    return final_temps_df, temp_evol_df


def compare_methods_for_different_step_sizes(step_spaces, step_times, method_str, period: int = 600, L: int = 100,):
    run_ids = []
    final_temps = {}
    temp_evols = {}
    method = methods[methods_str.index(method_str)]
    # max_step_sizes = np.max([step_spaces, step_times], axis=1)
    for step_space, step_time in zip(step_spaces, step_times):
        run_id = f'{method_str},Δt={step_time},Δx={step_space}'
        comp_matrix = method(step_time=step_time, step_space=step_space)
        final_temps[run_id] = comp_matrix[-1, :]
        step_20 = int(20 / step_space)
        temp_evols[run_id] = comp_matrix[:, step_20]
        run_ids.append(run_id)

    # trick the dataframe to store the data with biggest step sizes
    # final_temps_df = pd.DataFrame().from_dict(final_temps)
    # temp_evol_df = pd.DataFrame().from_dict(temp_evols)

    # Plot the distribution of temperature
    fig, ax = plt.subplots()
    plt.grid()
    plt.title(f'Distribution of temperature along the rod with `{method_str}` method')
    plt.xlabel('Length (cm)')
    plt.ylabel('Temperature')
    for i, (step_space, step_time) in enumerate(zip(step_spaces, step_times)):
        id = f' Δt={step_time},Δx={step_space}'
        x = np.arange(0, L + 1, step_space)
        y = final_temps[run_ids[i]]
        sns.lineplot(y=y, x=x, label=method_str + id, ax=ax, marker="v", linestyle=linestyles[i])
    ax.legend()
    plt.savefig(f'data/temp_dist_for_{method_str}.png')

    # Plot the evolution temperature
    fig, ax = plt.subplots()
    plt.grid()
    plt.title(f'Evolution of temperature along the rod with `{method_str}` method')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature')
    y = analytical_method()[:, 1]
    t = np.arange(0, period + 1, step_times[0])
    sns.lineplot(y=y, x=t, label=methods_str[0], ax=ax, marker="o", linestyle=linestyles[0])
    for i, (step_space, step_time) in enumerate(zip(step_spaces, step_times)):
        id = f' Δt={step_time},Δx={step_space}'
        t = np.arange(0, period + 1, step_time)
        y = temp_evols[run_ids[i]]
        sns.lineplot(y=y, x=t, label=method_str + id, ax=ax, marker="v",
                     linestyle=linestyles[i])
    ax.legend()
    plt.savefig(f'data/temp_evolution_for_{method_str}.png')
    # print(final_temps)


def main():
    k = 0.835
    HEAT_MAX = 500
    HEAT_MIN = 0

    L = 100
    step_space = 20

    period = 600
    step_time = 100

    step_times = [100, 50, 100]
    step_spaces = [20, 20, 10]

    # Analytical solution
    N = 10_000_001
    # experiment1(title='Analytical', method=analytical_method, step_times=step_times, step_spaces=step_spaces)
    #
    # # Explicit experiments
    # experiment1(title='Explicit', method=explicit_method, step_times=step_times, step_spaces=step_spaces)
    #
    # # Implicit experiments
    # experiment1(title='Implicit', method=implicit_method, step_times=step_times, step_spaces=step_spaces)

    temperature_distribution_and_evolution()
    compare_methods_for_different_step_sizes(step_spaces=step_spaces, step_times=step_times, method_str='Explicit')
    compare_methods_for_different_step_sizes(step_spaces=step_spaces, step_times=step_times, method_str='Implicit')


if __name__ == '__main__':
    main()
