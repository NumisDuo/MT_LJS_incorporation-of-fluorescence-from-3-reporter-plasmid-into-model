#!/usr/bin/env python

"""
This script searches for all RMSE values in the log file, sorts the n best in ascending order and then extracts
the line with the RMSE as well as the preceding 18 lines and stores then in the file best.txt for further usage.
Plots of histogram distribution of paramerter estimations and latin hypercube sampling distributions are created
"""

# String operations
import re
import glob

# Plotting
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

# Data handling
import numpy as np
import csv
import pandas as pd


def main(plot=False, LHS=False, fromcsv=False, xyzplot=False, n=10):
    """
    Main function for reading and extraction the n best RMSE values and corresponding parameter values.
    :param n number of best RMSE vals to be processed. int or 'all'
    """

    # Find log filename
    logfiles = glob.glob("*.log")
    logfile = [filename for filename in logfiles if "ignore" not in filename][0]  # ignore some logs

    if not fromcsv:
        # Read logfile
        with open(logfile, "r") as f:
            txt = f.read()

        # Find RMSE values, split the actual value, sort and store the best n results
        RMSE_vals = re.findall("RMSE: [0-9.]*", txt)
        if n == 'all':
            n = len(RMSE_vals)
        RMSE_vals = sorted([float(RMSE.split(" ")[-1]) for RMSE in RMSE_vals])[:n]

        # Create search pattern
        pattern = re.compile(r'Simulator - <<< Parameters: >>>.*?RMSE: [0-9.]*', flags=re.S)

        # Loop over every RMSE
        res_lst = [match for val in RMSE_vals for match in re.findall(pattern, txt) if str(val) in match]

        # Store data
        data_dict = {}
        variables_lst = ['RMSE', 'qS_max', 'K_S', 'Ko', 'qAp_max', 'K_qS', 'Y_ASof', 'qm', 'Y_XSem', 'qAc_max', 'K_A',
                         'Ki_AS', 'Y_XA', 'Y_OS', 'Y_OA', 'Y_PS', 'dSox_P', 'sip', 'kp1', 'kla', 'Y_EX', 'q_Ed',
                         'M_Rb_m', 'M_Rb_c', 'M_Rb_a', 'M_Rb_b', 'q_Rbdeg', 'q_Rbmax', 'K_Rbsyn',
                         'qHSLR_max', 'qHSLR_csyn', 'qHSLR_deg', 'K_HSLR', 'K_thresh']

        for i in range(len(res_lst)):
            split_lst = list(res_lst[i].split("\n"))
            for j in range(len(split_lst)):
                split_vals = re.split(':| ', split_lst[j])
                for item in enumerate(split_vals):
                    if item[1] in variables_lst:
                        if item[1] not in data_dict:
                            data_dict[item[1]] = []
                        data_dict[item[1]].append(float(split_vals[-1]))

        # Rearrange RMSE to first position
        last_key, last_value = data_dict.popitem()
        data_dict = {last_key: last_value, **data_dict}

        # Save data dictionary to csv for further analysis
        dict_to_csv(data_dict, logfile)

    else:
        csvfiles = glob.glob("*.csv")
        csvfile = [file for file in csvfiles if "ignore" not in file][0]  # ignore some logs and select latest
        data_dict = pd.read_csv(csvfile)

    # Plot data
    if plot:       
        # Create a subplot for each variable
        fig = make_subplots(rows=len(data_dict), cols=1, subplot_titles=list(data_dict.keys()))

        # Add a histogram for each variable with x bins
        for i, (variable, values) in enumerate(data_dict.items()):
            fig.add_trace(go.Histogram(x=values, nbinsx=50), row=i+1, col=1)
            fig.update_xaxes(title_text='Parameter value', row=i+1, col=1)
            fig.update_yaxes(title_text='Frequency', row=i+1, col=1)
            # Add best result
            fig.add_vline(x=data_dict[variable][0],
                          line=dict(color='black', width=1.5),
                          row=i+1,
                          col=1,
                          annotation_text=f'{data_dict[variable][0]:.4f}',
                          annotation_position="top right")

        # Update layout
        fig.update_layout(template="simple_white", height=300*len(data_dict), width=800, showlegend=False)

        # Save the plot as an HTML file
        fig.write_html(f'{logfile[:-4]}_Histogram_plots.html')

    # Plot data
    if LHS:

        # Selection for ribosomal content

        # selection = ['q_Rbmax', 'K_Rbsyn', 'q_Rbdeg']  # segmon model
        # var_plot = [fr"$q_{{Rb,max}}\,[h^{{-1}}]$",
        #             fr"$K_{{M_{{Rb}}}}\,[h^{{-1}}]$",
        #             fr"$q_{{Rb,deg}}\,[h^{{-1}}]$"]

        selection = ['M_Rb_a', 'M_Rb_b']  # tanh model
        var_plot = [fr"$q_{{Rb,max}}\,[h^{{-1}}]$",
                    fr"$K_{{M_{{Rb}}}}\,[h^{{-1}}]$"]

        # Selection for HSLR

        # selection = ['qHSLR_max', 'qHSLR_csyn', 'qHSLR_deg', 'K_HSLR', 'K_thresh']  # sig_opt model
        # var_plot = [fr"$q_{{HSLR,max}}\,[h^{{-1}}]$",
        #             fr"$q_{{HSLR,const}}\,[h^{{-1}}]$",
        #             fr"$q_{{HSLR,deg}}\,[h^{{-1}}]$",
        #             fr"$K_{{HSLR}}\,[h^{{-1}}]$",
        #             fr"$K_{{thresh}}\,[h^{{-1}}]$"]

        # Colormap
        color_range = np.linspace(0, 1, len(data_dict['RMSE']))

        # Add LHS distribution
        fig = make_subplots(rows=len(selection), cols=1)

        for i, variable in enumerate(selection):
            if variable in selection:
                fig.add_trace(go.Scatter(x=data_dict['RMSE'],
                                         y=data_dict[variable],
                                         mode='markers',
                                         marker=dict(size=3, color=color_range, colorscale='Viridis')),
                              row=i + 1, col=1)

                fig.update_xaxes(title_text='RMSE [-]', row=i + 1, col=1)
                fig.update_yaxes(title_text=var_plot[i], tickmode='array', row=i + 1, col=1)
                # Add best result
                #fig.add_trace(go.Scatter(x=[data_dict['RMSE'][0]], y=[data_dict[variable][0]], mode='markers',
                #                         marker=dict(color='red', size=4)), row=i + 1, col=1,)


        # Update layout
        fig.update_layout(template="simple_white",
                          height=250 * len(selection),
                          width=800,
                          showlegend=False,
                          font_family='Open-Sherif',
                          )

        # Save the plot as an HTML file
        #fig.write_html(f'{logfile[:-4]}_LHS_plots.html')
        fig.write_image("LHS_tanh_analysis_1000.pdf")

        if xyzplot:
            plot_3d_from_csv(data_dict, 'qHSLR_csyn', 'qHSLR_max', 'RMSE')



def dict_to_csv(data_dict, logfile):
    # Get keys
    keys = data_dict.keys()

    # Get the length of the lists for all keys
    list_length = len(data_dict['RMSE'])

    # Prepare data for writing
    rows = [keys]  # Header row
    for i in range(list_length):
        row = [data_dict[key][i] for key in keys]
        rows.append(row)

    # Write to CSV
    with open(f"{logfile[:-4]}_best.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(rows)


def plot_3d_from_csv(data, x_param, y_param, z_param):
    # Extract selected parameters
    x_data = data[x_param]
    y_data = data[y_param]
    z_data = data[z_param]

    # Create a meshgrid
    X, Y = np.meshgrid(np.unique(x_data), np.unique(y_data))

    # Interpolate z values for all points on the meshgrid
    points = np.column_stack((x_data, y_data))
    Z = griddata(points, z_data, (X, Y), method='linear')

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

    # Set labels and title
    ax.set_xlabel(x_param)
    ax.set_ylabel(y_param)
    ax.set_zlabel(z_param)
    ax.set_title('3D Surface Plot of {} vs {} vs {}'.format(x_param, y_param, z_param))

    # Add color bar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

    # Show plot
    plt.show()

if __name__ == "__main__":
    main(plot=False, LHS=True, fromcsv=False, xyzplot=False, n='all')  # 'all' or int
