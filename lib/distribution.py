from ssim_tool import *

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot(sets, labels, weighted=False, xlimit=(-10,0), kbw=0.2, tables=False, energy="Total Energy"):
    """ Function calculates the KDE of a  set of datasets and then plots them with labels.

        :param obj(s) sets  : will take list of objects to pass to plot
        :param str(s) labels: will take list of strings to plot. Make sure to use the order as objects passes.
        :param bool weighted: True for having data weighted based on surface area present between the two surfaces
        :param float xlimit : Takes list of floats for the range of the x axis .
        :param float kbw    : Sets the smoothness factor for the KDE
        :param bool tables  : True prints the table holding stats for all datasets
        :param str energy   : Takes string of energy type options: " Total Energy" , "Electrostatic" , "Van der Waals", "H-Bond"
        :return: object     : Distribution Graph
        """
    # This code creates a KDE plot for all the energies
    if sets is None:
        print("Please input objects containing ssim_tool.")
    table = pd.DataFrame(columns=labels)
    if weighted:
        energy = "Weighted " + str(energy)
    else:
        energy = energy
    x = 0
    sns.set(style="whitegrid", palette="pastel", color_codes=True)
    a4_dims = (8, 8)
    plt.subplots(figsize=a4_dims)
    full = sns.set(context='paper', font_scale=2.1)
    for x, data in enumerate(sets):
        exp = labels[x]
        if weighted:
            full = sns.kdeplot(data.weighted_data[energy], shade=True, bw=kbw, label=exp).set(xlim=xlimit)
            ttable = data.weighted_data[energy].describe()
            table[exp] = ttable
        else:
            full = sns.kdeplot(data.violin_data[energy], shade=True, bw=kbw, label=exp).set(xlim=xlimit)
            ttable = data.violin_data[energy].describe()
            table[exp] = ttable
    plt.legend(prop={'size': 16}, title='Surfaces')
    plt.title('Density Plot with Multiple Surfaces (' + energy + ')')
    plt.ylabel('Density', fontsize=12)
    plt.xlabel(energy, fontsize=12)
    if tables:
        return table
