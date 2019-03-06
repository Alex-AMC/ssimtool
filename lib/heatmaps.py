import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def drawheatmap(data, energy_comp, facet_list, rot_slide, slider):
    """
    Draws heatmaps of XYE datasets where is a chosen energy for a chosen set of facets. Also plots a distribution plot
    :param obj data : Takes the heatmap data from SSIM analysis module
    :param str energy_comp: Used to select the energy to plot the heat map
    :param str facet_list:  Used to select the facet combination to plot
    :param int rot_slide: Used to select the rotation of the facets
    :param int slider:  Used to select the roation of the facets
    :return: Returns HeatMap plus distribution plot
    """
    sns.set(style="whitegrid", palette="pastel", color_codes=True)
    sns.set(context='paper')
    temp_list = []
    col_names = [energy_comp]
    for col in col_names:
        #print("Col Value = " + str(col))
        df_new = data[[col, 'Facet', 'X axis Displacement', 'Y axis Displacement', 'Rotation']]
        df_new = df_new.rename(columns={col: 'Interaction Energy(mJ/m^2)'})
        df_new['Energy Contribution'] = col
        temp_list.append(df_new)
    df_data = pd.concat(temp_list)
    to_matrix = df_data[(df_data['Facet'] == facet_list) & (df_data['Rotation'] == slider)][['X axis Displacement',
                                                                                             'Y axis Displacement',
                                                                                             'Interaction Energy(mJ/m^2)',
                                                                                             'Energy Contribution']]
    to_hist = to_matrix['Interaction Energy(mJ/m^2)']
    to_plot = to_matrix.pivot(index='Y axis Displacement',
                              columns='X axis Displacement',
                              values='Interaction Energy(mJ/m^2)')

    a4_dims = (16, 8)
    g, axs = plt.subplots(figsize=a4_dims, ncols=2)
    sns.set(context='paper', font_scale=2.1)
    g = sns.heatmap(to_plot, cbar_kws={'label': energy_comp}, cmap="RdBu_r", ax=axs[0])
    g = sns.distplot(to_hist, color="r", ax=axs[1])
    axs[0].set_xlabel('X axis Displacement(A)', fontsize=20)
    axs[0].set_ylabel('Y axis Displacement(A)', fontsize=20)
    axs[0].tick_params(labelsize=15)
    axs[1].tick_params(labelsize=15)
    axs[1].set_xlabel(energy_comp, fontsize=20)
    axs[1].set_ylabel('Density', fontsize=20)
    # axs[1].set(xlim(to_hist['Interaction Energy(mJ/m^2)'.min()]))
    axs[0].invert_yaxis()
    return g
