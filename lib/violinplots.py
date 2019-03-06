import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def draw_violin(data, energy_comp, energy_comp2, facet_list, ylimit, title,bw,inner,orient):
    """
    Draws violin plots based on widget selections and data set from SSIM.
    If two different energy compmenets are selected the violin plots are split in two.

    :param obj data: Interaction data for multiple facets interacting
    :param str energy_comp: Uses Widget energy compnent number 1
    :param str energy_comp2: Uses Widget energy component number 2
    :param str facet_list: uses widget list of facet available
    :param floats ylimit: Takes list of two number for the y axis scale
    :param str title:  Takes the title of the graph.
    :return: Violin Graphs
    """

    sns.set(style="whitegrid", palette="pastel", color_codes=True, context='paper', font_scale=2.1)
    a4_dims = (11.7, 8.27)
    g = plt.subplots(figsize=a4_dims)

    if energy_comp == energy_comp2:
        split_val = False
    else:
        split_val = True
    temp_list = []
    col_names = [energy_comp, energy_comp2]
   # print("DATA FILE = "+ str(data))
    for col in col_names:
    #    print("Col Name = "+ str(col))
        df_new = data[[col, 'Facet']]
        df_new = df_new.rename(columns={col: 'Interaction Energy(mJ/m^2)'})
        df_new['Energy Contribution'] = col
        temp_list.append(df_new)
    df_data = pd.concat(temp_list)

    g = sns.violinplot(x='Facet', y='Interaction Energy(mJ/m^2)', data=df_data, hue="Energy Contribution",
                       split=split_val, order=facet_list,
                       inner=inner, bw=bw, orient=orient)

    g.set_xticklabels(g.get_xticklabels(), rotation=30)
    g.set(ylim=ylimit)
    g.set_title(title)
    g.set_xlabel('Facets', fontsize=20)
    g.set_ylabel('Interaction Energy(mJ/m^2)', fontsize=20)
    g.tick_params(labelsize=15)
    g.legend(loc='lower right')
