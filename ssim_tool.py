from __future__ import print_function
import pandas as _pd
from lib import violinplots, distribution, heatmaps, cabplots
from glob import glob as _glob
import re as _re
from ipywidgets import fixed, widgets
from lib import distribution as _dis


class SSIMAnalyse:

    """Class used to extract data from Morphology Files and Data Outputed by Surface-Surface Interaction Tool

    Attributes
    ----------
    violin_data : obj:'dataframe'
        Dataframe holding all the normalised interactions to the area of the smallest surface and converted to mJ/m^2 from kcal/mol.

    weighted_data: obj:'dataframe'
        Dataframe holding all the normalised interactions to the area of the smallest surface and converted to mJ/m^2 from kcal/mol.
        This has been weighted to the probabily of the two surfaces colliding. The probability is calculated by multiplying the surface %
        of each given surface with the one interacting

    whole_list : obj : 'str'
        List containing all the facet combinations from the morphology table supplied.

    col_options : obj : 'str'
        List containing all energy type calculations.

    data_all : obj : 'dataframe'
        Dataframe holding all the normalised interactions to the area of the smallest surface and converted to mJ/m^2 from kcal/mol.
        As well as extra information.

    """

    def __init__(self, morph_path, data_path, save_path=None,printing=False):
        """
        Initiates the Class

        :param str morph_path: Path required for the location of the Morphology File (example 'C:\Data\**\*LGA_Morphology*.csv')
        :param str data_path: Path required for the location of all the Data Files  (example ' C:\Data\**\*LGA*)_*.csv')
        :param str save_path: Path required for saving imags   (example ' C:\Data\') UNDER DEVELOPMENT
        """
        """ 

        Parameters:
        -----------
        morph_path (str) : Path required for the location of the Morphology File (example 'C:\Data\**\*LGA_Morphology*.csv')
        data_path (str) : Path required for the location of all the Data Files  (example ' C:\Data\**\*LGA*)_*.csv')

        """
        self.morph_path = morph_path
        self.data_path = data_path
        self.save_path = save_path
        self.violin_data = None
        self.data_all = None
        self.heatmap_data = None
        self.whole_list = None
        self.weighted_data = None
        self.col_options = None

        self.analyse(printing=printing)


    def morphology_extraction(self):
        """Extra MorphologyData and calculates combined probability.

        Extracts the Morphology Data from the morphology file by splitting the hkl into just the numbers. Cross references between
        two morphology files and creates every combination possible. If only 1 morphology file is given it duplicates it.

        :return: obj : Dictonary of all facet combinations and their associated probabilities
        """

        morph_files = [x for x in
                       _glob(self.morph_path, recursive=True)]  # Search through all folders and find all file names
        if len(morph_files) == 1:
            morph_files.append(morph_files[0])
        morph_data = []
        for m in morph_files:  # Go through all files
            morph_df1 = _pd.read_csv(m, usecols=['% Total facet area', 'hkl'])  # Read in every files into DataFrame
            morph_df1 = morph_df1[_pd.notnull(morph_df1['% Total facet area'])]  # Drops all rows which contain NaN
            morph_df1['hkl'] = morph_df1['hkl'].str.replace('({)|(})|(\s)',
                                                            '')  # Gets rid of white space and currly brackets
            morph_df1['Total facet area'] = morph_df1['% Total facet area'] / 100  # Convert % to absolute
            morph_df1 = morph_df1.set_index(morph_df1['hkl'])
            morph_data.append(morph_df1)  # Create list of dataframes
        columns = morph_data[0]['hkl'].tolist()  # Create lists of hkl for column and index
        index = morph_data[1]['hkl'].tolist()  #
        morph_df = _pd.DataFrame(index=index, columns=columns)  # Creates an empty dataframe
        # Create multiply all the probabilities together
        for i in columns:
            for j in index:
                product = morph_data[0].loc[i]['Total facet area'] * morph_data[1].loc[j]['Total facet area']
                morph_df.loc[j][i] = product
        morph_dict = morph_df.to_dict()  # convert to dictonary
        return morph_dict

    def analyse(self, printing=False):
        """
        Analyses the data from the files passes to the objects
        :param bool printing:  True will print all facets calculated plus the probability from surface areas
        :return: obj containg the analysed interaction data
        """
        def find_lowest(area1, area2):  # Find the lowest value of area
            if area1 == area2:
                return area1
            if area1 < area2:
                return area1
            else:
                return area2

        if self.data_path is not None:
            names = [x for x in _glob(self.data_path, recursive=True)]  # Parse All the file names into the names var
        else:
            print("Please point to correct folder as no files were found that end in .csv")
            return

        data_all = []
        print("Facets Processed")
        for n in names:
            t = _re.search('\(([0-9\-]+)\)\S+\(([0-9\-]+)\)', n)
            y = _re.search('\_([0-9\-]+)\S+\_([0-9\-]+)', n)
            if t:
                var = t.group(1) + "/" + t.group(2)
                temp_df = _pd.read_csv(n)
                if printing: print("Facets Interacting:", var)
                temp_df['Facet'] = var
                area = find_lowest(int(y.group(1)), int(y.group(2)))  # Finds out Area size based on file name
                if printing: print("Area used : " + str(area))
                temp_df['Total Energy'] = temp_df['Interaction Energy'].apply(lambda row: (row / area) * 6.94769E2)
                temp_df['Electrostatic'] = temp_df['Total ES Energy'].apply(lambda row: (row / area) * 6.94769E2)
                temp_df['Van der Waals'] = temp_df['Total VDW Energy'].apply(lambda row: (row / area) * 6.94769E2)
                temp_df['H-Bond'] = temp_df['Total HB Energy'].apply(lambda row: (row / area) * 6.94769E2)
                morph_dict = self.morphology_extraction()
                try:
                    probability = morph_dict[t.group(1)][t.group(2)]
                    if printing: print("       Probability:", probability)
                    temp_df['Weighted Total Energy'] = temp_df['Total Energy'].apply(lambda row: row * probability)
                    temp_df['Weighted Electrostatic'] = temp_df['Electrostatic'].apply(lambda row: row * probability)
                    temp_df['Weighted Van der Waals'] = temp_df['Van der Waals'].apply(lambda row: row * probability)
                    temp_df['Weighted H-Bond'] = temp_df['H-Bond'].apply(lambda row: row * probability)
                except:
                    if printing: print(
                        "Probability does not exist in table. Please check the current facets have been calculated")
                    data_all.append(temp_df)
                    continue
                data_all.append(temp_df)
            else:
                print("Failed on {} ".format(n))
        violin_data = _pd.concat(data_all)
        self.data_all = violin_data
        weighted_data = violin_data[
            ['Facet', 'Weighted Total Energy', 'Weighted Electrostatic', 'Weighted Van der Waals', 'Weighted H-Bond']]
        violin_data2 = violin_data[['Facet', 'Total Energy', 'Electrostatic', 'Van der Waals', 'H-Bond']]
        self.col_options = violin_data[['Total Energy', 'Electrostatic', 'Van der Waals', 'H-Bond',
                                        'Weighted Total Energy', 'Weighted Electrostatic', 'Weighted Van der Waals',
                                        'Weighted H-Bond']].columns
        self.whole_list = violin_data['Facet'].unique()
        self.weighted_data = weighted_data
        self.violin_data = violin_data2
        self.heatmap_data = violin_data
        print("Number of facet combinations : " + str(len(names)))

    def get_col_options(self, weighted=False):
        """Gets all col options based on if the weighted data has been enabled"""
        if weighted is False:
            col_options = self.col_options.drop(['Weighted Total Energy',
                                                 'Weighted Electrostatic',
                                                 'Weighted Van der Waals',
                                                 'Weighted H-Bond'])
        else:
            col_options = self.col_options.drop(['Total Energy', 'Electrostatic', 'Van der Waals', 'H-Bond'])
        return col_options

    def get_facet_list(self, weighted=False):
        """Sorts the facet list based on mean of each facet-facet interaction"""
        facet_list = self.violin_data.groupby('Facet')['Total Energy'].mean().sort_values().index
        if weighted:
            facet_list = self.weighted_data.groupby('Facet')['Weighted Total Energy'].mean().sort_values().index
        return facet_list

    def violinplots(self, sort=True, weighted=False, ylimit=(-40, 0), title="",bw=0.2, inner=None, orient="v"):
        """Function that generates the widgets required for plotting Violin plots.
         This Function passes the information to the graph drawer which generates the image

        :param bool sort: Passes True/False if the Plots need to be ordered by mean.
        :param bool weighted: Passes the True/False if the data is weighted by the size of the surface area.
        :param floats ylimit: Passes a list for scaling the y axis .
        :param str title: Passes the Title of the graph above the chart.
        :return: Violin plots
        """
        if weighted:
            data = self.weighted_data
            print("Weighted data selected!")
        else:
            print("Normal Data")
            data = self.violin_data
        facet_list = self.whole_list

        col_options = self.get_col_options(weighted=weighted)
        #print("Col options are :" + str(col_options))
        energy_comp = widgets.RadioButtons(options=col_options, description='Energy Components- Left Side')
        energy_comp2 = widgets.RadioButtons(options=col_options, description='Energy Components- Right Side')
        #save_clicker = widgets.Button(description="Save Graph")
        if sort:
            facet_list = self.get_facet_list(weighted=weighted)

        facet_list = widgets.SelectMultiple(options=facet_list, description='Facets Available', disabled=False)
        widgets.interact(violinplots.draw_violin,
                         data=fixed(data),
                         energy_comp=energy_comp,
                         energy_comp2=energy_comp2,
                         facet_list=facet_list,
                         ylimit=fixed(ylimit), title=fixed(title),
                         bw=fixed(bw),inner=fixed(inner), orient=fixed(orient));

    def heatmaps(self):
        """Function that generates the widgets required for plotting HeatMaps.
         This Function passes the information to the graph drawer which generates the image"""
        animation = 1
        slider=None
        if animation < 1:
            rot_slide = widgets.IntSlider(value=0, min=self.heatmap_data['Rotation'].min(),
                                          max=self.heatmap_data['Rotation'].max(),
                                          step=5, description='Set_Rotation', disabled=False,
                                          continuous_update=False, orientation='horizontal', readout=True,
                                          readout_format='d')
        else:
            rot_slide = widgets.Play(
                interval=525,
                value=0,
                min=self.heatmap_data['Rotation'].min(),
                max=self.heatmap_data['Rotation'].max(),
                step=5,
                description="Press play",
                disabled=False
            )
            slider = widgets.IntSlider(value=0,
                                       min=self.heatmap_data['Rotation'].min(),
                                       max=self.heatmap_data['Rotation'].max(),
                                       step=5, )
            widgets.jslink((rot_slide, 'value'), (slider, 'value'))
            widgets.HBox([rot_slide, slider])

        col_options = self.get_col_options(weighted=False)
        energy_comp = widgets.RadioButtons(options=col_options, description='Energy Components- Left Side',
                                           value='Total Energy')
        #save_clicker = widgets.Button(description="Save Graph")
        #display(save_clicker)
        facet_list = widgets.Dropdown(options=self.whole_list, description='Facets Avail', disabled=False)
        figt = widgets.interact(heatmaps.drawheatmap,
                                data=fixed(self.heatmap_data), energy_comp=energy_comp,
                                rot_slide=rot_slide, slider=slider,
                                facet_list=self.whole_list)

    def cab_extraction(self, excipient, weighted=False, exp_probe=False, median=False):
        """
        Method takes self to extract the data from the excipient substrates.

        :param obj excipient: This is passing the oject that the probe would be adhering to (Typically as the Excipeint)
        :param bool weighted: True - will pass the energy as a function of the surface area per crystal.
        :param bool exp_probe:  True - Switch on if using Excipient as probe and the self. is an excipient object.
        :param bool median: True - Will pass back data as a median and not mean
        :return: df_data : Dataframe holding statistical data for the given probe and excipients
        """

        # Checks to see if the weighted attr has been applied to weight for the surface area of each facet as a factor
        # of contribution.
        if weighted:
            energy = 'Weighted Total Energy'
            tempA = excipient.weighted_data
            facet_list1 = excipient.whole_list
            tempB = self.weighted_data
            facet_list2 = self.whole_list
        else:
            energy = 'Total Energy'
            tempA = excipient.violin_data
            facet_list1 = excipient.whole_list
            tempB = self.violin_data
            facet_list2 = self.whole_list

        probes = []
        adh_energy = []
        adh_std_energy = []
        coh_energy = []
        coh_std_energy = []
        for i in self.whole_list:
            t = i.split('/')
            if t[0] not in probes:
                probes.append(t[0])
        # If the plots are set as excipient probe than the reverse of the probes are recorded
        for x, layer in enumerate(probes):
            if exp_probe:
                tempT = tempA[tempA['Facet'].str.contains("/" + layer)][['Facet', energy]]
            else:
                tempT = tempA[tempA['Facet'].str.contains(layer + "/")][['Facet', energy]]
            tempY = tempB[tempB['Facet'].str.contains(layer)][['Facet', energy]]

            # The mean is defaulted as the descriptor for the energies if median is passed it uses that.
            if median:
                adh_energy.append(tempT[energy].median())
                coh_energy.append(tempY[energy].median())
            else:
                adh_energy.append(tempT[energy].mean())
                coh_energy.append(tempY[energy].mean())

            adh_std_energy.append(tempT[energy].std())
            coh_std_energy.append(tempY[energy].std())

        df_hold = _pd.DataFrame({'Probes': probes,
                                'Adhesion': adh_energy,
                                'Adhesion STD': adh_std_energy,
                                'Cohesion': coh_energy,
                                'Cohesion STD': coh_std_energy})
        return df_hold

    def cabplots(self, excipient, weighted=False, exp_probe=False, median=False,
                 title=None, label='APIvExp', xlim=(0, -10), ylim=(0, -10)):
        """
        This method takes the self as the probe and it uses the excipient (obj) provided to run the cab_extraction
        function which extacts the information from the probe/excipients for the statistical descriptor
        selected (Mean/Media). A graph is produced.
        :param (obj)    excipient   : This is passing the oject that the probe would be adhering to (Typically as the Excipeint)
        :param (bool)   weighted    : True - will pass the energy as a function of the surface area per crystal.
        :param (bool)   exp_probe   : True - Switch on if using Excipient as probe and the self. is an excipient object.
        :param (bool)   median      : True - Will pass back data as a median and not mean
        :param (str)    title       : Passes the name of the title
        :param (str)    label       : Passes the labels for the excipients
        :param (floats) xlim        : Passes a list of two floats to set the x axis
        :param (floats) ylim        : Passes a list of two floats to set the y axis
        :return: Graph of CAB Plot
        """
        data = self.cab_extraction(excipient=excipient, weighted=weighted, exp_probe=exp_probe, median=median)
        cabplots.drawcabplots(df_data=data, title=title, label=label, xlim=xlim, ylim=ylim)


def distributions(sets, labels, weighted=False, xlimit=(-10,0), kbw=0.2, tables=False, energy="Total Energy"):
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
    _dis.plot(sets=sets, labels=labels,
                      weighted=weighted, xlimit=xlimit,
                      kbw=kbw, tables=tables, energy=energy)


if __name__ == "__main__":
    print("ssim_tool")
