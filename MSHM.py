# Import all the handy libraries.
import pandas as pd
import os
import plotly.express as px
import pickle

#pd.options.mode.copy_on_write = True
class MSHM:
    def __init__(self, path):
        """
        Mass Spect Heat Map (MSHM) plotter

        Consist of the following function:
            alias
            set_alias
            working_area
            set_working_area
            prep_dfs_for_plotting
            save_to_file
            load_from_file
            pfunction
        """
        # Set path
        self.path = path

        # Creat self.__dfs
        self.__dfs = self.__prep_dfs()
        self.__dfs_temp = {}

    def __get_folders(self):
        path = self.path
        folders = []
        for item in os.listdir(path):
            if os.path.isdir(os.path.join(path, item)):
                folders.append(item)
        return folders
        
        
        
    # This is a function that generate the array of intensity across mass for all deconvoluted samples. 
    def __prep_array(self):
        # Obtain the folder names, aka "sample names", from the directory.
        folders = self.__get_folders()
        # Create a empty dataframe to store data
        array = pd.DataFrame()
        for n in range(len(folders)):
            # Create an index.
            index = [folders[n]]
            # Read the data from _mass.txt file.
            data = pd.read_table(self.path+folders[n]+'/'+folders[n][:-12]+'_mass.txt',delim_whitespace=True, header=None, index_col=0).transpose().set_axis(index)
            # Concat to the array.
            array =pd.concat([array,data])
            #### End of the loop.
        array = array.divide(array.max(axis=1), axis = 0)
        array = array.sort_index(axis=1)
        return array 
    
    # This is a function that generate the array of intensity and molecular weight of all peaks.
    def __prep_peaks(self):
        # Obtain the folder names, aka "sample names", from the directory.
        folders = self.__get_folders()
        # Create a empty dataframe to store data
        array = pd.DataFrame()
        # A loop to obtain the peak data from each folder.
        for n in range(len(folders)):
            # Create a hierarchical index.
            tuples=[(folders[n][:-12], "Mass"), (folders[n][:-12], "Intensity")]
            index = pd.MultiIndex.from_tuples(tuples, names=["Sample", "Value"])
            # Read the "*_peaks.dat" file, which containing the "Mass" and "Intensity" data.
            peak_info = pd.read_table(self.path+folders[n]+'/'+folders[n][:-12]+'_peaks.dat', header=None, delim_whitespace=True).transpose().set_axis(index)
            # Add the above dataframe to the array.
            array =pd.concat([array,peak_info])
            #### End of the loop.
        # return result as array
        return array
    
    # Assemble array and peaks into a dictionary.
    def __prep_dfs(self):
        array = self.__prep_array()
        array.columns = array.columns.astype(int)
        peaks = self.__prep_peaks()
        # An extra one-column dataframe with alias.
        alias = pd.DataFrame(array.index, index = array.index, columns= ['Alias'])
        # Make a dictionary with working area
        mwrange= [array.columns[0],array.columns[-1]]
        aindex = array.index.to_list()
        working_area = {"mw_limit":mwrange,
                        "mw_range":mwrange,
                        "sample_alias":aindex,
                        "selected_alias":aindex
                        }
        # Intergrate everything into a dictionary.
        self.__dfs = {"array":array,"peaks":peaks,"alias":alias, "working_area":working_area}
        return self.__dfs
    
    # aquire the alias and return as a dictionary.
    def alias(self):
        """
        Display sample alias
        """
        return self.__dfs['alias'].to_dict()['Alias']
    
    # Input new alias as a dictionary.
    def set_alias(self,d):
        """
        Rename sample alias
        """
        # Pull out the old alias as dictionary
        d_old = self.alias()
        # Loop through all the keys in the new dictionary
        for i in d.keys():
            # For each key, replace the old value with the new value
            d_old[i] = d[i]
        # Turn the dictionary back to a dataframe
        alias = pd.DataFrame.from_dict(d_old,orient="index", columns=["Alias"])
        # Update the alias dataframe
        self.__dfs['alias'] = alias
        # Update sample_alias
        self.__dfs["working_area"]["sample_alias"] = alias.iloc[:,0].to_list()
        # Acquire the old selected alias (a list)
        selected_alias = self.__dfs["working_area"]["selected_alias"]
        # Replace the old alias with new ones
        temp = [
            alias.iloc[alias.index.get_loc(x), 0]
            if x in alias.index
            else x
            for x in selected_alias
            ]
        self.__dfs["working_area"]["selected_alias"] = temp
        
    # Get current working area.    
    def working_area(self):
        """
        Display current working area (mass range & selected sample alias)
        """
        return self.__dfs["working_area"]
    
    # Set current working area.
    def set_working_area(self,mw_limit = None,sample_name = None):
        """
        Set working area (mass range & selected sample alias)
        """
        working_area = self.__dfs["working_area"]
        if isinstance(mw_limit, list):
            if len(mw_limit) == 2 and all(working_area["mw_limit"][0] <= n <= working_area["mw_limit"][1] for n in mw_limit):
                working_area['mw_range'] = mw_limit
            else:
                print('Input is not a valid list. Ignoring.')
                pass
        else:
            print("Input is not a list. Ignoring.")
            pass
        
        if isinstance(sample_name, list):
            alias = self.__dfs['alias'].iloc[:,0].to_list()
            common_elements = [x for x in sample_name if x in alias]
            working_area["selected_alias"]= common_elements
        else:
            print("Input is not a list. Ignoring.")
            pass       
        
        self.__dfs["working_area"] = working_area
    
    # Assemble a new dictionay. Has array (with alias as index) and peaks. Both are sliced with working area.
    def prep_dfs_for_plotting(self, normalisation = True, deduction = 0, rank = 5, threhold = 10):
        """
        Prepare dataframe used for plotting figures. Can be use to generate the result table.
        """
        sample_name = self.__dfs["working_area"]["selected_alias"]
        mw_range = self.__dfs["working_area"]["mw_range"]
        
        
        array = self.__dfs['array']
        peaks = self.__dfs['peaks']
        alias = self.__dfs['alias']
        
        array.rename(index = {o:n for o, n in zip(array.index.to_list(), alias.iloc[:,0].to_list())}, inplace = True)
        peaks.rename(index = {o:n for o, n in zip(peaks.index.get_level_values(0).unique(), alias.iloc[:,0].to_list())}, level=0, inplace = True)
        
        array= array.loc[sample_name,:]
        peaks= peaks.loc[sample_name,:]
        
        left = mw_range[0]
        right = mw_range[1]
        # Set range for 'array'
        array_columns=array.columns[(array.columns >= left) & (array.columns <= right)]
        array = array[array_columns]
        # Set range for 'peaks'
        peaks_temp = pd.DataFrame()
        peaks_temp_slice = pd.DataFrame()
        for i in range(0,len(peaks.index),2):
            peaks_temp_slice = peaks.iloc[i:i+2,:]
            mask = (peaks_temp_slice.iloc[0,:] <= right) & (peaks_temp_slice.iloc[0,:] >= left)
            peaks_temp_slice = peaks_temp_slice.loc[:,mask]
            # remove column names
            peaks_temp_slice.columns = range(peaks_temp_slice.shape[1])
            # assemble the df
            peaks_temp =pd.concat([peaks_temp,peaks_temp_slice])
        peaks = peaks_temp
        
        # Picking peaksdata = MSHM("./Trouble shooting/")
        peaks_temp = pd.DataFrame()
        if normalisation == True:
            # list all the sample names into a list.
            samples_alias = array.index.to_list()
            # Start a loop to pick peaks
            for i in samples_alias:
                peakslice = peaks.loc[[i],:]
                # Pick whats greater than threhold
                greater = peakslice.loc[:, peakslice.iloc[1,:]>=threhold]
                # Pick top x
                top = greater.iloc[1].nlargest(rank).index
                peakslice = peakslice[top]
                sliceindex = peakslice.iloc[0].sort_values().index
                peakslice = peakslice[sliceindex]
                # Remove column names
                peakslice.columns = range(peakslice.shape[1])
                # Assemble the Dataframe
                peaks_temp = pd.concat([peaks_temp,peakslice])
            # Normolising
            # Sum of select peaks
            peaks_sum = peaks_temp[1::2].sum(axis=1).to_list()
            array_temp = array.div(peaks_sum,axis=0)*100
            peaks_temp[1::2] = peaks_temp[1::2].div(peaks_sum,axis=0)*100
            # do deduction
            array_temp.columns = array_temp.columns - deduction
            peaks_temp[0::2] = peaks_temp[0::2] - deduction
            # asign to array
            peaks = peaks_temp
            array = array_temp
        else:
            pass
        
        self.__dfs_temp = {"array":array,"peaks":peaks}

    def save_to_file(self, file_path):
        """
        Input a file path/name
        Not quite finished but good enough for exporting the database and analyse it on another machine
        """
        with open(file_path, 'wb') as output_file:
            pickle.dump(self.__dfs, output_file)
    def load_from_file(self, file_path):
        """
        Input a file path/name
        Not quite finished but good enough to open the database export by this script from another machine
        """
        with open(file_path, 'rb') as input_file:
            self.__dfs = pickle.load(input_file)
        
    # This is a function to plot the figure in plotly. 
    def pfunction(self, title="", color=[(0, "white"),(1, "black")], annotation = True, deduction = 0, normalisation = True, rank = 5, threhold = 10):
        """
        The plotter.
        Has 4 arguments: title, color, annotation, deduction and normalisation.

        title = "the title you want"
        color is black and white by default. Please check plotly documentation in order to change it.
        annotation = True by default.
        deduction = 0 by default. You can input your the base line value and annotate the peak by mass adduct.
        normalisation = True by default. It recalculates the relative peak hight.
        """

        self.prep_dfs_for_plotting(normalisation = normalisation, deduction = deduction, rank = rank, threhold = threhold)
        
        if deduction == 0:
            xtitle = 'Molecular Weight (Da)'
        else:
            xtitle = 'Molecular Weight Additive (Da)'
        
        df = self.__dfs_temp['array']
        pk = self.__dfs_temp['peaks']
        fig = px.imshow(df, color_continuous_scale=color)
        
        fig.add_shape(
            type="rect",
            x0=df.columns[0]-0.5,  # Start x-coordinate (adjusted to align with heatmap)
            y0=0-0.5,  # Start y-coordinate
            x1=df.columns[-1]+0.5,  # End x-coordinate
            y1=df.shape[0]-0.5,  # End y-coordinate
            line=dict(
                color="black",
                width=2
            ),
            fillcolor="rgba(0,0,0,0)"  # Transparent fill
        )
        for i in range(len(df.index)-1):
            fig.add_shape(
                type="line",
                x0=df.columns[0]-0.5,  # Start x-coordinate (adjusted to align with heatmap)
                y0=i+0.5,  # Start y-coordinate
                x1=df.columns[-1]+0.5,  # End x-coordinate
                y1=i+0.5,  # End y-coordinate
                line=dict(
                    color="black",
                    width=1
                )
            )
        if annotation == True:
            for i in range(len(df.index)):
                sample_name = df.index[i]
                peak = pk.loc[sample_name,:].dropna(axis = 1)
                for j in range(len(peak.columns)):
                    mass = peak.iloc[0,j].astype(int)
                    intensity = str(peak.iloc[1,j].round(1))+"%"
                    fig.add_annotation(
                            x=mass,  # x-coordinate where the annotation should be placed
                            y=i-0.2,
                            text= f"<span style='letter-spacing: -1px;'><b>{mass}</b></span>",
                            showarrow=False,
                            textangle=25,
                            font=dict(
                            color="#4A90E2",
                            size=9)
                        )
                    fig.add_annotation(
                            x=mass,  # x-coordinate where the annotation should be placed
                            y=i+0.2,
                            text=f"<span style='letter-spacing: -1px;'><b>{intensity}</b></span>",
                            showarrow=False,
                            textangle=25,
                            font=dict(
                            color="#ff66cc",
                            size=9)
                        )    
        fig.update_layout(
            width= 1200,
           # width=3000,
            height= len(df.index)*50+150,
            margin=dict(
                l=200,
                r=100,
                b=75,
                t=75
            ),
            autosize=False,
            title={
                'text': f'<b>{title}</b>', # Bold Text
                'xanchor': 'center',  # Anchor the title at its center
                'yanchor': 'middle',  # Anchor the title at its middle
                'y':(len(df.index)*50+120)/(len(df.index)*50+150),  # Put the center of the title 50 pixel away from the figure top border
                'x':0.5,  # Center the title
                'font': {
                    'size':24,
                    'color':"black"
                }
            },

            xaxis=dict(
                showline=True,  # Show the x-axis line
                linewidth=3,
                linecolor='black',
                mirror=True,     # Mirror the line to the top
                title=dict(
                    text=xtitle,
                    font=dict(size=20)
                ),
                title_standoff=20
            ),
            yaxis=dict(
                showline=True,  # Show the y-axis line
                linewidth=3,
                linecolor='black',
                mirror=True,     # Mirror the line to the right
                title=dict(
                    text='Samples',
                    font=dict(size=20)
                )
            ),
            coloraxis=dict(
                cmin=0,
                cmax=1,
                colorbar=dict(
                    tickfont=dict(size=18),
                    tickvals=[0, 1],  # Define specific tick positions
                    ticktext=["0%", "100%"]  # Corresponding labels
                    )
                )
            )
        fig.update_layout(
            xaxis_tickfont=dict(size=15),  # Set the font size for x-axis ticks
            yaxis_tickfont=dict(size=15)   # Set the font size for y-axis ticks
        )
        return fig
