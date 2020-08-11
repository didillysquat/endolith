import os
import pandas as pd
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import sys
import numpy as np


class Plotting:
    def __init__(self):
        self.seq_data_dir = sys.argv[1]
        self.all_sample_names = sorted(self._get_all_sample_names())
        self.abundance_df = self._make_abundance_df()
        self.fig = plt.figure(figsize=(10,10))
        gs = self.fig.add_gridspec(2, 2)
        self.ax_array_all = [self.fig.add_subplot(gs[0,i]) for i in range(2)]
        self.ax_array_minor = [self.fig.add_subplot(gs[1,i]) for i in range(2)]
        self.seq_c_dict = self._make_seq_c_dict()
        
    
    def _make_seq_c_dict(self):
        col_list, grey_list = self._get_colour_lists()
        return {seq: col_list[i] if i < len(col_list) else grey_list[i % len(grey_list)] for i, seq in enumerate(list(self.abundance_df))}

    def _plot_one_ax(self, ax, spec_char, minor_only=False):
        print(f'Processing "{spec_char}" samples')
        rect_patch_list = []
        rect_c_list = []
        x_coord = 0
        labs = []
        for sample in [_ for _ in self.all_sample_names if spec_char in _]: # Species specific sample list
            if sample in self.abundance_df.index:
                print(f'making Rectangle objects for {sample}')
                bottom = 0
                ser = self.abundance_df.loc[sample]
                seq_abunds = ser[ser > 0]
                if minor_only:
                    # We want to remove the most abundant sequence
                    seq_abunds = seq_abunds.reindex(seq_abunds.sort_values(ascending=False).index)[1:]
                    seq_abunds = seq_abunds.div(sum(seq_abunds))
                seq_abund_tup_list = [(seq, seq_abunds[seq]) for seq in [_ for _ in list(self.abundance_df) if _ in seq_abunds]]
                for seq, abund in seq_abund_tup_list:
                    rect_patch_list.append(Rectangle((x_coord, bottom), 10, abund, color=self.seq_c_dict[seq]))
                    rect_c_list.append(self.seq_c_dict[seq])
                    bottom += abund
                x_coord += 10
            else:
                # blank column for failed sample
                x_coord += 10
            labs.append(sample)
        # Here we have the rectangle patches done
        this_cmap = ListedColormap(rect_c_list)
        
        # Here we have a list of Rectangle patches
        # Create the PatchCollection object from the patches_list
        print('creating patch collection')
        patches_collection = PatchCollection(rect_patch_list, cmap=this_cmap)
        patches_collection.set_array(np.arange(len(rect_patch_list)))
        ax.set_ylim(0,1)
        ax.set_xlim(0, x_coord + 10)
        ax.set_xticks(range(5,x_coord + 5, 10))
        ax.set_xticklabels(labs, rotation='vertical')
        print('adding patch collection to axis')
        ax.add_collection(patches_collection)
        print('autoscaling')
        ax.autoscale_view()
    
    def plot(self):
        """
        Generate two plots, one for each of the samples. We will organise the samples
        according to their first character i.e. taxa they are supposed to be G or P
        Later we can conduct a dissimilarity for the samples, and then we can plot them according 
        to hierarchical clustering. For the time being though, let's simply plot them in an order
        sorted by name.
        Before doing the individual plotting create an abundance df where the sequences are the columns
        and the samples are the rows.
        Get the order of the sequences.
        Then on a persample basis per plot, make rectangles
        Then plot the rectangles
        """

        for ax, spec_char in zip(self.ax_array_all, ['G', 'P']):
            self._plot_one_ax(ax=ax, spec_char=spec_char)
        
        for ax, spec_char in zip(self.ax_array_minor, ['G', 'P']):
            self._plot_one_ax(ax=ax, spec_char=spec_char, minor_only=True)
        # Here we have both of the ax objects plotted up
        # now write out.
        plt.tight_layout()
        plt.savefig('stacked_all_seq.png', dpi=1200)
                
    def _make_abundance_df(self):
        """ Create an abundance dictionary that has sample as key and a dict as value
        The dict will have the nuc sequence as key and the absolute abundance as value
        This can then be converted into a df that we can use for the plotting
        """
        abund_df = {}
        for fasta_file in [_ for _ in os.listdir('.') if _.endswith('c.fasta')]:
            # The dictionary that will be the value for the given sample key
            sample_abund_dict = {}
            sample_name = fasta_file.split('.')[0].replace('_R1', '')
            names_file = fasta_file.replace('.unique.filtered.c.fasta', '.filtered.c.names')
            with open(names_file, 'r') as f:
                abs_abund_dict = {line.rstrip().split('\t')[0]: len(line.rstrip().split('\t')[1].split(',')) for line in f}
            with open(fasta_file, 'r') as f:
                fasta_as_list = [_.rstrip() for _ in f]

            for i in range(0, len(fasta_as_list), 2):
                sample_abund_dict[fasta_as_list[i+1]] = abs_abund_dict[fasta_as_list[i][1:]]
            
            # Here we have the sample_abund dict populated with absolute abundances
            abund_df[sample_name] = sample_abund_dict
        
        # Here we have the abundance dictionary populated.
        # Now generate a df from it
        df = pd.DataFrame.from_dict(abund_df, orient="index")
        df[pd.isna(df)] = 0

        # Now normalise the dict
        df = df.divide(df.sum(axis='columns'), axis='rows')
        # sort cols
        df = df.reindex(df.sum(axis='rows').sort_values(ascending=False).index, axis='columns')
        # sort rows
        df = df.sort_index()
        return df

    def _get_all_sample_names(self):
        return [_.replace('_R1.fastq.gz', '') for _ in os.listdir(self.seq_data_dir) if (('R1' in _) and (_.startswith('G') or _.startswith('P')))]

    def _get_colour_lists(self):
        colour_palette = self._get_colour_list()
        grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
        return colour_palette, grey_palette

    def _get_colour_list(self):
        colour_list = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5",
                "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693",
                "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900",
                "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
                "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744",
                "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68",
                "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED",
                "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
                "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
                "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
                "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
                "#C8D0F6", "#A3A489", "#806C66", "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59",
                "#8ADBB4", "#1E0200", "#5B4E51", "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94",
                "#7ED379", "#012C58", "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393",
                "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5", "#E773CE",
                "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4", "#00005F", "#A97399",
                "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058", "#A45B02",
                "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
                "#F4D749", "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9",
                "#FFFFFE", "#C6DC99", "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527",
                "#8BB400", "#797868", "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C",
                "#B88183", "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F", "#003109",
                "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E", "#1A3A2A", "#494B5A",
                "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F", "#BDC9D2", "#9FA064", "#BE4700",
                "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71", "#868E7E", "#98D058",
                "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F",
                "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
        return colour_list

p = Plotting()
p.plot()
