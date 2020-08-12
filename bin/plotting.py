#!/usr/bin/env python3
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
import scipy.spatial.distance
import scipy.cluster.hierarchy


class Plotting:
    def __init__(self):
        self.seq_data_dir = sys.argv[1]
        self.all_sample_names = sorted(self._get_all_sample_names())
        self.abundance_df = self._make_abundance_df()
        self.fig = plt.figure(figsize=(10,10))
        gs = self.fig.add_gridspec(4, 2)
        self.ax_array_all = [self.fig.add_subplot(gs[0,i]) for i in range(2)]
        self.ax_array_minor = [self.fig.add_subplot(gs[1,i]) for i in range(2)]
        self.dendro_ax = self.fig.add_subplot(gs[2,:])
        self.similarity_sorted_stack_ax = self.fig.add_subplot(gs[3,:])
        self.seq_c_dict = self._make_seq_c_dict()
        # Dynamics
        self.rect_patch_list = [] 
        self.rect_c_list = []
        self.x_coord = 0
        self.labs = []
        # The dendro
        self.dendro = None
    
    def _make_seq_c_dict(self):
        col_list, grey_list = self._get_colour_lists()
        return {seq: col_list[i] if i < len(col_list) else grey_list[i % len(grey_list)] for i, seq in enumerate(list(self.abundance_df))}

    def _make_rectangles(self, minor_only, sample, dendro=False):
        if sample in self.abundance_df.index:
            print(f'making Rectangle objects for {sample}')
            bottom = 0
            
            if minor_only:
                # We want to remove the most abundant sequence
                if dendro:
                    # If we are working with the dendro stacked bars, we remove the two most abundant
                    # seqs from the dataset as this is what the dendro was calculated on. This has
                    # to be done before sorting as the second most abundant sequence in a given sample
                    # is not necessarily one of the most abundant sequence in the dataset
                    ser = self.abundance_df.loc[sample][2:]
                    seq_abunds = ser[ser > 0]
                    seq_abunds = seq_abunds.reindex(seq_abunds.sort_values(ascending=False).index)
                    seq_abunds = seq_abunds.div(sum(seq_abunds))
                else:
                    ser = self.abundance_df.loc[sample]
                    seq_abunds = ser[ser > 0]
                    seq_abunds = seq_abunds.reindex(seq_abunds.sort_values(ascending=False).index)[1:]
                    seq_abunds = seq_abunds.div(sum(seq_abunds))
            else:
                ser = self.abundance_df.loc[sample]
                seq_abunds = ser[ser > 0]
            seq_abund_tup_list = [(seq, seq_abunds[seq]) for seq in [_ for _ in list(self.abundance_df) if _ in seq_abunds]]
            for seq, abund in seq_abund_tup_list:
                self.rect_patch_list.append(Rectangle((self.x_coord, bottom), 10, abund, color=self.seq_c_dict[seq]))
                self.rect_c_list.append(self.seq_c_dict[seq])
                bottom += abund
            self.x_coord += 10
        else:
            # blank column for failed sample
            self.x_coord += 10
        self.labs.append(sample)

    def _apply_rectangles(self, ax):
        # Here we have the rectangle patches done
        this_cmap = ListedColormap(self.rect_c_list)
        
        # Here we have a list of Rectangle patches
        # Create the PatchCollection object from the patches_list
        print('creating patch collection')
        self.patches_collection = PatchCollection(self.rect_patch_list, cmap=this_cmap)
        self.patches_collection.set_array(np.arange(len(self.rect_patch_list)))
        ax.set_ylim(0,1)
        
    def _format_ax_one_plot(self, ax):
        ax.set_xlim(0, self.x_coord)
        ax.set_xticks(range(5, self.x_coord + 5, 10))
        ax.set_xticklabels(self.labs, rotation='vertical')
        print('adding patch collection to axis')
        ax.add_collection(self.patches_collection)
        print('autoscaling')
        ax.autoscale_view()

    def _format_ax_all_samples(self, ax):
        ax.set_xlim(self.dendro_ax.get_xlim())
        ax.set_xticks([])
        # ax.set_xticks(range(5,int(self.dendro_ax.get_xlim()[1]),10))
        # ax.set_xticklabels(self.labs, rotation='vertical')
        print('adding patch collection to axis')
        ax.add_collection(self.patches_collection)
        print('autoscaling')
        ax.autoscale_view()

    def _apply_rectangles_all_samples(self, ax):
        # Here we have the rectangle patches done
        this_cmap = ListedColormap(self.rect_c_list)
        
        # Here we have a list of Rectangle patches
        # Create the PatchCollection object from the patches_list
        print('creating patch collection')
        patches_collection = PatchCollection(self.rect_patch_list, cmap=this_cmap)
        patches_collection.set_array(np.arange(len(self.rect_patch_list)))
        ax.set_ylim(0,1)
        ax.set_xlim(0, self.x_coord + 10)
        ax.set_xticks(range(5, self.x_coord + 5, 10))
        ax.set_xticklabels(self.labs, rotation='vertical')
        print('adding patch collection to axis')
        ax.add_collection(patches_collection)
        print('autoscaling')
        ax.autoscale_view()

    def _reset_dynamics(self):
        self.rect_patch_list = []
        self.rect_c_list = []
        self.x_coord = 0
        self.labs = []

    def _plot_one_ax(self, ax, spec_char, minor_only=False):
        print(f'Processing "{spec_char}" samples')
        self._reset_dynamics()
        for sample in [_ for _ in self.all_sample_names if spec_char in _]: # Species specific sample list
            self._make_rectangles(minor_only=minor_only, sample=sample)
        self._apply_rectangles(ax=ax)
        self._format_ax_one_plot(ax=ax)
    
    def _plot_one_ax_all_samples(self):
        print(f'Processing all samples by similarity')
        self._reset_dynamics()
        for sample in self.dendro['ivl']:
            self._make_rectangles(minor_only=True, sample=sample, dendro=True)
        self._apply_rectangles(ax=self.similarity_sorted_stack_ax)
        self._format_ax_all_samples(ax=self.similarity_sorted_stack_ax)

    def _dendro_plot(self):
        # Remove 2 most abund seqs
        dist_abund_df = self.abundance_df.iloc[:,2:]
        # remove samples that only contained the most abund seq
        dist_abund_df = dist_abund_df.loc[(dist_abund_df != 0).any(axis=1)]
        dist_condensed = scipy.spatial.distance.pdist(dist_abund_df, metric='braycurtis')
        dist_square = scipy.spatial.distance.squareform(dist_condensed)
        dist_df = pd.DataFrame(dist_square, index=dist_abund_df.index, columns=dist_abund_df.index)
        # dist_df.to_csv(self.dist_out_path, index=True, header=True, compression='infer')
        # condensed_dist = scipy.spatial.distance.squareform(dist_df)
        linkage = scipy.cluster.hierarchy.linkage(y=dist_condensed, optimal_ordering=True)
        self.dendro = scipy.cluster.hierarchy.dendrogram(linkage, ax=self.dendro_ax, labels=list(dist_df.index), link_color_func=lambda k: 'black')
        self.dendro_ax.collections[0]._linewidths=np.array([0.5])
        self.dendro_ax.spines['right'].set_visible(False)
        self.dendro_ax.spines['top'].set_visible(False)
        self.dendro_ax.spines['bottom'].set_visible(False)
        self.dendro_ax.spines['left'].set_visible(False)
        self.dendro_ax.set_yticks([])

    def plot(self):
        """
        Generate 6 plots.
        Row one, all seqs, left and right split by species
        Row two, minor seqs, left and right split by species
        Row three hierarchical clustering
        Row four plotted in the order of the hierarchical clustering all together
        """

        # Row one
        for ax, spec_char in zip(self.ax_array_all, ['G', 'P']):
            self._plot_one_ax(ax=ax, spec_char=spec_char)
            ax.set_title(f'All sequences "{spec_char}" samples', fontsize='small')
        
        # Row two
        for ax, spec_char in zip(self.ax_array_minor, ['G', 'P']):
            self._plot_one_ax(ax=ax, spec_char=spec_char, minor_only=True)
            ax.set_title(f'Maj seq removed "{spec_char}" samples', fontsize='small')
        
        # Row three
        # Do braycurtis distances
        self._dendro_plot()
        self.dendro_ax.set_title(f'Bray-Curtis dissimilarity (computed on all sequences less the 2 most abund of dataset) hierarchical clustering', fontsize='small')

        # Row four
        self._plot_one_ax_all_samples()
        self.similarity_sorted_stack_ax.set_xlabel('2 most abund data set seqs removed all samples same order as hierarchical clustring', fontsize='small')

        # Here we have both of the ax objects plotted up
        # now write out.
        plt.tight_layout()
        plt.savefig('stacked_all_seq.png', dpi=1200)
        plt.savefig('stacked_all_seq.svg', dpi=1200)
                
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
