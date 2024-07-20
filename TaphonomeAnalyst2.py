import os
import copy
import venn
import random
import argparse
import warnings
import matplotlib
import numpy as np
import pandas as pd
import igraph as ig
import networkx as nx
import seaborn as sns
from scipy import stats
from pathlib import Path
from matplotlib import cm
from functools import reduce
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import community as community_louvain
from itertools import combinations, permutations
from skbio.diversity.alpha import chao1, chao1_ci
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import braycurtis, euclidean
from sparcc import SparCC_MicNet

warnings.filterwarnings("ignore")

sns.reset_orig()
sns.set_style('white')
matplotlib.use('Agg')
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['lines.linewidth'] = 2


# plt.rcParams['figure.figsize'] = (24, 16)

def str2list(s):
    return s.split(',')


def str2dictlist(s):
    group_dict = {}
    group_list = s.split(',')
    for i in group_list:
        k, v = i.split(':')
        group_dict[k] = v.split('/')
    return group_dict


# M5
# 水生
def clusterenv(args):
    all_level_set = set()
    for i in data_df_dict.keys():
        all_level_set |= set(data_df_dict[i][args.level])
    all_level_set.remove('unknown')
    all_otu_df_list = list()
    for i in data_df_dict.keys():
        level_list = [m for m, n in data_df_dict[i][args.level].value_counts().items() if m != 'unknown']
        count_list = [n for m, n in data_df_dict[i][args.level].value_counts().items() if m != 'unknown']
        diff_set = all_level_set - set(level_list)
        level_list.extend(list(diff_set))
        count_list.extend([0] * len(diff_set))
        all_otu_df_list.append(pd.DataFrame(count_list, index=level_list))
    otu_level_df = pd.concat([i.T for i in all_otu_df_list], ignore_index=True)
    otu_level_df.index = data_df_dict.keys()
    del_columns_list = list()
    for i in otu_level_df.columns:
        if all([j < 5 for j in otu_level_df[i]]):
            del_columns_list.append(i)
    otu_level_df_filter = otu_level_df.drop(columns=del_columns_list)
    if args.aquatic:
        otu_level_df_filter = otu_level_df_filter.loc[:,
                              otu_level_df_filter.columns.isin(args.aquatic)]
    otu_level_df_filter_standardscale = otu_level_df_filter.apply(
        lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
    plt.figure(figsize=otu_level_df_filter.shape, dpi=200)
    cluster_map = sns.clustermap(otu_level_df_filter_standardscale.T, cmap='Blues', row_cluster=False,
                                 method='average',
                                 metric='braycurtis')
    col_order = cluster_map.dendrogram_col.reordered_ind
    x_axis_labels = cluster_map.ax_heatmap.get_xticklabels()
    plt.savefig(fname=f'{args.output}_aquatic.{args.format}', bbox_inches='tight')
    plt.clf()
    otu_level_df_filter_standardscale_sort = otu_level_df_filter_standardscale.iloc[col_order, :]
    mergings = linkage(otu_level_df_filter_standardscale_sort, method='average', metric='braycurtis',
                       optimal_ordering=False)
    plt.figure(figsize=otu_level_df_filter.shape, dpi=200)
    dendrogram(Z=mergings, labels=otu_level_df_filter_standardscale_sort.index, leaf_rotation=90)
    plt.savefig(fname=f'{args.output}_tree.{args.format}', bbox_inches='tight')
    x_labels = [label.get_text() for label in x_axis_labels]

    if args.geochem:
        geochem = pd.read_excel(f'{args.geochem}', index_col=[0])
        geochem.dropna(axis=1, how='all', inplace=True)
        geochem.sort_index(axis=1, inplace=True)
        geochem = geochem.reindex(labels=x_labels)
        geochem_standardscale = geochem.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
        plt.figure(figsize=geochem.shape, dpi=200)
        ax = sns.heatmap(geochem_standardscale.T, cmap='Blues', annot=False, center=0, linewidths=.5, fmt='.2g')
        ax.set_xlabel('')
        plt.savefig(fname=f'{args.output}_geochem.{args.format}', bbox_inches='tight')
    print('Finished!')


# M4
# 自定义分组
def divvenn(args):
    groups_list = [i.split('/') for i in args.groups]
    all_level_set = set()
    for i in data_df_dict.keys():
        all_level_set |= set(data_df_dict[i][args.level])
    all_level_set.remove('unknown')
    all_otu_df_list = list()
    for i in data_df_dict.keys():
        level_list = [m for m, n in data_df_dict[i][args.level].value_counts().items() if m != 'unknown']
        count_list = [n for m, n in data_df_dict[i][args.level].value_counts().items() if m != 'unknown']
        diff_set = all_level_set - set(level_list)
        level_list.extend(list(diff_set))
        count_list.extend([0] * len(diff_set))
        all_otu_df_list.append(pd.DataFrame(count_list, index=level_list))
    otu_level_df = pd.concat([i.T for i in all_otu_df_list], ignore_index=True)
    otu_level_df.index = data_df_dict.keys()
    labels_data_list = []
    for i in groups_list:
        tmp = otu_level_df.loc[i, :].sum()
        labels_data_list.append(set(tmp[tmp != 0].index))
    data = {i: j for i, j in zip(args.groups, labels_data_list)}
    plt.figure(figsize=(20, 20), dpi=200)
    if len(groups_list) == 6:
        venn.pseudovenn(data, fmt='{size} ({percentage:.1f}%)', legend_loc="upper center")
    else:
        venn.venn(data, fmt='{size} ({percentage:.1f}%)', legend_loc="upper center")
    plt.gca().get_legend().set_bbox_to_anchor((0.5, 1.1))
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


# M3
def TGotus(args):
    plot_df_dict = {
        i: data_df_dict[i][args.level].value_counts().rename_axis(args.level).reset_index(name=i).set_index(
            args.level)
        for i in data_df_dict}
    otu_level_df_T = reduce(
        lambda left, right: pd.merge(left, right, how='outer', left_index=True, right_index=True),
        list(plot_df_dict.values()))
    otu_level_df_T['total'] = otu_level_df_T.apply(lambda x: x.sum(), axis=1)
    level_count_df = otu_level_df_T[['total']].reset_index()

    level_grade_df = pd.DataFrame(columns=grade_list)
    for i in level_count_df[args.level]:
        grade_A2E_count_list = list()
        tmp_grade_tuple_list = list()
        for j in data_df_dict.keys():
            tmp_grade_tuple_list.extend([(m, n) for m, n in data_df_dict[j][data_df_dict[j][args.level] == i][
                'taphonomic grade'].value_counts().items()])
        for k in grade_list:
            try:
                grade_A2E_count_list.append(sum([j[1] for j in tmp_grade_tuple_list if j[0] == k]))
            except Exception as e:
                grade_A2E_count_list.append(0)
        level_grade_df.loc[level_count_df[args.level].tolist().index(i)] = grade_A2E_count_list
    level_grade_count_df = pd.concat([level_count_df, level_grade_df], axis=1).set_index(args.level)
    for i in grade_list:
        level_grade_count_df[i] /= level_grade_count_df['total']
    level_grade_count_df.sort_values(by=['A'], inplace=True)
    level_grade_count_df_ = level_grade_count_df.drop('total', axis=1)
    fig = level_grade_count_df_.plot.barh(stacked=True, color=color_list[:len(grade_list)])
    fig.set_ylabel('')
    ax2 = fig.twinx()
    ax2.set_ylim(fig.get_ylim())
    ax2.set_yticks([i for i in range(level_grade_count_df.shape[0])])
    ax2.set_yticklabels([' n = ' + str(int(i)) for i in level_grade_count_df['total']])
    fig.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), fontsize=12)
    fig.tick_params(labelsize=12)
    ax2.tick_params(labelsize=12)
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


def TGplots(args):
    all_plot_grade_count_df = pd.DataFrame(index=data_df_dict.keys(), columns=grade_list)
    for i in data_df_dict.keys():
        grade_and_count_list = [(m, n) for m, n in
                                data_df_dict[i]['taphonomic grade'].value_counts().sort_index().items()]
        tmp_grade_list = list(map(itemgetter(0), grade_and_count_list))
        tmp_grade_count_list = list(map(itemgetter(1), grade_and_count_list))
        if len(tmp_grade_list) != 5:
            tmp_grade_count_list = [0 if j not in tmp_grade_list else tmp_grade_count_list[tmp_grade_list.index(j)]
                                    for j in grade_list]
        all_plot_grade_count_df.loc[i] = tmp_grade_count_list

    all_plot_grade_count_df['total'] = all_plot_grade_count_df.apply(lambda x: x.sum(), axis=1)
    all_plot_grade_count_df.fillna(0, inplace=True)
    all_plot_grade_count_df['A+B'] = all_plot_grade_count_df['A'] + all_plot_grade_count_df['B']
    for i in grade_list:
        all_plot_grade_count_df[i] /= all_plot_grade_count_df['total']
    all_plot_grade_count_df['A+B'] /= all_plot_grade_count_df['total']
    if args.groups:
        groups_list = [i.split('/') for i in args.groups]
        all_plot_grade_count_df_sort = pd.concat(
            [all_plot_grade_count_df.loc[i, :].sort_values(by=['A+B'], ascending=True) for i in groups_list[::-1]])
    else:
        all_plot_grade_count_df_sort = all_plot_grade_count_df.sort_values(by=['A+B'], ascending=True)
    all_plot_grade_count_df_sort_ = all_plot_grade_count_df_sort.drop(['total', 'A+B'], axis=1)
    fig = all_plot_grade_count_df_sort_.plot.barh(stacked=True, color=color_list[:len(grade_list)])
    ax2 = fig.twinx()
    ax2.set_ylim(fig.get_ylim())
    ax2.set_yticks([i for i in range(all_plot_grade_count_df_sort.shape[0])])
    ax2.set_yticklabels([' n = ' + str(i) for i in all_plot_grade_count_df_sort['total']])
    fig.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), fontsize=12)
    fig.tick_params(labelsize=12)
    ax2.tick_params(labelsize=12)
    if args.groups:
        plt.savefig(fname=f'{args.output}_groups.{args.format}', bbox_inches='tight')
    else:
        plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


# M2
def abundplots(args):
    all_plot_grade_count_df = pd.DataFrame(index=data_df_dict.keys(), columns=grade_list)
    for i in data_df_dict.keys():
        grade_and_count_list = [(m, n) for m, n in
                                data_df_dict[i]['taphonomic grade'].value_counts().sort_index().items()]
        tmp_grade_list = list(map(itemgetter(0), grade_and_count_list))
        tmp_grade_count_list = list(map(itemgetter(1), grade_and_count_list))
        if len(tmp_grade_list) != 5:
            tmp_grade_count_list = [0 if j not in tmp_grade_list else tmp_grade_count_list[tmp_grade_list.index(j)]
                                    for j in grade_list]
        all_plot_grade_count_df.loc[i] = tmp_grade_count_list
    all_plot_grade_count_df['total'] = all_plot_grade_count_df.apply(lambda x: x.sum(), axis=1)
    for i in grade_list:
        all_plot_grade_count_df[i] /= all_plot_grade_count_df['total']
    plot_normalize_df_dict = {
        i: data_df_dict[i][args.level].value_counts(normalize=True).rename_axis(args.level).reset_index(
            name=i).set_index(args.level)
        for i in data_df_dict}
    otu_level_normalize_df = reduce(
        lambda left, right: pd.merge(left, right, how='outer', left_index=True, right_index=True),
        list(plot_normalize_df_dict.values())).T
    fig = otu_level_normalize_df.plot.barh(stacked=True, color=color_list[:len(otu_level_normalize_df.columns)])
    ax2 = fig.twinx()
    ax2.set_ylim(fig.get_ylim())
    ax2.set_yticks([i for i in range(otu_level_normalize_df.shape[0])])
    ax2.set_yticklabels([' n = ' + str(int(i)) for i in all_plot_grade_count_df['total']])
    fig.legend(loc='center right', bbox_to_anchor=(1.35, 0.5), fontsize=12)
    fig.tick_params(labelsize=12)
    ax2.tick_params(labelsize=12)
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


def cooccurnet(args):
    plot_df_dict = {
        i: (data_df_dict[i][args.level].value_counts(normalize=True) * 3000).rename_axis(args.level).reset_index(
            name=i).set_index(args.level) for i in data_df_dict}
    otu_level_df = reduce(lambda left, right: pd.merge(left, right, how='outer', left_index=True, right_index=True),
                          list(plot_df_dict.values())).fillna(0).T
    try:
        del otu_level_df['unknown']
    except Exception as e:
        pass
    all_comb = list(combinations(list(otu_level_df), 2))
    limit_list = []
    nodes = []
    edges = []
    if args.corr in ['pearson', 'spearman', 'kendall']:
        corr_dict = {'pearson': 'pearsonr', 'spearman': 'spearmanr', 'kendall': 'kendalltau'}
        for comb in all_comb:
            corr_and_p_value = getattr(stats, corr_dict[args.corr])(otu_level_df[f'{comb[0]}'],
                                                                    otu_level_df[f'{comb[1]}'])
            if corr_and_p_value[1] < args.p_value and corr_and_p_value[0] > args.corr_coef:
                limit_list.append(comb[0])
                limit_list.append(comb[1])
        limit_list = list(set(limit_list))
        for i in limit_list:
            node = (limit_list.index(i), i)
            nodes.append(node)
        for comb in all_comb:
            corr_and_p_value = getattr(stats, corr_dict[args.corr])(otu_level_df[f'{comb[0]}'],
                                                                    otu_level_df[f'{comb[1]}'])
            if corr_and_p_value[1] < args.p_value and corr_and_p_value[0] > args.corr_coef:
                edge = (limit_list.index(comb[0]), limit_list.index(comb[1]), corr_and_p_value[0])
                edges.append(edge)
    elif args.corr == 'sparcc':
        # set parameters for SparCC
        n_iteractions = 30
        x_iteractions = 30
        low_abundance = False
        threshold = 0.1
        normalization = 'dirichlet'
        log_transform = True
        num_simulate_data = 100
        type_pvalues = 'two_sided'

        # Create Sparcc object
        SparCC_MN = SparCC_MicNet(n_iteractions=n_iteractions,
                                  x_iteractions=x_iteractions,
                                  low_abundance=low_abundance,
                                  threshold=threshold,
                                  normalization=normalization,
                                  log_transform=log_transform,
                                  num_simulate_data=num_simulate_data,
                                  type_pvalues=type_pvalues,
                                  )
        data = otu_level_df.loc[:, (otu_level_df != 0).any(axis=0)].T
        SparCC_MN.run_all(data_input=data)
        DF_SparCC = pd.read_csv(Path(SparCC_MN.save_corr_file).resolve(), index_col=0)
        DF_PValues = pd.read_csv(Path(SparCC_MN.outfile_pvals).resolve(), index_col=0)
        for i in combinations([i for i in range(len(DF_SparCC))], 2):
            if DF_PValues.iloc[i[0], i[1]] < args.p_value and DF_SparCC.iloc[i[0], i[1]] > args.corr_coef:
                limit_list.append(i[0])
                limit_list.append(i[1])
        limit_list = list(set(limit_list))
        for i in limit_list:
            node = (limit_list.index(i), data.index[i])
            nodes.append(node)
        for i in combinations([i for i in range(len(DF_SparCC))], 2):
            if DF_PValues.iloc[i[0], i[1]] < args.p_value and DF_SparCC.iloc[i[0], i[1]] > args.corr_coef:
                edge = (limit_list.index(i[0]), limit_list.index(i[1]), DF_SparCC.iloc[i[0], i[1]])
                edges.append(edge)
    with open('./cooccurnet.gml', 'w') as f:
        f.write("graph\n")
        f.write("[\n")
        for i in nodes:
            f.write("node\n")
            f.write("[\n")
            f.write("id " + str(i[0]) + "\n")
            f.write("label \"" + i[1] + "\"\n")
            f.write("]\n")
        for i in edges:
            f.write("edge\n")
            f.write("[\n")
            f.write("source " + str(i[0]) + "\n")
            f.write("target " + str(i[1]) + "\n")
            f.write("value " + str(i[2]) + "\n")
            f.write("]\n")
        f.write("]")
    f.close()
    G = nx.read_gml('./cooccurnet.gml')
    partition = community_louvain.best_partition(G)
    pos = nx.kamada_kawai_layout(G)
    fig, ax = plt.subplots(figsize=(40, 24), frameon=False)
    for key, spine in ax.spines.items():
        if key == 'right' or key == 'top' or key == 'left' or key == 'bottom':
            spine.set_visible(False)
    cmap = cm.get_cmap(None, max(partition.values()) + 1)
    node_color_list = []
    for i in list(partition.values()):
        node_color_list.append(color_list[i])
    nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=400,
                           cmap=cmap, node_color=node_color_list, alpha=0.8)
    nx.draw_networkx_edges(G, pos, alpha=0.2)
    nx.draw_networkx_labels(G, pos, alpha=0.8, font_size=20, verticalalignment='bottom')
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')
    plt.clf()
    network_dict = dict()
    for i in sorted(set(partition.values())):
        network_list = []
        for j in zip(partition.values(), partition.keys()):
            if j[0] == i:
                network_list.append(j[1])
        network_dict[i] = network_list
    network_dict_ = copy.deepcopy(network_dict)
    for i in network_dict.keys():
        for j in network_dict[i]:
            for m in list(set(partition.keys()) - set(network_dict[i])):
                if G.has_edge(j, m):
                    network_dict_[i].append(m)
    network_dict_ = {i: set(j) for i, j in network_dict_.items()}
    color_dict_venn = {i: j for i, j in zip(partition.values(), node_color_list)}
    color_list_venn = [j for i, j in sorted(color_dict_venn.items())]
    plt.figure(figsize=(20, 20), dpi=200)
    cmap = [mcolors.hex2color(i) for i in color_list_venn]
    if len(network_dict_) == 6:
        venn.pseudovenn(network_dict_, cmap=cmap, fmt='{size}', legend_loc="best")
    else:
        venn.venn(network_dict_, cmap=cmap, fmt='{size}', legend_loc="best")
    # plt.gca().get_legend().set_bbox_to_anchor((0.5, 1.1))
    plt.savefig(fname=f'{args.output}_venn.{args.format}', bbox_inches='tight')


# M1 chao1
def chao(args):
    plot_dict = args.groups
    otu_df_dict = dict()
    for i in data_df_dict.keys():
        result_dict = dict()
        for j in range(50, data_df_dict[i].shape[0] + 1, 50):
            result_list = [k for k in data_df_dict[i][:j][args.level] if not k in [np.nan, 'unknown']]
            result_dict[j] = [(k, result_list.count(k)) for k in set(result_list)]
        otu_df_dict[i] = result_dict
    plot_otu_df_dict = dict()
    for i in plot_dict.keys():
        plot_df_list = list()
        for j in plot_dict[i]:
            plot_df_list.extend([k[0] for k in otu_df_dict[j][max(otu_df_dict[j].keys())]])
        plot_otu_df_dict[i] = set(plot_df_list)
    otu_df_dict_ = copy.deepcopy(otu_df_dict)
    for i in plot_dict.keys():
        for j in otu_df_dict.keys():
            if j in plot_dict[i]:
                for k in range(50, max(otu_df_dict[j].keys()) + 1, 50):
                    diff = set(plot_otu_df_dict[i]).difference(set([l[0] for l in otu_df_dict[j][k]]))
                    zip_diff = zip(list(diff), [0] * len(diff))
                    otu_df_dict_[j][k].extend([l for l in zip_diff])
    step_df_list = []
    for i in plot_dict.keys():
        for j in plot_dict[i]:
            for k in range(50, max(otu_df_dict[j].keys()) + 1, 50):
                columns = [l[0] for l in otu_df_dict_[j][k]]
                line = [l[1] for l in otu_df_dict_[j][k]]
                df = pd.DataFrame(line, columns)
                step_df_list.append(df)
    step_idx = []
    step_len_list = [i.shape[0] for i in step_df_list]
    step_len_set_list = list(set(step_len_list))
    step_len_set_list.sort(key=step_len_list.index)
    for i in step_len_set_list:
        step_idx.append([k for k, l in enumerate([j.shape[0] for j in step_df_list]) if l == i][-1])
    step_dict = dict()
    count = 0
    for i in plot_dict.keys():
        step_dict[i] = [i for i in zip([0] + step_idx[:-1], step_idx)][count]
        count += 1
    step_dict_ = dict()
    base = 0
    for i in step_dict.keys():
        step_dict_[i] = dict()
        count_ = 0
        plot_list = plot_dict[i]
        step_start, step_end = step_dict[i][0], step_dict[i][1]
        plot_len = len(plot_list)
        interval = (step_end - step_start + 1) // plot_len
        for j in plot_list:
            if base == 0:
                start = step_start + count_ * interval
            else:
                start = step_start + count_ * interval + 1
            end = start + interval - 1
            step_dict_[i][j] = (start, end)
            count_ += 1
        base += 1
    step_df_concat = pd.concat([i.T for i in step_df_list], ignore_index=True).fillna(0)
    plot_chao_dict = dict()
    plot_chao_ci_dict = dict()
    for i in step_dict_.keys():
        for j in step_dict_[i]:
            plot_chao_dict[j] = list()
            plot_chao_ci_dict[j] = list()
            for k in range(step_dict_[i][j][0], step_dict_[i][j][1] + 1):
                plot_chao_dict[j].append(chao1(step_df_concat.T[k]))
                plot_chao_ci_dict[j].append(chao1_ci(step_df_concat.T[k]))
    ncols = 2
    nrows = (len(plot_dict.keys()) + 1) // 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, figsize=(nrows * 5, ncols * 15))
    ax_list = [axes[i, j] for i in range(nrows) for j in range(ncols)]
    loc_list = [i for i in step_dict_.keys()]
    sns.set_style('ticks')
    sns.set_palette(palette=color_list[:len(loc_list)], n_colors=len(loc_list), desat=None, color_codes=True)
    for i in loc_list:
        plot_list = [n for n in plot_chao_dict.keys() if n in plot_dict[i]]
        for j in plot_list:
            sns.regplot([n for n in range(50, len(plot_chao_dict[j]) * 50 + 1, 50)], plot_chao_dict[j], logx=True,
                        label=j, ax=ax_list[loc_list.index(i)], scatter_kws={'s': 12})
            ci_0 = np.array([n[0] for n in plot_chao_ci_dict[j]]) - np.array(plot_chao_dict[j])
            ci_0 = [abs(n) for n in ci_0]
            ci_1 = np.array([n[1] for n in plot_chao_ci_dict[j]]) - np.array(plot_chao_dict[j])
            ci_1 = [abs(n) for n in ci_1]
            ci = np.array([ci_0, ci_1])
            ax_list[loc_list.index(i)].errorbar([n for n in range(50, len(plot_chao_dict[j]) * 50 + 1, 50)],
                                                plot_chao_dict[j], yerr=ci, fmt='none')
            ax_list[loc_list.index(i)].tick_params(labelsize=16)
            ax_list[loc_list.index(i)].xaxis.set_tick_params(labelbottom=True)
            ax_list[loc_list.index(i)].set_xlabel('Number of fossil individuals', fontsize=18)
            ax_list[loc_list.index(i)].set_ylabel(f'Number of {args.level}', fontsize=18)
            ax_list[loc_list.index(i)].legend(markerscale=3, fontsize=18)
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


# M1 sobs
def samplecurve(args):
    plot_dict = args.groups
    otu_df_dict = dict()
    for i in data_df_dict.keys():
        result_dict = dict()
        for j in range(50, data_df_dict[i].shape[0] + 1, 50):
            result_list = [k for k in data_df_dict[i][:j][args.level] if k not in [np.nan, 'unknown']]
            result_dict[j] = [(k, result_list.count(k)) for k in set(result_list)]
        otu_df_dict[i] = result_dict
    plot_otu_df_dict = dict()
    for i in plot_dict.keys():
        plot_df_list = list()
        for j in plot_dict[i]:
            plot_df_list.extend([k[0] for k in otu_df_dict[j][list(otu_df_dict[j].keys())[-1]]])
        plot_otu_df_dict[i] = set(plot_df_list)
    otu_df_dict_ = copy.deepcopy(otu_df_dict)
    for i in plot_dict.keys():
        for j in otu_df_dict.keys():
            if j in plot_dict[i]:
                for k in range(50, list(otu_df_dict[j].keys())[-1] + 1, 50):
                    diff = set(plot_otu_df_dict[i]).difference(set([l[0] for l in otu_df_dict[j][k]]))
                    zip_diff = zip(list(diff), [0] * len(diff))
                    otu_df_dict_[j][k].extend([l for l in zip_diff])
    step_df_list = []
    for i in plot_dict.keys():
        for j in plot_dict[i]:
            for k in range(50, list(otu_df_dict[j].keys())[-1] + 1, 50):
                columns = [l[0] for l in otu_df_dict_[j][k]]
                line = [l[1] for l in otu_df_dict_[j][k]]
                df = pd.DataFrame(line, columns)
                step_df_list.append(df)
    step_idx = []
    step_len_list = [i.shape[0] for i in step_df_list]
    step_len_set_list = list(set(step_len_list))
    step_len_set_list.sort(key=step_len_list.index)
    for i in step_len_set_list:
        step_idx.append([k for k, l in enumerate([j.shape[0] for j in step_df_list]) if l == i][-1])
    step_dict = dict()
    count = 0
    for i in plot_dict.keys():
        step_dict[i] = [i for i in zip([0] + step_idx[:-1], step_idx)][count]
        count += 1
    step_dict_ = dict()
    base = 0
    for i in step_dict.keys():
        step_dict_[i] = dict()
        count_ = 0
        plot_list = plot_dict[i]
        step_start, step_end = step_dict[i][0], step_dict[i][1]
        plot_len = len(plot_list)
        interval = (step_end - step_start + 1) // plot_len
        for j in plot_list:
            if base == 0:
                start = step_start + count_ * interval
            else:
                start = step_start + count_ * interval + 1
            end = start + interval - 1
            step_dict_[i][j] = (start, end)
            count_ += 1
        base += 1
    step_df_concat = pd.concat([i.T for i in step_df_list], ignore_index=True).fillna(0)
    plot_real_dict = dict()
    for i in step_dict_.keys():
        for j in step_dict_[i]:
            plot_real_dict[j] = list()
            for k in range(step_dict_[i][j][0], step_dict_[i][j][1] + 1):
                plot_real_dict[j].append(len([l for l in step_df_concat.T[k] if int(l) != 0]))
    ncols = 2
    nrows = (len(plot_dict.keys()) + 1) // 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, figsize=(nrows * 5, ncols * 15))
    ax_list = [axes[i, j] for i in range(nrows) for j in range(ncols)]
    loc_list = [i for i in step_dict_.keys()]
    sns.set_style('ticks')
    sns.set_palette(palette=color_list[:len(loc_list)], n_colors=len(loc_list), desat=None, color_codes=True)
    for i in loc_list:
        plot_list = [n for n in plot_real_dict.keys() if n in plot_dict[i]]
        for j in plot_list:
            sns.regplot([n for n in range(50, len(plot_real_dict[j]) * 50 + 1, 50)], plot_real_dict[j], logx=True,
                        label=j, ax=ax_list[loc_list.index(i)], scatter_kws={'s': 12})
            ax_list[loc_list.index(i)].errorbar([n for n in range(50, len(plot_real_dict[j]) * 50 + 1, 50)],
                                                plot_real_dict[j], fmt='none')
            ax_list[loc_list.index(i)].tick_params(labelsize=16)
            ax_list[loc_list.index(i)].xaxis.set_tick_params(labelbottom=True)
            ax_list[loc_list.index(i)].set_xlabel('Number of fossil individuals', fontsize=18)
            ax_list[loc_list.index(i)].set_ylabel(f'Number of {args.level}', fontsize=18)
            ax_list[loc_list.index(i)].legend(markerscale=3, fontsize=18)
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


def _otus_rare(freq_counts, rare_threshold):
    return freq_counts[1:rare_threshold + 1].sum()


def _otus_abundant(freq_counts, rare_threshold):
    return freq_counts[rare_threshold + 1:].sum()


def _number_rare(freq_counts, rare_threshold, gamma=False):
    n_rare = 0
    if gamma:
        for i, j in enumerate(freq_counts[:rare_threshold + 1]):
            n_rare = n_rare + (i * j) * (i - 1)
    else:
        for i, j in enumerate(freq_counts[:rare_threshold + 1]):
            n_rare = n_rare + (i * j)
    return n_rare


def _ace(counts, rare_threshold=10):
    counts = np.asarray(counts)
    freq_counts = np.bincount(counts)
    s_rare = _otus_rare(freq_counts, rare_threshold)
    singles = freq_counts[1]
    s_abun = _otus_abundant(freq_counts, rare_threshold)
    if s_rare == 0:
        return s_abun
    n_rare = _number_rare(freq_counts, rare_threshold)
    c_ace = 1 - singles / n_rare
    top = s_rare * _number_rare(freq_counts, rare_threshold, gamma=True)
    bottom = c_ace * n_rare * (n_rare - 1)
    gamma_ace = (top / bottom) - 1
    if gamma_ace < 0:
        gamma_ace = 0
    return s_abun + (s_rare / c_ace) + ((singles / c_ace) * gamma_ace)


# M1 ace
def ace(args):
    plot_dict = args.groups
    otu_df_dict = dict()
    for i in data_df_dict.keys():
        result_dict = dict()
        for j in range(50, data_df_dict[i].shape[0] + 1, 50):
            result_list = [k for k in data_df_dict[i][:j][args.level] if not k in [np.nan, 'unknown']]
            result_dict[j] = [(k, result_list.count(k)) for k in set(result_list)]
        otu_df_dict[i] = result_dict
    plot_otu_df_dict = dict()
    for i in plot_dict.keys():
        plot_df_list = list()
        for j in plot_dict[i]:
            plot_df_list.extend([k[0] for k in otu_df_dict[j][max(otu_df_dict[j].keys())]])
        plot_otu_df_dict[i] = set(plot_df_list)
    otu_df_dict_ = copy.deepcopy(otu_df_dict)
    for i in plot_dict.keys():
        for j in otu_df_dict.keys():
            if j in plot_dict[i]:
                for k in range(50, max(otu_df_dict[j].keys()) + 1, 50):
                    diff = set(plot_otu_df_dict[i]).difference(set([l[0] for l in otu_df_dict[j][k]]))
                    zip_diff = zip(list(diff), [0] * len(diff))
                    otu_df_dict_[j][k].extend([l for l in zip_diff])
    step_df_list = []
    for i in plot_dict.keys():
        for j in plot_dict[i]:
            for k in range(50, max(otu_df_dict[j].keys()) + 1, 50):
                columns = [l[0] for l in otu_df_dict_[j][k]]
                line = [l[1] for l in otu_df_dict_[j][k]]
                df = pd.DataFrame(line, columns)
                step_df_list.append(df)
    step_idx = []
    step_len_list = [i.shape[0] for i in step_df_list]
    step_len_set_list = list(set(step_len_list))
    step_len_set_list.sort(key=step_len_list.index)
    for i in step_len_set_list:
        step_idx.append([k for k, l in enumerate([j.shape[0] for j in step_df_list]) if l == i][-1])
    step_dict = dict()
    count = 0
    for i in plot_dict.keys():
        step_dict[i] = [i for i in zip([0] + step_idx[:-1], step_idx)][count]
        count += 1
    step_dict_ = dict()
    base = 0
    for i in step_dict.keys():
        step_dict_[i] = dict()
        count_ = 0
        plot_list = plot_dict[i]
        step_start, step_end = step_dict[i][0], step_dict[i][1]
        plot_len = len(plot_list)
        interval = (step_end - step_start + 1) // plot_len
        for j in plot_list:
            if base == 0:
                start = step_start + count_ * interval
            else:
                start = step_start + count_ * interval + 1
            end = start + interval - 1
            step_dict_[i][j] = (start, end)
            count_ += 1
        base += 1
    step_df_concat = pd.concat([i.T for i in step_df_list], ignore_index=True).fillna(0)
    plot_ace_dict = dict()
    for i in step_dict_.keys():
        for j in step_dict_[i]:
            plot_ace_dict[j] = list()
            for k in range(step_dict_[i][j][0], step_dict_[i][j][1] + 1):
                plot_ace_dict[j].append(_ace(step_df_concat.T[k].astype(np.int64), rare_threshold=10))
    ncols = 2
    nrows = (len(plot_dict.keys()) + 1) // 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, figsize=(nrows * 5, ncols * 15))
    ax_list = [axes[i, j] for i in range(nrows) for j in range(ncols)]
    loc_list = [i for i in step_dict_.keys()]
    sns.set_style('ticks')
    sns.set_palette(palette=color_list[:len(loc_list)], n_colors=len(loc_list), desat=None, color_codes=True)
    for i in loc_list:
        plot_list = [n for n in plot_ace_dict.keys() if n in plot_dict[i]]
        for j in plot_list:
            sns.regplot([n for n in range(50, len(plot_ace_dict[j]) * 50 + 1, 50)], plot_ace_dict[j], logx=True,
                        label=j, ax=ax_list[loc_list.index(i)], scatter_kws={'s': 12})
            ax_list[loc_list.index(i)].errorbar([n for n in range(50, len(plot_ace_dict[j]) * 50 + 1, 50)],
                                                plot_ace_dict[j], fmt='none')
            ax_list[loc_list.index(i)].tick_params(labelsize=16)
            ax_list[loc_list.index(i)].xaxis.set_tick_params(labelbottom=True)
            ax_list[loc_list.index(i)].set_xlabel('Number of fossil individuals', fontsize=18)
            ax_list[loc_list.index(i)].set_ylabel(f'Number of {args.level}', fontsize=18)
            ax_list[loc_list.index(i)].legend(markerscale=3, fontsize=18)
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


# M8
def corrotus(args):
    plot_df_dict = {
        i: (data_df_dict[i][args.level].value_counts(normalize=True) * 3000).rename_axis(args.level).reset_index(
            name=i).set_index(args.level) for i in data_df_dict}
    otu_level_df = reduce(lambda left, right: pd.merge(left, right, how='outer', left_index=True, right_index=True),
                          list(plot_df_dict.values())).fillna(0).T
    try:
        del otu_level_df['unknown']
    except Exception as e:
        pass
    all_perm = list(permutations(list(otu_level_df), 2))

    name_list = list(otu_level_df.columns)
    data = pd.DataFrame(columns=name_list, index=name_list).sort_index().sort_index(axis=1)
    for perm in all_perm:
        if args.corr in ['pearson', 'spearman', 'kendall']:
            corr_dict = {'pearson': 'pearsonr', 'spearman': 'spearmanr', 'kendall': 'kendalltau'}
            corr_and_p_value = getattr(stats, corr_dict[args.corr])(otu_level_df[f'{perm[0]}'],
                                                                    otu_level_df[f'{perm[1]}'])
            if corr_and_p_value[1] < args.p_value:
                data.loc[perm[0], perm[1]] = corr_and_p_value[0]
        elif args.corr == 'sparcc':
            filter = otu_level_df.loc[:, (otu_level_df != 0).any(axis=0)]
            n_iteractions = 30
            x_iteractions = 30
            low_abundance = False
            threshold = 0.1
            normalization = 'dirichlet'
            log_transform = True
            num_simulate_data = 100
            type_pvalues = 'two_sided'
            SparCC_MN = SparCC_MicNet(n_iteractions=n_iteractions,
                                      x_iteractions=x_iteractions,
                                      low_abundance=low_abundance,
                                      threshold=threshold,
                                      normalization=normalization,
                                      log_transform=log_transform,
                                      num_simulate_data=num_simulate_data,
                                      type_pvalues=type_pvalues,
                                      )
            SparCC_MN.run_all(data_input=filter.T)
            corr = pd.read_csv(Path(SparCC_MN.save_corr_file).resolve(), index_col=0)
            corr.index = filter.columns
            corr.columns = filter.columns
            p = pd.read_csv(Path(SparCC_MN.outfile_pvals).resolve(), index_col=0)
            p.index = filter.columns
            p.columns = filter.columns
            for col1 in name_list:
                for col2 in name_list:
                    if col1 != col2:
                        try:
                            _p = p.at[col1, col2]
                            _corr = corr.at[col1, col2]
                            if p > args.p_value or abs(_corr) < args.corr_coef:
                                data.at[col1, col2] = 0
                            else:
                                data.at[col1, col2] = corr.at[col1, col2]
                        except Exception as e:
                            data.at[col1, col2] = np.NaN
    data.fillna(np.nan, inplace=True)
    mask = np.zeros_like(data, dtype=np.bool_)
    mask[np.triu_indices_from(mask)] = True
    fig, ax = plt.subplots(figsize=data.shape)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=-360)
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.heatmap(data, cmap=cmap, annot=True, mask=mask, center=0, linewidths=.5, fmt='.2g')
    data.to_csv(f'{args.output}.csv')
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


def mantel(args):
    os.environ['R_HOME'] = args.rhome
    from rpy2 import robjects
    from rpy2.robjects.vectors import IntVector, StrVector
    pd.read_excel(args.geochem, header=0, index_col=[0]).sort_index().to_csv('./geochem.csv')
    aquatic_level_list = args.aquatic
    aquatic_tmp_list = []
    for i in data_df_dict.keys():
        for j in aquatic_level_list:
            aquatic_tmp_list.extend(
                data_df_dict[i][data_df_dict[i][args.level_aquatic] == j][args.level_terrestrial].tolist())
    aquatic_tmp_list = list(set(aquatic_tmp_list))
    terrestrial_level_list = []
    for i in data_df_dict.keys():
        for j in aquatic_tmp_list:
            terrestrial_level_list.extend(
                data_df_dict[i][data_df_dict[i][args.level_terrestrial] != j][args.level_terrestrial].tolist())
    terrestrial_level_list = [i for i in set(terrestrial_level_list) if i != 'unknown']
    all_level_set = set()
    for i in data_df_dict.keys():
        all_level_set |= set(data_df_dict[i][args.level_aquatic])
    all_level_set.remove('unknown')
    all_otu_df_list = list()
    for i in data_df_dict.keys():
        level_list = [m for m, n in data_df_dict[i][args.level_aquatic].value_counts().items() if m != 'unknown']
        count_list = [n for m, n in data_df_dict[i][args.level_aquatic].value_counts().items() if m != 'unknown']
        diff_set = all_level_set - set(level_list)
        level_list.extend(list(diff_set))
        count_list.extend([0] * len(diff_set))
        all_otu_df_list.append(pd.DataFrame(count_list, index=level_list))
    otu_level_df = pd.concat([i.T for i in all_otu_df_list], ignore_index=True)
    otu_level_df.index = data_df_dict.keys()
    del_columns_list = list()
    for i in otu_level_df.columns:
        if all([j < 5 for j in otu_level_df[i]]):
            del_columns_list.append(i)
    otu_level_df_filter = otu_level_df.drop(columns=del_columns_list)
    otu_level_df_filter = otu_level_df_filter.loc[:, otu_level_df_filter.columns.isin(aquatic_level_list)]
    otu_level_df_filter_standardscale_aquatic = otu_level_df_filter.apply(
        lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
    otu_level_df_filter_standardscale_aquatic.sort_index(inplace=True)
    all_level_set = set()
    for i in data_df_dict.keys():
        all_level_set |= set(data_df_dict[i][args.level_terrestrial])
    all_level_set.remove('unknown')
    all_otu_df_list = list()
    for i in data_df_dict.keys():
        level_list = [m for m, n in data_df_dict[i][args.level_terrestrial].value_counts().items() if m != 'unknown']
        count_list = [n for m, n in data_df_dict[i][args.level_terrestrial].value_counts().items() if m != 'unknown']
        diff_set = all_level_set - set(level_list)
        level_list.extend(list(diff_set))
        count_list.extend([0] * len(diff_set))
        all_otu_df_list.append(pd.DataFrame(count_list, index=level_list))
    otu_level_df = pd.concat([i.T for i in all_otu_df_list], ignore_index=True)
    otu_level_df.index = data_df_dict.keys()
    del_columns_list = list()
    for i in otu_level_df.columns:
        if all([j < 5 for j in otu_level_df[i]]):
            del_columns_list.append(i)
    otu_level_df_filter = otu_level_df.drop(columns=del_columns_list)
    otu_level_df_filter = otu_level_df_filter.loc[:, otu_level_df_filter.columns.isin(terrestrial_level_list)]
    otu_level_df_filter_standardscale_terrestrial = otu_level_df_filter.apply(
        lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
    otu_level_df_filter_standardscale_terrestrial.sort_index(inplace=True)
    aquatic_and_terrestrial = pd.concat(
        [otu_level_df_filter_standardscale_aquatic, otu_level_df_filter_standardscale_terrestrial], axis=1)
    aquatic_and_terrestrial.to_csv('./aquatic_and_terrestrial.csv')
    r_aquatic_end = IntVector([otu_level_df_filter_standardscale_aquatic.shape[1]])
    r_terrestrial_start = IntVector([otu_level_df_filter_standardscale_aquatic.shape[1] + 1])
    r_terrestrial_end = IntVector([aquatic_and_terrestrial.shape[1]])
    r_corr = StrVector([args.corr])
    r_title = StrVector([args.corr.capitalize()])
    r_output = StrVector([args.output])
    r_format = StrVector([args.format])
    robjects.globalenv['r_aquatic_end'] = r_aquatic_end
    robjects.globalenv['r_terrestrial_start'] = r_terrestrial_start
    robjects.globalenv['r_terrestrial_end'] = r_terrestrial_end
    robjects.globalenv['r_corr'] = r_corr
    robjects.globalenv['r_title'] = r_title
    robjects.globalenv['r_output'] = r_output
    robjects.globalenv['r_format'] = r_format
    robjects.r('Sys.setlocale("LC_CTYPE", "en_US.UTF-8")')
    robjects.r(f'''
        library(dplyr)
        library(linkET)  
        library(ggplot2)

        
        geochem <- read.csv('./geochem.csv', header = TRUE, row.names = 1, fileEncoding = "UTF-8")
        aquatic_and_terrestrial <- read.csv('./aquatic_and_terrestrial.csv', header = TRUE, row.names = 1, fileEncoding = "UTF-8")

        mantel <- mantel_test(aquatic_and_terrestrial, geochem,
                              spec_select = list(Aquatic = 1:r_aquatic_end, Terrestrial = r_terrestrial_start:r_terrestrial_end)) %>% 
          mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                          labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                 pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                          labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

        qcorrplot(correlate(geochem, method = r_corr), type = "lower", diag = FALSE) +
          geom_square() +
          geom_couple(aes(colour = pd, size = rd), 
                      data = mantel, 
                      curvature = nice_curvature()) +
          scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),  # 反转颜色顺序
                               limits = c(-1.0, 1.0),  # 修改颜色范围
                               breaks = seq(-1.0, 1.0, 0.5)) +  # 设置刻度
          scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
          scale_colour_manual(values = color_pal(3)) +
          guides(size = guide_legend(title = "Mantel's r",
                                     override.aes = list(colour = "grey35"), 
                                     order = 2),
                 colour = guide_legend(title = "Mantel's p", 
                                       override.aes = list(size = 3), 
                                       order = 1),
                 fill = guide_colorbar(title = paste0(r_title, "'s r"), order = 3)) 
        ggsave(file = paste0(r_output, ".", r_format), width = 12, height = 8)
    ''')


def dissenvtest(args):
    geochem_df = pd.read_excel(args.geochem, header=0, index_col=[0]).sort_index()
    environmental_distances = pd.DataFrame(index=geochem_df.index, columns=geochem_df.index)
    for i in geochem_df.index:
        for j in geochem_df.index:
            environmental_distances.loc[i, j] = euclidean(geochem_df.loc[i], geochem_df.loc[j])
    environmental_distances.sort_index(inplace=True)
    environmental_distances.sort_index(axis=1, inplace=True)
    result_matrix = np.full_like(environmental_distances, -1)
    for i in range(environmental_distances.shape[0]):
        for j in range(i + 1, environmental_distances.shape[1]):
            result_matrix[i, j] = environmental_distances.iloc[i, j]
    environmental_distances_half = pd.DataFrame(result_matrix, index=environmental_distances.index,
                                                columns=environmental_distances.columns)
    environmental_distances_list = [i for i in environmental_distances_half.values.flatten() if i != -1]
    aquatic_level_list = args.aquatic
    aquatic_tmp_list = []
    for i in data_df_dict.keys():
        for j in aquatic_level_list:
            aquatic_tmp_list.extend(
                data_df_dict[i][data_df_dict[i][args.level_aquatic] == j][args.level_terrestrial].tolist())
    aquatic_tmp_list = list(set(aquatic_tmp_list))
    terrestrial_level_list = []
    for i in data_df_dict.keys():
        for j in aquatic_tmp_list:
            terrestrial_level_list.extend(
                data_df_dict[i][data_df_dict[i][args.level_terrestrial] != j][args.level_terrestrial].tolist())
    terrestrial_level_list = [i for i in set(terrestrial_level_list) if i != 'unknown']
    all_level_set = set()
    for i in data_df_dict.keys():
        all_level_set |= set(data_df_dict[i][args.level_aquatic])
    all_level_set.remove('unknown')
    all_otu_df_list = list()
    for i in data_df_dict.keys():
        level_list = [m for m, n in data_df_dict[i][args.level_aquatic].value_counts().items() if m != 'unknown']
        count_list = [n for m, n in data_df_dict[i][args.level_aquatic].value_counts().items() if m != 'unknown']
        diff_set = all_level_set - set(level_list)
        level_list.extend(list(diff_set))
        count_list.extend([0] * len(diff_set))
        all_otu_df_list.append(pd.DataFrame(count_list, index=level_list))
    otu_level_df = pd.concat([i.T for i in all_otu_df_list], ignore_index=True)
    otu_level_df.index = data_df_dict.keys()
    del_columns_list = list()
    for i in otu_level_df.columns:
        if all([j < 5 for j in otu_level_df[i]]):
            del_columns_list.append(i)
    otu_level_df_filter = otu_level_df.drop(columns=del_columns_list)
    otu_level_df_filter = otu_level_df_filter.loc[:, otu_level_df_filter.columns.isin(aquatic_level_list)]
    otu_level_df_filter_standardscale_aquatic = otu_level_df_filter.apply(
        lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
    otu_level_df_filter_standardscale_aquatic.sort_index(inplace=True)
    bray_curtis_matrix_aquatic = pd.DataFrame(index=otu_level_df_filter_standardscale_aquatic.index,
                                              columns=otu_level_df_filter_standardscale_aquatic.index)
    for i in otu_level_df_filter_standardscale_aquatic.index:
        for j in otu_level_df_filter_standardscale_aquatic.index:
            bray_curtis_matrix_aquatic.loc[i, j] = braycurtis(otu_level_df_filter_standardscale_aquatic.loc[i],
                                                              otu_level_df_filter_standardscale_aquatic.loc[j])
    bray_curtis_matrix_aquatic.sort_index(inplace=True)
    bray_curtis_matrix_aquatic.sort_index(axis=1, inplace=True)
    result_matrix = np.full_like(bray_curtis_matrix_aquatic, -1)
    for i in range(bray_curtis_matrix_aquatic.shape[0]):
        for j in range(i + 1, bray_curtis_matrix_aquatic.shape[1]):
            result_matrix[i, j] = bray_curtis_matrix_aquatic.iloc[i, j]
    bray_curtis_matrix_aquatic_half = pd.DataFrame(result_matrix, index=bray_curtis_matrix_aquatic.index,
                                                   columns=bray_curtis_matrix_aquatic.columns)
    bray_curtis_matrix_aquatic_list = [i for i in bray_curtis_matrix_aquatic_half.values.flatten() if i != -1]
    all_level_set = set()
    for i in data_df_dict.keys():
        all_level_set |= set(data_df_dict[i][args.level_terrestrial])
    all_level_set.remove('unknown')
    all_otu_df_list = list()
    for i in data_df_dict.keys():
        level_list = [m for m, n in data_df_dict[i][args.level_terrestrial].value_counts().items() if m != 'unknown']
        count_list = [n for m, n in data_df_dict[i][args.level_terrestrial].value_counts().items() if m != 'unknown']
        diff_set = all_level_set - set(level_list)
        level_list.extend(list(diff_set))
        count_list.extend([0] * len(diff_set))
        all_otu_df_list.append(pd.DataFrame(count_list, index=level_list))
    otu_level_df = pd.concat([i.T for i in all_otu_df_list], ignore_index=True)
    otu_level_df.index = data_df_dict.keys()
    del_columns_list = list()
    for i in otu_level_df.columns:
        if all([j < 5 for j in otu_level_df[i]]):
            del_columns_list.append(i)
    otu_level_df_filter = otu_level_df.drop(columns=del_columns_list)
    otu_level_df_filter = otu_level_df_filter.loc[:, otu_level_df_filter.columns.isin(terrestrial_level_list)]
    otu_level_df_filter_standardscale_terrestrial = otu_level_df_filter.apply(
        lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
    otu_level_df_filter_standardscale_terrestrial.sort_index(inplace=True)
    bray_curtis_matrix_terrestrial = pd.DataFrame(index=otu_level_df_filter_standardscale_terrestrial.index,
                                                  columns=otu_level_df_filter_standardscale_terrestrial.index)
    for i in otu_level_df_filter_standardscale_terrestrial.index:
        for j in otu_level_df_filter_standardscale_terrestrial.index:
            bray_curtis_matrix_terrestrial.loc[i, j] = braycurtis(
                otu_level_df_filter_standardscale_terrestrial.loc[i],
                otu_level_df_filter_standardscale_terrestrial.loc[j])
    bray_curtis_matrix_terrestrial.sort_index(inplace=True)
    bray_curtis_matrix_terrestrial.sort_index(axis=1, inplace=True)
    result_matrix = np.full_like(bray_curtis_matrix_terrestrial, -1)
    for i in range(bray_curtis_matrix_terrestrial.shape[0]):
        for j in range(i + 1, bray_curtis_matrix_terrestrial.shape[1]):
            result_matrix[i, j] = bray_curtis_matrix_terrestrial.iloc[i, j]
    bray_curtis_matrix_terrestrial_half = pd.DataFrame(result_matrix, index=bray_curtis_matrix_terrestrial.index,
                                                       columns=bray_curtis_matrix_terrestrial.columns)
    bray_curtis_matrix_terrestrial_list = [i for i in bray_curtis_matrix_terrestrial_half.values.flatten() if i != -1]
    slope_a, intercept_a, r_value_a, p_value_a, std_err_a = stats.linregress(bray_curtis_matrix_aquatic_list,
                                                                             environmental_distances_list)
    slope_t, intercept_t, r_value_t, p_value_t, std_err_t = stats.linregress(bray_curtis_matrix_terrestrial_list,
                                                                             environmental_distances_list)
    sns.regplot(x=bray_curtis_matrix_aquatic_list, y=environmental_distances_list, scatter_kws={'alpha': 0.4},
                label='Aquatic')
    sns.regplot(x=bray_curtis_matrix_terrestrial_list, y=environmental_distances_list, scatter_kws={'alpha': 0.4},
                label='Terrestrial')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    text_w = f'Slope = {slope_a:.3f}, R² = {r_value_a ** 2:.3f} {"***" if p_value_a < 0.001 else "**" if p_value_a < 0.01 else "*" if p_value_a < 0.05 else ""}'
    text_g = f'Slope = {slope_t:.3f}, R² = {r_value_t ** 2:.3f} {"***" if p_value_t < 0.001 else "**" if p_value_t < 0.01 else "*" if p_value_t < 0.05 else ""}'
    plt.text(0.2, 0.75, text_w, ha='center', va='center', transform=plt.gca().transAxes, fontsize=14,
             color=sns.color_palette()[0], alpha=1, weight='bold')
    plt.text(0.2, 0.70, text_g, ha='center', va='center', transform=plt.gca().transAxes, fontsize=14,
             color=sns.color_palette()[1], alpha=1, weight='bold')
    plt.xlabel('Environmental Distance', fontsize=14)
    plt.ylabel('Assembling Dissimilarity', fontsize=14)
    plt.savefig(fname=f'{args.output}.{args.format}', bbox_inches='tight')


def netVC(args):
    groups_list = [i.split('/') for i in args.groups]
    df = pd.DataFrame(
        index=[f"{j} {chr(i + 65)}" for j in ['Aquatic', 'Terrestrial'] for i in range(len(groups_list))],
        columns=['Degree', 'Complexity', 'Total nodes', 'Total links', 'Total edges', 'Density', 'Modularity'])
    aquatic_level_list = args.aquatic
    aquatic_tmp_list = []
    for i in data_df_dict.keys():
        for j in aquatic_level_list:
            aquatic_tmp_list.extend(
                data_df_dict[i][data_df_dict[i][args.level_aquatic] == j][args.level_terrestrial].tolist())
    aquatic_tmp_list = list(set(aquatic_tmp_list))
    terrestrial_level_list = []
    for i in data_df_dict.keys():
        for j in aquatic_tmp_list:
            terrestrial_level_list.extend(
                data_df_dict[i][data_df_dict[i][args.level_terrestrial] != j][args.level_terrestrial].tolist())
    terrestrial_level_list = [i for i in set(terrestrial_level_list) if i != 'unknown']
    all_level_set = set()
    for i in data_df_dict.keys():
        all_level_set |= set(data_df_dict[i][args.level_aquatic])
    all_level_set.remove('unknown')
    all_otu_df_list = list()
    for i in data_df_dict.keys():
        level_list = [m for m, n in data_df_dict[i][args.level_aquatic].value_counts().items() if m != 'unknown']
        count_list = [n for m, n in data_df_dict[i][args.level_aquatic].value_counts().items() if m != 'unknown']
        diff_set = all_level_set - set(level_list)
        level_list.extend(list(diff_set))
        count_list.extend([0] * len(diff_set))
        all_otu_df_list.append(pd.DataFrame(count_list, index=level_list))
    otu_level_df = pd.concat([i.T for i in all_otu_df_list], ignore_index=True)
    otu_level_df.index = data_df_dict.keys()
    del_columns_list = list()
    for i in otu_level_df.columns:
        if all([j < 5 for j in otu_level_df[i]]):
            del_columns_list.append(i)
    otu_level_df_filter = otu_level_df.drop(columns=del_columns_list)
    otu_level_df_filter = otu_level_df_filter.loc[:, otu_level_df_filter.columns.isin(aquatic_level_list)]
    otu_level_df_filter.sort_index(inplace=True)
    otu_level_df_filter_standardscale_aquatic = otu_level_df_filter.apply(
        lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
    otu_level_df_filter_standardscale_aquatic.sort_index(inplace=True)

    data_list = []
    for group_list in groups_list:
        index_and_columns = otu_level_df_filter.loc[group_list, :].columns
        data = pd.DataFrame(index=index_and_columns, columns=index_and_columns).sort_index().sort_index(axis=1)
        if args.corr in ['pearson', 'spearman', 'kendall']:
            corr_dict = {'pearson': 'pearsonr', 'spearman': 'spearmanr', 'kendall': 'kendalltau'}
            corr = otu_level_df_filter_standardscale_aquatic.loc[group_list, :].corr(method=args.corr)
            p = otu_level_df_filter_standardscale_aquatic.loc[group_list, :].corr(
                method=lambda x, y: getattr(stats, corr_dict[args.corr])(x, y)[1])
        elif args.corr == 'sparcc':
            filter = otu_level_df_filter.loc[group_list, :]
            filter = filter.loc[:, (filter != 0).any(axis=0)]
            n_iteractions = 30
            x_iteractions = 30
            low_abundance = False
            threshold = 0.1
            normalization = 'dirichlet'
            log_transform = True
            num_simulate_data = 100
            type_pvalues = 'two_sided'
            SparCC_MN = SparCC_MicNet(n_iteractions=n_iteractions,
                                      x_iteractions=x_iteractions,
                                      low_abundance=low_abundance,
                                      threshold=threshold,
                                      normalization=normalization,
                                      log_transform=log_transform,
                                      num_simulate_data=num_simulate_data,
                                      type_pvalues=type_pvalues,
                                      )
            SparCC_MN.run_all(data_input=filter.T)
            corr = pd.read_csv(Path(SparCC_MN.save_corr_file).resolve(), index_col=0)
            corr.index = filter.columns
            corr.columns = filter.columns
            p = pd.read_csv(Path(SparCC_MN.outfile_pvals).resolve(), index_col=0)
            p.index = filter.columns
            p.columns = filter.columns
        for col1 in index_and_columns:
            for col2 in index_and_columns:
                if col1 != col2:
                    try:
                        _p = p.at[col1, col2]
                        _corr = corr.at[col1, col2]
                        if _p > args.p_value or abs(_corr) < args.corr_coef:
                            data.at[col1, col2] = 0
                        else:
                            data.at[col1, col2] = corr.at[col1, col2]
                    except Exception as e:
                        data.at[col1, col2] = np.NaN
        data_list.append(data)

        G = ig.Graph(directed=False)
        G.add_vertices(data.columns)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if not pd.isna(data.iloc[i, j]) and i != j and data.iloc[i, j] != 0:
                    if not G.are_connected(data.columns[i], data.index[j]):
                        edge_alpha = 1 if abs(data.iloc[i, j]) > args.corr_coef else 0
                        edge_color = f"rgba(254, 129, 125, {edge_alpha})" if data.iloc[
                                                                                 i, j] > 0 else f"rgba(129, 184, 223, {edge_alpha})"
                        G.add_edge(data.columns[i], data.index[j], weight=abs(data.iloc[i, j]), color=edge_color,
                                   curved=True)
        layout = G.layout_kamada_kawai(maxiter=2000, epsilon=1e-8)
        ig.plot(G, layout=layout, vertex_color='rgba(69, 189, 120, 1)', vertex_size=10, vertex_label_dist=1.5,
                edge_width=[1 * edge["weight"] for edge in G.es],
                edge_color=[edge["color"] for edge in G.es],
                #         vertex_label=G.vs["name"],
                curved=True,
                autocurve=True,
                bbox=(500, 500), margin=50,
                target=f'{args.output}_Aquatic_{chr(groups_list.index(group_list) + 65)}.{args.format}')
        plt.clf()
        degree = G.degree()
        complexity = [i / G.vcount() for i in degree]
        linked_nodes = set()
        for edge in G.es:
            linked_nodes.add(edge.tuple[0])
            linked_nodes.add(edge.tuple[1])
        total_linked_nodes = len(linked_nodes)
        total_nodes = G.vcount()
        total_edges = G.ecount()
        total_existed_nodes = len(
            otu_level_df_filter.loc[group_list, :].columns[otu_level_df_filter.loc[group_list, :].any()])
        density = 2 * total_edges / (total_existed_nodes * (total_existed_nodes - 1))
        if G.ecount() == 0:
            modularity = 0
        else:
            modularity = G.modularity(G.community_multilevel(), weights='weight')
        df.loc[f'Aquatic {chr(groups_list.index(group_list) + 65)}'] = [degree, complexity, total_nodes,
                                                                        total_linked_nodes, total_edges, density,
                                                                        modularity]

    color_list = ['#ff7f0e', '#2ca02c', '#db4142', '#9467bd', '#1f77b4', '#2c9f9f', '#7f7f7f', '#ffbb78', '#d62728',
                  '#bcbd22']
    ratio_list = [i / 100 for i in range(81)]
    count = 0
    for i in data_list:
        x_list = []
        y_list = []
        node_count_old = len(i)
        for ratio in ratio_list:
            nodes_to_remove = int(node_count_old * ratio)
            for _ in range(100):
                random_nodes_to_remove = random.sample(i.index.tolist(), nodes_to_remove)
                df_new = i.drop(random_nodes_to_remove, axis=0).drop(random_nodes_to_remove, axis=1)
                x_list.append(ratio)
                y_list.append(len(set(
                    [element for pair in df_new[df_new.notna() & (df_new != 0)].stack().index for element in
                     pair])) / node_count_old)
        ax = sns.regplot(x=x_list, y=y_list, scatter_kws={'alpha': 0.05, 'color': color_list[count]},
                         line_kws={'color': color_list[count], 'linewidth': 3},
                         label=f'Aquatic {chr(count + 65)}')
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_list, y_list)
        text = f'Slope = {slope:.3f}, R² = {r_value ** 2:.3f} {"***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""}'
        plt.text(0.2, 0.81 + 0.05 * count, text, ha='center', va='center', transform=plt.gca().transAxes, fontsize=14,
                 color=color_list[count], alpha=1, weight='bold')
        count += 1
    ax.set_ylim(bottom=-0.1, top=0.7)
    ax.set_yticks([i / 10 for i in range(8)])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    plt.xlabel('Proportion of removed nodes', fontsize=14)
    plt.ylabel('Natural connectivity', fontsize=14)
    plt.savefig(fname=f'{args.output}_Aquatic_Robustness.{args.format}', bbox_inches='tight')
    plt.clf()

    all_level_set = set()
    for i in data_df_dict.keys():
        all_level_set |= set(data_df_dict[i][args.level_terrestrial])
    all_level_set.remove('unknown')
    all_otu_df_list = list()
    for i in data_df_dict.keys():
        level_list = [m for m, n in data_df_dict[i][args.level_terrestrial].value_counts().items() if m != 'unknown']
        count_list = [n for m, n in data_df_dict[i][args.level_terrestrial].value_counts().items() if m != 'unknown']
        diff_set = all_level_set - set(level_list)
        level_list.extend(list(diff_set))
        count_list.extend([0] * len(diff_set))
        all_otu_df_list.append(pd.DataFrame(count_list, index=level_list))
    otu_level_df = pd.concat([i.T for i in all_otu_df_list], ignore_index=True)
    otu_level_df.index = data_df_dict.keys()
    del_columns_list = list()
    for i in otu_level_df.columns:
        if all([j < 5 for j in otu_level_df[i]]):
            del_columns_list.append(i)
    otu_level_df_filter = otu_level_df.drop(columns=del_columns_list)
    otu_level_df_filter = otu_level_df_filter.loc[:, otu_level_df_filter.columns.isin(terrestrial_level_list)]
    otu_level_df_filter.sort_index(inplace=True)
    otu_level_df_filter_standardscale_terrestrial = otu_level_df_filter.apply(
        lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
    otu_level_df_filter_standardscale_terrestrial.sort_index(inplace=True)

    data_list = []
    for group_list in groups_list:
        index_and_columns = otu_level_df_filter.loc[group_list, :].columns
        data = pd.DataFrame(index=index_and_columns, columns=index_and_columns).sort_index().sort_index(axis=1)
        if args.corr in ['pearson', 'spearman', 'kendall']:
            corr_dict = {'pearson': 'pearsonr', 'spearman': 'spearmanr', 'kendall': 'kendalltau'}
            corr = otu_level_df_filter_standardscale_terrestrial.loc[group_list, :].corr(method=args.corr)
            p = otu_level_df_filter_standardscale_terrestrial.loc[group_list, :].corr(
                method=lambda x, y: getattr(stats, corr_dict[args.corr])(x, y)[1])
        elif args.corr == 'sparcc':
            filter = otu_level_df_filter.loc[group_list, :]
            filter = filter.loc[:, (filter != 0).any(axis=0)]
            n_iteractions = 30
            x_iteractions = 30
            low_abundance = False
            threshold = 0.1
            normalization = 'dirichlet'
            log_transform = True
            num_simulate_data = 100
            type_pvalues = 'two_sided'
            SparCC_MN = SparCC_MicNet(n_iteractions=n_iteractions,
                                      x_iteractions=x_iteractions,
                                      low_abundance=low_abundance,
                                      threshold=threshold,
                                      normalization=normalization,
                                      log_transform=log_transform,
                                      num_simulate_data=num_simulate_data,
                                      type_pvalues=type_pvalues,
                                      )
            SparCC_MN.run_all(data_input=filter.T)
            corr = pd.read_csv(Path(SparCC_MN.save_corr_file).resolve(), index_col=0)
            corr.index = filter.columns
            corr.columns = filter.columns
            p = pd.read_csv(Path(SparCC_MN.outfile_pvals).resolve(), index_col=0)
            p.index = filter.columns
            p.columns = filter.columns
        for col1 in index_and_columns:
            for col2 in index_and_columns:
                if col1 != col2:
                    try:
                        _p = p.at[col1, col2]
                        _corr = corr.at[col1, col2]
                        if _p > args.p_value or abs(_corr) < args.corr_coef:
                            data.at[col1, col2] = 0
                        else:
                            data.at[col1, col2] = corr.at[col1, col2]
                    except Exception as e:
                        data.at[col1, col2] = np.NaN
        data_list.append(data)

        G = ig.Graph(directed=False)
        G.add_vertices(data.columns)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if not pd.isna(data.iloc[i, j]) and i != j and data.iloc[i, j] != 0:
                    if not G.are_connected(data.columns[i], data.index[j]):
                        edge_alpha = 1 if abs(data.iloc[i, j]) > args.corr_coef else 0
                        edge_color = f"rgba(254, 129, 125, {edge_alpha})" if data.iloc[
                                                                                 i, j] > 0 else f"rgba(129, 184, 223, {edge_alpha})"
                        G.add_edge(data.columns[i], data.index[j], weight=abs(data.iloc[i, j]), color=edge_color,
                                   curved=True)
        layout = G.layout_kamada_kawai(maxiter=2000, epsilon=1e-8)
        ig.plot(G, layout=layout, vertex_color='rgba(69, 189, 120, 1)', vertex_size=10, vertex_label_dist=1.5,
                edge_width=[1 * edge["weight"] for edge in G.es],
                edge_color=[edge["color"] for edge in G.es],
                #         vertex_label=G.vs["name"],
                curved=True,
                autocurve=True,
                bbox=(500, 500), margin=50,
                target=f'{args.output}_Terrestrial_{chr(groups_list.index(group_list) + 65)}.{args.format}')
        plt.clf()

        degree = G.degree()
        complexity = [i / G.vcount() for i in degree]
        linked_nodes = set()
        for edge in G.es:
            linked_nodes.add(edge.tuple[0])
            linked_nodes.add(edge.tuple[1])
        total_linked_nodes = len(linked_nodes)
        total_nodes = G.vcount()
        total_edges = G.ecount()
        total_existed_nodes = len(
            otu_level_df_filter.loc[group_list, :].columns[otu_level_df_filter.loc[group_list, :].any()])
        density = 2 * total_edges / (total_existed_nodes * (total_existed_nodes - 1))
        if G.ecount() == 0:
            modularity = 0
        else:
            modularity = G.modularity(G.community_multilevel(), weights='weight')
        df.loc[f'Terrestrial {chr(groups_list.index(group_list) + 65)}'] = [degree, complexity, total_nodes,
                                                                            total_linked_nodes, total_edges, density,
                                                                            modularity]

    color_list = ['#ff7f0e', '#2ca02c', '#db4142', '#9467bd', '#1f77b4', '#2c9f9f', '#7f7f7f', '#ffbb78', '#d62728',
                  '#bcbd22']
    ratio_list = [i / 100 for i in range(81)]
    count = 0
    for i in data_list:
        x_list = []
        y_list = []
        node_count_old = len(i)
        for ratio in ratio_list:
            nodes_to_remove = int(node_count_old * ratio)
            for _ in range(100):
                random_nodes_to_remove = random.sample(i.index.tolist(), nodes_to_remove)
                df_new = i.drop(random_nodes_to_remove, axis=0).drop(random_nodes_to_remove, axis=1)
                x_list.append(ratio)
                y_list.append(len(set(
                    [element for pair in df_new[df_new.notna() & (df_new != 0)].stack().index for element in
                     pair])) / node_count_old)
        ax = sns.regplot(x=x_list, y=y_list, scatter_kws={'alpha': 0.05, 'color': color_list[count]},
                         line_kws={'color': color_list[count], 'linewidth': 3},
                         label=f'Terrestrial {chr(count + 65)}')
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_list, y_list)
        text = f'Slope = {slope:.3f}, R² = {r_value ** 2:.3f} {"***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""}'
        plt.text(0.2, 0.81 + 0.05 * count, text, ha='center', va='center', transform=plt.gca().transAxes, fontsize=14,
                 color=color_list[count], alpha=1, weight='bold')
        count += 1
    ax.set_ylim(bottom=-0.1, top=0.7)
    ax.set_yticks([i / 10 for i in range(8)])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    plt.xlabel('Proportion of removed nodes', fontsize=14)
    plt.ylabel('Natural connectivity', fontsize=14)
    plt.savefig(fname=f'{args.output}_Terrestrial_Robustness.{args.format}', bbox_inches='tight')
    plt.clf()

    df['Level'] = ['Aquatic' if 'Aquatic' in index else 'Terrestrial' for index in df.index]
    df['Group'] = [i[-1] for i in df.index]

    Degree = df[['Degree', 'Level', 'Group']].explode('Degree').reset_index()
    ax = sns.boxplot(data=Degree, x='Group', order=[chr(i + 65) for i in range(len(groups_list))], y='Degree',
                     hue='Level', hue_order=['Aquatic', 'Terrestrial'])
    count = 0
    for i in [chr(i + 65) for i in range(len(groups_list))]:
        p_value = stats.mannwhitneyu(Degree[Degree['index'] == f'Aquatic {i}']['Degree'],
                                     Degree[Degree['index'] == f'Terrestrial {i}']['Degree'])[1]
        text = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "n.s."
        plt.text(0.12 + 0.25 * count, 0.82, text, ha='center', va='center', transform=plt.gca().transAxes, fontsize=22,
                 color='black', alpha=1, weight='bold')
        count += 1
    ax.set_ylim(bottom=-1, top=35)
    ax.set_yticks([i for i in range(0, 40, 5)])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('Degree', fontsize=14)
    plt.savefig(fname=f'{args.output}_Degree.{args.format}', bbox_inches='tight')
    plt.clf()

    Complexity = df[['Complexity', 'Level', 'Group']].explode('Complexity').reset_index()
    ax = sns.boxplot(data=Complexity, x='Group', order=[chr(i + 65) for i in range(len(groups_list))], y='Complexity',
                     hue='Level', hue_order=['Aquatic', 'Terrestrial'])
    count = 0
    for i in [chr(i + 65) for i in range(len(groups_list))]:
        p_value = stats.mannwhitneyu(Complexity[Complexity['index'] == f'Aquatic {i}']['Complexity'],
                                     Complexity[Complexity['index'] == f'Terrestrial {i}']['Complexity'])[1]
        text = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "n.s."
        plt.text(0.12 + 0.25 * count, 0.82, text, ha='center', va='center', transform=plt.gca().transAxes, fontsize=22,
                 color='black', alpha=1, weight='bold')
        count += 1
    # ax.set_ylim(bottom=-0.0125, top=0.25)
    # ax.set_yticks([i / 1000 for i in range(0, 275, 25)])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('Complexity', fontsize=14)
    plt.savefig(fname=f'{args.output}_Complexity.{args.format}', bbox_inches='tight')
    plt.clf()

    Total_nodes = df[['Total nodes', 'Level', 'Group']].explode('Total nodes').reset_index()
    ax = sns.barplot(data=Total_nodes, x='Group', order=[chr(i + 65) for i in range(len(groups_list))], y='Total nodes',
                     hue='Level', hue_order=['Aquatic', 'Terrestrial'])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('Total nodes', fontsize=14)
    plt.savefig(fname=f'{args.output}_Total_nodes.{args.format}', bbox_inches='tight')
    plt.clf()

    Total_links = df[['Total links', 'Level', 'Group']].explode('Total links').reset_index()
    ax = sns.barplot(data=Total_links, x='Group', order=[chr(i + 65) for i in range(len(groups_list))], y='Total links',
                     hue='Level', hue_order=['Aquatic', 'Terrestrial'])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('Total links', fontsize=14)
    plt.savefig(fname=f'{args.output}_Total_links.{args.format}', bbox_inches='tight')
    plt.clf()

    Total_edges = df[['Total edges', 'Level', 'Group']].explode('Total edges').reset_index()
    ax = sns.barplot(data=Total_edges, x='Group', order=[chr(i + 65) for i in range(len(groups_list))], y='Total edges',
                     hue='Level', hue_order=['Aquatic', 'Terrestrial'])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('Total edges', fontsize=14)
    plt.savefig(fname=f'{args.output}_Total_edges.{args.format}', bbox_inches='tight')
    plt.clf()

    Density = df[['Density', 'Level', 'Group']].explode('Density').reset_index()
    ax = sns.barplot(data=Density, x='Group', order=[chr(i + 65) for i in range(len(groups_list))], y='Density',
                     hue='Level', hue_order=['Aquatic', 'Terrestrial'])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('Density', fontsize=14)
    plt.savefig(fname=f'{args.output}_Density.{args.format}', bbox_inches='tight')
    plt.clf()

    Modularity = df[['Modularity', 'Level', 'Group']].explode('Modularity').reset_index()
    ax = sns.barplot(data=Modularity, x='Group', order=[chr(i + 65) for i in range(len(groups_list))], y='Modularity',
                     hue='Level', hue_order=['Aquatic', 'Terrestrial'])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('Modularity', fontsize=14)
    plt.savefig(fname=f'{args.output}_Modularity.{args.format}', bbox_inches='tight')
    plt.clf()


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--input', type=str, required=True,
                               help='Absolute or relative path file.(e.g. "./data.xlsx")')
    parent_parser.add_argument('--format', type=str, default='png', choices=['png', 'svg', 'pdf'],
                               help='Output format.(default: %(default)s)')
    parser = argparse.ArgumentParser(description='A comprehensive visual software for study taphonome.')
    parser.add_argument('-v', '--version', action='version', version='TaphonomeAnalyst 2.0')
    subparsers = parser.add_subparsers(help='commands')

    clusterenv_parser = subparsers.add_parser(name='clusterenv',
                                              parents=[parent_parser],
                                              help='Hierarchical clustering-sedimentary environment. (Module Ⅳ and Ⅴ)\t[clustermap]')
    clusterenv_parser.add_argument('--level', type=str, required=True,
                                   choices=['order', 'family', 'genera', 'species'],
                                   help='Taxonomic level.(For both statistical and aquatic OTUs.)')
    clusterenv_parser.add_argument('--aquatic', type=str2list, default=None,
                                   help='Aquatic OTUs.(default: all OTUs)\t[e.g. "OTU1,OTU2,OTU3"]')
    clusterenv_parser.add_argument('--geochem', type=str, required=False,
                                   help='Absolute or relative path geochemical file.(e.g. "./geochem.xlsx")')
    clusterenv_parser.add_argument('--output', type=str, default='./clusterenv',
                                   help='Absolute path or relative path and filename.(default: %(default)s)')
    clusterenv_parser.set_defaults(func=clusterenv)

    divvenn_parser = subparsers.add_parser(name='divvenn', parents=[parent_parser],
                                           help='Venn diagram-sampling locations or environments. (Module Ⅳ)\t[venn]')
    divvenn_parser.add_argument('--level', type=str, default='family',
                                choices=['order', 'family', 'genera', 'species'],
                                help='Taxonomic level.(default: %(default)s)')
    divvenn_parser.add_argument('--groups', type=str2list, required=True,
                                help='Custom Groups.(Recommend to group the plots by environments or locations)\t[e.g. "plotA1/plotB2,plotB1/plotA2,plotC1/plotC2"]')
    divvenn_parser.add_argument('--output', type=str, default='./divvenn',
                                help='Absolute path or relative path and filename.(default: %(default)s)')
    divvenn_parser.set_defaults(func=divvenn)

    TGotus_parser = subparsers.add_parser(name='TGotus', parents=[parent_parser],
                                          help='Taphonomic grades-taxa. (Module Ⅲ)\t[barh]')
    TGotus_parser.add_argument('--level', type=str, default='order',
                               choices=['order', 'family', 'genera', 'species'],
                               help='Taxonomic level.(default: %(default)s)')
    TGotus_parser.add_argument('--output', type=str, default='./TGotus',
                               help='Absolute path or relative path and filename.(default: %(default)s)')
    TGotus_parser.set_defaults(func=TGotus)

    TGplots_parser = subparsers.add_parser(name='TGplots', parents=[parent_parser],
                                           help='Taphonomic grades-sampling plots (in customized order). (Module Ⅲ)\t[barh]')
    TGplots_parser.add_argument('--groups', type=str2list, default=None,
                                help='Environment groups.(Recommend to group the plots by different aquatic and terrestrial environments)\t[e.g. "plotA1/plotB2,plotB1/plotA2,plotC1/plotC2"]')
    TGplots_parser.add_argument('--output', type=str, default='./TGplots',
                                help='Absolute path or relative path and filename.(default: %(default)s)')
    TGplots_parser.set_defaults(func=TGplots)

    abundplots_parser = subparsers.add_parser(name='abundplots', parents=[parent_parser],
                                              help='Abundance-sampling plots. (Module Ⅱ)\t[barh]')
    abundplots_parser.add_argument('--level', type=str, default='order',
                                   choices=['order', 'family', 'genera', 'species'],
                                   help='Taxonomic level.(default: %(default)s)')
    abundplots_parser.add_argument('--output', type=str, default='./abundplots',
                                   help='Absolute path or relative path and filename.(default: %(default)s)')
    abundplots_parser.set_defaults(func=abundplots)

    cooccurnet_parser = subparsers.add_parser(name='cooccurnet', parents=[parent_parser],
                                              help='Co-occurrence networks. (Module Ⅸ)\t[network]')
    cooccurnet_parser.add_argument('--level', type=str, default='family',
                                   choices=['order', 'family', 'genera', 'species'],
                                   help='Taxonomic level.(default: %(default)s)')
    cooccurnet_parser.add_argument('--corr', type=str, default='pearson',
                                   choices=['pearson', 'spearman', 'kendall', 'sparcc'],
                                   help='Correlation algorithm.(default: %(default)s)')
    cooccurnet_parser.add_argument('--corr_coef', type=float, default=0.7,
                                   help='Minimum threshold of correlation coefficient.(default: %(default)s)')
    cooccurnet_parser.add_argument('--p_value', type=float, default=0.1,
                                   help='Maximum threshold of p-value.(default: %(default)s)')
    cooccurnet_parser.add_argument('--output', type=str, default='./cooccurnet',
                                   help='Absolute path or relative path and filename.(default: %(default)s)')
    cooccurnet_parser.set_defaults(func=cooccurnet)

    samplecurve_parser = subparsers.add_parser(name='samplecurve', parents=[parent_parser],
                                               help='Sampling coverage curve. (Module Ⅰ)\t[regplot]')
    samplecurve_parser.add_argument('--level', type=str, default='family',
                                    choices=['order', 'family', 'genera', 'species'],
                                    help='Taxonomic level.(default: %(default)s)')
    samplecurve_parser.add_argument('--groups', type=str2dictlist, required=True,
                                    help='Grouping plots (Sheet names) with customized names.\t[e.g. "plotA:plotA1/plotA2,plotB:plotB1/plotB2"]')
    samplecurve_parser.add_argument('--output', type=str, default='./samplecurve',
                                    help='Absolute path or relative path and filename.(default: %(default)s)')
    samplecurve_parser.set_defaults(func=samplecurve)

    chao_parser = subparsers.add_parser(name='chao', parents=[parent_parser],
                                        help='Chao1 potential diversity curve. (Module Ⅰ)\t[regplot]')
    chao_parser.add_argument('--level', type=str, default='family',
                             choices=['order', 'family', 'genera', 'species'],
                             help='Taxonomic level.(default: %(default)s)')
    chao_parser.add_argument('--groups', type=str2dictlist, required=True,
                             help='Grouping plots (Sheet names) with customized names.\t[e.g. "plotA:plotA1/plotA2,plotB:plotB1/plotB2"]')
    chao_parser.add_argument('--output', type=str, default='./chao',
                             help='Absolute path or relative path and filename.(default: %(default)s)')
    chao_parser.set_defaults(func=chao)

    ace_parser = subparsers.add_parser(name='ace', parents=[parent_parser],
                                       help='ACE potential diversity curve. (Module Ⅰ)\t[regplot]')
    ace_parser.add_argument('--level', type=str, default='family', choices=['order', 'family', 'genera', 'species'],
                            help='Taxonomic level.(default: %(default)s)')
    ace_parser.add_argument('--groups', type=str2dictlist, required=True,
                            help='Grouping plots (Sheet names) with customized names.\t[e.g. "plotA:plotA1/plotA2,plotB:plotB1/plotB2"]')
    ace_parser.add_argument('--output', type=str, default='./ace',
                            help='Absolute path or relative path and filename.(default: %(default)s)')
    ace_parser.add_argument('--rare', type=int, default=10, help='ACE rare threshold.(default: %(default)s)')
    ace_parser.set_defaults(func=ace)

    corrotus_parser = subparsers.add_parser(name='corrotus', parents=[parent_parser],
                                            help='Heatmap-OTUs correlation analysis. (Module Ⅷ)\t[heatmap]')
    corrotus_parser.add_argument('--level', type=str, default='family',
                                 choices=['order', 'family', 'genera', 'species'],
                                 help='Taxonomic level.(default: %(default)s)')
    corrotus_parser.add_argument('--corr', type=str, default='pearson',
                                 choices=['pearson', 'spearman', 'kendall', 'sparcc'],
                                 help='Correlation algorithm.(default: %(default)s)')
    corrotus_parser.add_argument('--p_value', type=float, default=0.1,
                                 help='Maximum threshold of p-value.(default: %(default)s)')
    corrotus_parser.add_argument('--output', type=str, default='./corrotus',
                                 help='Absolute path or relative path and filename.(default: %(default)s)')
    corrotus_parser.set_defaults(func=corrotus)

    mantel_parser = subparsers.add_parser(name='mantel', parents=[parent_parser],
                                          help='Mantel Test between species abundance and ecological environmental variables. (Module Ⅶ)\t[multiplot]')
    mantel_parser.add_argument('--rhome', type=str, required=True,
                               help='Absolute path of R_HOME.(e.g. "C:\Program Files\R\R-4.3.2")')
    mantel_parser.add_argument('--geochem', type=str, required=True,
                               help='Absolute or relative path geochemical file.(e.g. "./geochem.xlsx")')
    mantel_parser.add_argument('--aquatic', type=str2list, required=True,
                               help='Aquatic OTUs.\t[e.g. "OTU1,OTU2,OTU3"]')
    mantel_parser.add_argument('--level_aquatic', type=str, required=True,
                               choices=['order', 'family', 'genera', 'species'],
                               help="Taxonomic level for aquatic OTUs.")
    mantel_parser.add_argument('--level_terrestrial', type=str, required=True,
                               choices=['order', 'family', 'genera', 'species'],
                               help="Taxonomic level for terrestrial OTUs.")
    mantel_parser.add_argument('--corr', type=str, default='pearson',
                               choices=['pearson', 'spearman', 'kendall'],
                               help='Correlation algorithm for geochem.(default: %(default)s)')
    mantel_parser.add_argument('--output', type=str, default='./mantel',
                               help='Absolute path or relative path and filename.(default: %(default)s)')
    mantel_parser.set_defaults(func=mantel)

    dissenvtest_parser = subparsers.add_parser(name='dissenvtest', parents=[parent_parser],
                                               help='Assembling dissimilarity- environmental distance test. (Module Ⅵ)\t[regplot]')
    dissenvtest_parser.add_argument('--geochem', type=str, required=True,
                                    help='Absolute or relative path geochemical file.(e.g. "./geochem.xlsx")')
    dissenvtest_parser.add_argument('--aquatic', type=str2list, required=True,
                                    help='Aquatic OTUs.\t[e.g. "OTU1,OTU2,OTU3"]')
    dissenvtest_parser.add_argument('--level_aquatic', type=str, required=True,
                                    choices=['order', 'family', 'genera', 'species'],
                                    help="Taxonomic level for aquatic OTUs.")
    dissenvtest_parser.add_argument('--level_terrestrial', type=str, required=True,
                                    choices=['order', 'family', 'genera', 'species'],
                                    help="Taxonomic level for terrestrial OTUs.")
    dissenvtest_parser.add_argument('--output', type=str, default='./dissenvtest',
                                    help='Absolute path or relative path and filename.(default: %(default)s)')
    dissenvtest_parser.set_defaults(func=dissenvtest)

    netVC_parser = subparsers.add_parser(name='netVC', parents=[parent_parser],
                                         help='Generate a unified layout network for comparison. (Module Ⅹ)\t[network/boxplot/barplot]')
    netVC_parser.add_argument('--aquatic', type=str2list, required=True,
                              help='Species-level aquatic OTUs.\t[e.g. "OTU1,OTU2,OTU3"]')
    netVC_parser.add_argument('--level_aquatic', type=str, required=True,
                              choices=['order', 'family', 'genera', 'species'],
                              help="Taxonomic level for aquatic OTUs.")
    netVC_parser.add_argument('--level_terrestrial', type=str, required=True,
                              choices=['order', 'family', 'genera', 'species'],
                              help="Taxonomic level for terrestrial OTUs.")
    netVC_parser.add_argument('--groups', type=str2list, required=True,
                              help='Environment groups.(Grouping the plots by different aquatic and terrestrial environments)\t[e.g. "plotA1/plotB2,plotB1/plotA2,plotC1/plotC2"]')
    netVC_parser.add_argument('--corr', type=str, default='pearson',
                              choices=['pearson', 'spearman', 'kendall', 'sparcc'],
                              help='Correlation algorithm.(default: %(default)s)')
    netVC_parser.add_argument('--corr_coef', type=float, default=0.7,
                              help='Minimum threshold of correlation coefficient.(default: %(default)s)')
    netVC_parser.add_argument('--p_value', type=float, default=0.1,
                              help='Maximum threshold of p-value.(default: %(default)s)')
    netVC_parser.add_argument('--output', type=str, default='./netVC',
                              help='Absolute path or relative path and filename.(default: %(default)s)')
    netVC_parser.set_defaults(func=netVC)

    args = parser.parse_args()

    data_file = args.input
    data_xls = pd.ExcelFile(data_file)
    all_sheet_name_list = [i for i in data_xls.sheet_names if not i.startswith('Sheet')]
    data_df_dict = dict()
    for i in all_sheet_name_list:
        data_df_dict[i] = pd.read_excel(data_xls, sheet_name=i, header=0,
                                        usecols=['order', 'family', 'genera', 'species', 'taphonomic grade'],
                                        nrows=3000)
    # if data_df_dict[i].shape[0] < 3000:
    #     data_df_dict[i] = pd.concat([data_df_dict[i]] * ceil(3000 / data_df_dict[i].shape[0]), ignore_index=True)
    #     data_df_dict[i] = data_df_dict[i][:3000]

    color_list = ['#e64b35', '#e18727', '#ffdc91', '#0072b5', '#4dbbd5', '#00a087', '#925e9f', '#ee4c97', '#f39b7f',
                  '#7e6148',
                  '#19b4ca', '#298af0', '#00236e', '#ff8d4a', '#b2442a', '#ff5f78', '#6da160', '#11b368', '#0c6480',
                  '#819eb7',
                  '#fffc00', '#ff6c00', '#ff0000', '#9cff00', '#00ffd2', '#166700', '#bcadff', '#af1a5d', '#9c00ff',
                  '#9c9c9c',
                  '#1e00ff', '#a9ffea', '#7d6c4b', '#ffd409', '#502e03', '#e998ff', '#435200', '#50e5a2', '#bbff90',
                  '#000000']
    grade_list = ['A', 'B', 'C', 'D', 'E']

    args.func(args)
