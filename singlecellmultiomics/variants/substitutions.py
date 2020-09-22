#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from collections import OrderedDict
from itertools import product

mpl.rcParams['figure.dpi'] = 300


def conversion_dict():
    conversions_single_nuc = ("CA", "CG", "CT", "TA", "TC", "TG")
    pattern_counts = OrderedDict()
    for ref, to in conversions_single_nuc:
        for context in product('ACGT',repeat=2 ):
            pattern_counts[(f'{context[0]}{ref}{context[1]}', to)] = 0
    return pattern_counts


def conversion_dict_stranded():
    conversions_single_nuc = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')
    pattern_counts = OrderedDict()
    for ref, to in conversions_single_nuc:
        for context in product('ACGT',repeat=2 ):
            pattern_counts[(f'{context[0]}{ref}{context[1]}', to)] = 0
    return pattern_counts


def substitution_plot(pattern_counts: dict,
                      figsize: tuple = (10, 4),
                      conversion_colors: tuple = None,

                      ylabel: str = '# conversions per molecule',
                      add_main_group_labels: bool = True,
                      fig=None,
                      ax=None,
                      stranded=False,
                      **plot_args
                      ):
    """
    Create 3bp substitution plot

    Args:
        pattern_counts(OrderedDict) : Dictionary containing the substitutions to plot.
            Use variants.vcf_to_variant_contexts to create it.
            Format:
            ```OrderedDict([(('ACA', 'A'), 0),
                 (('ACC', 'A'), 1),
                 (('ACG', 'A'), 0),
                 ...
                 (('TTG', 'G'), 0),
                 (('TTT', 'G'), 0)])```

        figsize(tuple) : size of the figure to create

        conversion_colors(tuple) : colors to use for the conversion groups

        ylabel(str) : y axis label

        add_main_group_labels(bool) : Add conversion group labels to top of plot

         **plot_args : Additional argument to pass to .plot()

    Returns
        fig : handle to the figure
        ax : handle to the axis

    Example:
        >>> from singlecellmultiomics.variants import vcf_to_variant_contexts, substitution_plot
        >>> import matplotlib.pyplot as plt
        >>> pobs = vcf_to_variant_contexts('variants.vcf.gz', 'reference.fasta')
        >>> for sample, conversions in pobs.items():
        >>>     fig, ax = substitution_plot(conversions)
        >>>     ax.set_title(sample)
        >>>     plt.show()

    """
    if stranded:
        conversions_single_nuc = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')
        conversion_colors = conversion_colors if conversion_colors is not None else ('b', 'k', 'r', 'grey', 'g', 'pink', 'b', 'k', 'r', 'k', 'w', 'g'),
    else:
        conversions_single_nuc = ("CA", "CG", "CT", "TA", "TC", "TG")
        conversion_colors = conversion_colors if conversion_colors is not None else ('b', 'k', 'r', 'grey', 'g', 'pink')

    # Colors for the conversion groups:
    color_d = dict(zip(conversions_single_nuc, conversion_colors))
    colors = [color_d.get(f'{context[1]}{to}') for context, to in pattern_counts.keys()]

    if fig is None or ax is None:
        assert fig is None and ax is None
        fig, ax = plt.subplots(figsize=figsize)

    substitution_dataframe = pd.DataFrame(pattern_counts.values(), index=list(pattern_counts.keys())).T
    substitution_dataframe.plot(kind='bar', color=colors, legend=False, width=1.0, ax=ax, edgecolor='k', **plot_args)
    offset = (1 / len(pattern_counts)) * 0.5  # Amount of distance for a half bar

    # Add 3bp context ticks:
    ax.set_xticks(np.linspace(-0.5 + offset, 0.5 - offset, len(pattern_counts)))
    ax.set_xticklabels( [context for context, to in pattern_counts.keys()], rotation=90, size=6)
    ax.set_ylabel(ylabel)
    ax.set_xlim((-0.5, 0.5))

    sns.despine()
    if add_main_group_labels:
        for i, (u, v) in enumerate(conversions_single_nuc):
            ax.text(  # position text relative to Axes
                (i + 0.5) / len(conversions_single_nuc), 1.0, f'{u}>{v}', fontsize=8,
                ha='center', va='top',
                transform=ax.transAxes, bbox=dict(facecolor='white', alpha=1, lw=0)
            )

    ax.set_axisbelow(True)
    ax.grid(axis='y')

    return fig, ax
