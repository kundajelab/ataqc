#! user/bin/env python

import matplotlib
matplotlib.use('Agg')

import gzip
import pandas as pd
import numpy as np

from base64 import b64encode
from io import BytesIO
from collections import OrderedDict
from matplotlib import pyplot as plt

class Peaks():

    def __init__(self, peak_file, peak_file_name):
        self.peak_file = peak_file
        self.peak_file_name = peak_file_name
        self.metrics = {}



    def get_region_size_metrics(self):
        peak_size_summ = OrderedDict([# Create default quartile metrics (0 for everything)
            ('Min size', 0),
            ('25 percentile', 0),
            ('50 percentile (median)', 0),
            ('75 percentile', 0),
            ('Max size', 0),
            ('Mean', 0),
        ])

        # Load peak file. If it fails, return nothing as above
        try:
            peak_df = pd.read_table(self.peak_file, compression='gzip', header=None)# Attempt to read peak file into a panda dataframe
        except:
            return peak_size_summ, ''# If peak file cannot be read, return empty quartile metrics
        
        # Subtract third column from second to get summary
        region_sizes = peak_df.ix[:,2] - peak_df.ix[:,1]
        
        # Summarize and store in ordered dict
        peak_summary_stats = region_sizes.describe()
        
        peak_size_summ = OrderedDict([
            ('Min size', peak_summary_stats['min']),
            ('25 percentile', peak_summary_stats['25%']),
            ('50 percentile (median)', peak_summary_stats['50%']),
            ('75 percentile', peak_summary_stats['75%']),
            ('Max size', peak_summary_stats['max']),
            ('Mean', peak_summary_stats['mean']),
        ])

        # Plot density diagram using matplotlib
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        y, binEdges = np.histogram(region_sizes, bins=100)
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])

        # density = gaussian_kde(y) # from scipy.stats import gaussian_kde
        # density.covariance_factor = lambda : .25
        # density._compute_covariance()

        plt.plot(bincenters, y, '-')
        filename = self.peak_file.split('/')[-1]
        ax.set_title('Peak width distribution for {0}'.format(filename))
        #ax.set_yscale('log')
        
        plot_img = BytesIO()
        fig.savefig(plot_img, format='png')
        
        return peak_size_summ, b64encode(plot_img.getvalue())



    def count(self):
        return sum(1 for line in gzip.open(self.peak_file))



    def run_metrics(self, mode='all_metrics'):
        self.metrics['name'] = self.peak_file_name
        self.metrics['sizes'] = self.get_region_size_metrics()
        self.metrics['peak_count'] = self.count()
           
        return self.metrics

    def get_name(self):
        return self.peak_file_name

    def get_sizes(self):
        return self.metrics['sizes']

    def get_peak_count(self):
        return self.metrics['peak_count']
