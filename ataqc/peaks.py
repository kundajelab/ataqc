

import matplotlib
matplotlib.use('Agg')

import gzip

from collections import OrderedDict
from matplotlib import pyplot as plt



class Peaks():


    def __init__(self, peak_file, peak_file_name):
        self.peak_file = peak_file
        self.peak_file_name = peak_file_name
        

    def get_region_size_metrics(self):
        '''
        From the peak file, return a plot of the region size distribution and
        the quartile metrics (summary from R)
        '''

        peak_size_summ = OrderedDict([
            ('Min size', 0),
            ('25 percentile', 0),
            ('50 percentile (median)', 0),
            ('75 percentile', 0),
            ('Max size', 0),
            ('Mean', 0),
        ])
            
        # Load peak file. If it fails, return nothing as above
        try:
            peak_df = pd.read_table(peak_file, compression='gzip', header=None)
        except:
            return peak_size_summ, ''
        
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
        filename = peak_file.split('/')[-1]
        ax.set_title('Peak width distribution for {0}'.format(filename))
        #ax.set_yscale('log')
        
        plot_img = BytesIO()
        fig.savefig(plot_img, format='png')
        
        return peak_size_summ, b64encode(plot_img.getvalue())

    
    def count(self):
        '''Return peak count'''
        
        return sum(1 for line in gzip.open(self.peak_file))
    

    def run_metrics(self):
        """Run QC metrics"""

        metrics = {}
        metrics['name'] = self.peak_file_name
        metrics['sizes'] = self.get_region_size_metrics()
        metrics['peak_count'] = self.count()
        
        
        return metrics
