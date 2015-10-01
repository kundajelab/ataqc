from base64 import b64encode
from collections import namedtuple
from collections import OrderedDict
from io import BytesIO

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks_cwt

from jinja2 import Template


def preseq_plot(data_file):
    data = np.loadtxt(data_file, skiprows=1)
    data /= 1e6  # scale to millions of reads

    fig = plt.figure()

    # Plot the average expected yield
    plt.plot(data[:, 0], data[:, 1], 'r-')

    # Plot confidence intervals
    ci_lower, = plt.plot(data[:, 0], data[:, 2], 'b--')
    ci_upper, = plt.plot(data[:, 0], data[:, 3], 'b--')
    plt.legend([ci_lower], ['95% confidence interval'], loc=4)

    plt.title('Preseq estimated yield')
    plt.xlabel('Sequenced fragments [ millions ]')
    plt.ylabel('Expected distinct fragments [ millions ]')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return plot_img.getvalue()


def read_picard_histogram(data_file):
    with open(data_file) as fp:
        for line in fp:
            if line.startswith('## HISTOGRAM'):
                break
        data = np.loadtxt(fp, skiprows=1)

    return data


QCResult = namedtuple('QCResult', ['metric', 'qc_pass', 'message'])
INF = float("inf")


class QCCheck(object):
    def __init__(self, metric):
        self.metric = metric

    def check(self, value):
        return True

    def message(self, value, qc_pass):
        return ('{} - OK'.format(value) if qc_pass
                else '{} - Failed'.format(value))

    def __call__(self, value):
        qc_pass = self.check(value)
        return QCResult(self.metric, qc_pass, self.message(value, qc_pass))


class QCIntervalCheck(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCIntervalCheck, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, value):
        return self.lower <= value <= self.upper

    def message(self, value, qc_pass):
        return ('{} - OK'.format(value) if qc_pass else
                '{} out of range [{}, {}]'.format(value, self.lower,
                                                  self.upper))


class QCLessThanEqualCheck(QCIntervalCheck):
    def __init__(self, metric, upper):
        super(QCLessThanEqualCheck, self).__init__(metric, -INF, upper)


class QCGreaterThanEqualCheck(QCIntervalCheck):
    def __init__(self, metric, lower):
        super(QCGreaterThanEqualCheck, self).__init__(metric, lower, INF)


class QCHasElementInRange(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCHasElementInRange, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, elems):
        return (len([elem for elem in elems
                    if self.lower <= elem <= self.upper]) > 0)

    def message(self, elems, qc_pass):
        return ('OK' if qc_pass else
                'Cannot find element in range [{}, {}]'.format(
                    self.lower, self.upper))


def fragment_length_qc(data):
    results = []

    NFR_UPPER_LIMIT = 150
    MONO_NUC_LOWER_LIMIT = 150
    MONO_NUC_UPPER_LIMIT = 300

    # % of NFR vs res
    percent_nfr = data[:NFR_UPPER_LIMIT].sum() / data.sum()
    results.append(
        QCGreaterThanEqualCheck('Fraction of reads in NFR', 0.4)(percent_nfr))

    # % of NFR vs mononucleosome
    percent_nfr_vs_mono_nuc = (
        data[:NFR_UPPER_LIMIT].sum() /
        data[MONO_NUC_LOWER_LIMIT:MONO_NUC_UPPER_LIMIT + 1].sum())
    results.append(
        QCGreaterThanEqualCheck('NFR / mono-nuc reads', 2.5)(
            percent_nfr_vs_mono_nuc))

    # peak locations
    peaks = find_peaks_cwt(data[:, 1], np.array([25]))
    nuc_range_metrics = [('Presence of NFR peak', 20, 90),
                         ('Presence of Mono-Nuc peak', 120, 250),
                         ('Presence of Di-Nuc peak', 300, 500)]
    for range_metric in nuc_range_metrics:
        results.append(QCHasElementInRange(*range_metric)(peaks))

    return results


def fragment_length_plot(data_file, peaks=None):
    data = read_picard_histogram(data_file)

    fig = plt.figure()
    plt.bar(data[:, 0], data[:, 1])
    plt.xlim((0, 1000))

    if peaks:
        peak_vals = [data[peak_x, 1] for peak_x in peaks]
        plt.plot(peaks, peak_vals, 'ro')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return plot_img.getvalue()


html_template = Template("""
{% macro inline_img(base64_img, img_type='png') -%}
    <img src="data:image/{{ img_type }};base64,{{ base64_img }}">
{%- endmacro %}

<html>

<head>
  <title>{{ sample['name'] }} - ATAqC report</title>
</head>

<body>
  <h2>Basic Information</h2>
  <table>
    <tbody>
      {% for field, value in sample['basic_info'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>

  <h2>Alignment statistics</h2>


  <h2>Fragment length distribution</h2>
  {{ inline_img(sample['fraglen_dist']) }}

  <h2>Enrichment plots</h2>
  <h3>TSS enrichment plot</h3>
  {{ inline_img(sample['enrichment_plots']['tss']) }}

  <h2>Library complexity</h2>

  <h2>Yield prediction</h2>
  {{ inline_img(sample['yield_prediction']) }}
</body>

</html>
""")

TEST_BASIC_INFO = OrderedDict([
    ('Filename', 'some_atac_sample'),
    ('Genome', 'hg19'),
])

TEST_ENRICHMENT_PLOTS = {
    'tss': b64encode(open('test.png', 'rb').read())
}

TEST_SAMPLE = {
    'name': 'some_atac_sample',
    'basic_info': TEST_BASIC_INFO,
    'enrichment_plots': TEST_ENRICHMENT_PLOTS,
    'yield_prediction': b64encode(preseq_plot('test.preseq.dat')),
    'fraglen_dist': b64encode(fragment_length_plot('test_hist_data.log')),
}

print html_template.render(sample=TEST_SAMPLE)
