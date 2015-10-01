from base64 import b64encode
from collections import OrderedDict
from io import BytesIO

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

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


def fragment_length_plot(data_file, peaks=None):
    data = read_picard_histogram(data_file)

    fig = plt.figure()
    plt.bar(data[:, 0], data[:, 1])
    plt.xlim((0, 1000))

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
