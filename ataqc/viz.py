# functions to take QC objects and write them to html, txt, multiQC

from jinja2 import Template


def html_header(sample_name):
    template = Template(
        """
        <html>

        <head>
        <title>ATAqC - {{ sample_name }}</title>
        <style>
        .qc_table{
            font-family:"Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
            font-size:12px;
            width:480px;
            text-align:left;
            border-collapse:collapse;
            margin:20px;
        }

        .qc_table th{
            font-size:14px;
            font-weight:normal;
            background:#8c1515;
            border-top:4px solid #700000;
            border-bottom:1px solid #fff;
            color:white;
            padding:8px;
        }

        .qc_table td{
            background:#f2f1eb;
            border-bottom:1px solid #fff;
            color:black;
            border-top:1px solid transparent;
            padding:8px;
        }

        .qc_table .fail{
            color:#ff0000;
            font-weight:bold;
        }
        </style>
        </head>

        <body>
        """)

    return template.render(sample_name=sample_name)


def html_footer():
    template = Template(
        """
        </body>
        </html>
        """)
    return template.render()

def html_section_header(header):
    template = Template(
        """
        {% if header is not none %}
            <h2>{{ header }}</h2>
        {% endif %}
        """)
    return template.render(header=header)


def html_text_description(description):
    template = Template(
        """
        {% if description is not none %}
        <pre>{{ description }}</pre>
        {% endif %}
        """)
    return template.render(description=description)


def html_metric_and_fraction_table(header, description, metric_dict):
    template = Template(
        """
        {% if header is not none %}
            <h3>{{ header }}</h3>
        {% endif %}

        <table class='qc_table'>
            <thead>
            <tr>
                <th scope='col'>Metric</th>
                <th scope='col'>Value</th>
                <th scope='col'>Fraction</th>
            </tr>
            </thead>

            <tbody>
            {% for field, value in metric_dict.iteritems() %}
            <tr>
                <td>{{ field }}</td>
                <td>{{ '{0:,}'.format(value[0]) }}</td>
                <td>{{ '{0:.3f}'.format(value[1]) }}</td>
            </tr>
            {% endfor %}
            </tbody>
        </table>

        {% if description is not none %}
        <pre>{{ description }}</pre>
        {% endif %}
        """)
    return template.render(header=header,
                           description=description,
                           metric_dict=metric_dict)


def html_metric_table(header, description, metric_dict):
    template = Template(
        """
        {% if header is not none %}
            <h3>{{ header }}</h3>
        {% endif %}

        <table class='qc_table'>
            <thead>
            <tr>
                <th scope='col'>Metric</th>
                <th scope='col'>Value</th>
            </tr>
            </thead>

            <tbody>
            {% for field, value in metric_dict.iteritems() %}
            <tr>
                <td>{{ field }}</td>
                <td>{{ value }}</td>
            </tr>
            {% endfor %}
            </tbody>
        </table>

        {% if description is not none %}
        <pre>{{ description }}</pre>
        {% endif %}
        """)
    return template.render(header=header,
                           description=description,
                           metric_dict=metric_dict)


def qc_to_html(qc_object):
    """Take a qc dict and write into html format"""

    if qc_object['type'] == 'plot':
        pass
    
    elif qc_object['type'] == 'qc_table':
        pass
    
    elif qc_object['type'] == 'table':
        rendered = html_metric_table(qc_object['header'],
                                     qc_object['description'],
                                     qc_object['qc'])
        
    elif qc_object['type'] == 'metric_fraction_table':
        rendered = html_metric_and_fraction_table(qc_object['header'],
                                                  qc_object['description'],
                                                  qc_object['qc'])
        
    else:
        rendered = ''

    return rendered



def write_html(qc_groups, args):
    """Take the QC and write it all out into HTML"""

    rendered = ""
    
    # write title and other useful stuff
    rendered += html_header('test')

    for qc_group in qc_groups:
        
        # TODO write header
        rendered += html_section_header(qc_group.get_name())
        
        qc_ordered_dict = qc_group.get_qc()
        
        for qc_object in qc_ordered_dict.items():
            print qc_object
            rendered += qc_to_html(qc_object[1])
            
        rendered += html_text_description(qc_group.get_description())

    rendered += html_footer()

    with open('test_qc.html', 'w') as out:
        out.write(rendered)

    return None
