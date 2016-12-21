def generate_metric_list():
    # Generate the set of metric objects to call based on parameters passed in
    # If a required parameter is None, skip (or return a warning)
    pass


def generate_report(metric_list):
    pass

def generate_dict(metric_list):
    pass

metric_list = generate_metric_list()

for metric in metric_list:
    metric.compute()
    metric.html_display()
    metric.to_dict()

