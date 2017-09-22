# description: helper functions

import os
import json

from collections import OrderedDict

def sample_sheet_to_dict(sample_sheet):
    """Take in a sample sheet and convert to a list of sample dicts
    """
    samples = []
    line_num = 0
    with open(sample_sheet, 'r') as fp:
        for line in fp:
            if line_num == 0:
                header = line.strip().split("\t")
                line_num += 1
                continue

            sample_info = line.strip().split("\t")

            sample_info_clean = []
            for elem in sample_info:
                if elem == "":
                    sample_info_clean.append(None)
                elif "=" in elem:
                    # make into dict
                    elem_dict = OrderedDict()
                    elem_fields = elem.split(";")
                    for elem_field in elem_fields:
                        key, val = elem_field.split("=")
                        elem_dict[key] = val
                    sample_info_clean.append(elem_dict)
                else:
                    sample_info_clean.append(elem)    
            sample_dict = dict(zip(header, sample_info_clean))
            samples.append(sample_dict)

    return samples


def unicode_to_byte(data, ignore_dicts=False):
    # if data is a unicode string, return its string representation
    if isinstance(data, unicode):
        return data.encode('utf-8')
    
    # if data is a list of unicode values, return list of byte-encoded values
    if isinstance(data, list):
        return [unicode_to_byte(element, ignore_dicts=True) for element in data]

    # if data is a dictionary, return dictionary of byte-encoded keys and values,
    # but only if key-value pairs have not already been converted
    if isinstance(data, dict) and not ignore_dicts:
        return {
            unicode_to_byte(key, ignore_dicts=True): unicode_to_byte(value, ignore_dicts=True) for key, value in data.iteritems()
        }
    
    return data


def load_json_to_string_dict(json_file):
    """Load json (unicode) and further convert to string
    """
    with open(json_file, 'r') as fp:
        unicode_json_dict = json.load(fp, object_pairs_hook=OrderedDict)
        json_dict = unicode_to_byte(unicode_json_dict)
        
    return json_dict
