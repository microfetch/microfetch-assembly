#! /usr/bin/env python

# This program filters genomes based on QC metrics.
# quality report from Assembly pipeline
# The results are saved to a tsv file.

import pandas as pd
import json
import pathlib


def apply_filter(df: pd.DataFrame, metric: str, f_type: str, f_value: [int, float, str]) -> pd.DataFrame:
    if f_type == "min_max":
        metric_min, metric_max = f_value
        return df[(df[metric] >= metric_min) & (df[metric] <= metric_max)]
    elif f_type == "max":
        return df[df[metric] <= f_value]
    elif f_type == "min":
        return df[df[metric] >= f_value]
    elif f_type == "contains":
        f_value = [x.lower() for x in f_value]
        return df[df[metric].str.contains('|'.join(f_value), case=False)]
    raise ValueError("Only 'min_max', 'min', 'max', and 'contains' can be used as f_type")


def get_qc_filters(qc_file: str) -> list:
    with open(qc_file, 'r') as f:
        data = json.load(f)
    filters = data["filters"]
    result = []
    for f in filters:
        temp = []
        temp.append(f)
        temp.append(filters[f]["type"])
        value = filters[f]["value"] if isinstance(filters[f]["value"], int) else tuple(filters[f]["value"])
        temp.append(value)
        result.append(temp)
    return result

print("QC file:\t${qc_file[1]}")  # test_output/quality_reports/qualifyr_report.tsv

# Set of QC filters -- should be provided by the API server
print("QC filters:\t${api_response}")  # a JSON file that varys by organism

# Output file
print("OUT file:\tpost_assembly_filter.tsv")

# Load quality_summary table
qc_file = "${qc_file[1]}"
df = pd.read_csv(qc_file, "\t")
with open("${api_response}", "rb") as f:
    content = json.load(f)
filters = content['post_assembly_filters']

out_file = "post_assembly_filter.tsv"
if not filters:
    print("No post_assembly_filters to apply.")
    pathlib.Path(out_file).touch(exist_ok=True)
else:
    filtering = get_qc_filters(filters)
    print("# Original size: " + str(len(df)))
    for f in filtering:
        df = apply_filter(df, f[0], f[1], f[2])
        print("# Size after '" + f[0] + "' filter: " + str(len(df)))
    df.to_csv(out_file, sep="\t", index=False)
