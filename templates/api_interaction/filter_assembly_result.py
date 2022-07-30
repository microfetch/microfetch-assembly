#! /usr/bin/env python

# This program filters genomes based on QC metrics.
# quality report from Assembly pipeline

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


print("QC file:\t${qc_file[1]}")  # test_output/quality_reports/qualifyr_report.tsv

# Set of QC filters -- should be provided by the API server
print("API Response file:\t${api_response}")  # a JSON file that varys by organism

# Load quality_summary table
qc_file = "${qc_file[1]}"
df = pd.read_csv(qc_file, delimiter="\t")
with open("${api_response}", "rb") as f:
    content = json.load(f)
filters = content['post_assembly_filters']
record_id = content['id']

if not filters or len(filters) == 0:
    print("Filters found:\tNone")
else:
    print(f"Filters found:\t{len(filters)}")

    # List comprehension to map filters to a dictionary of filter names: True/False (True = passed)
    results = {k: len(apply_filter(df, k, v['type'], v['value'])) > 0 for k, v in filters.items()}
    failed_filters = []
    for k, v in results.items():
        print(f"Filter '{k}':\t{'pass' if v else 'FAIL'}")
        if not v:
            failed_filters.append(k)

    if len(failed_filters):
        fails = ", ".join([f"'{f}'" for f in failed_filters])
        raise ValueError(f"Post assembly filters failed: {fails}.")
