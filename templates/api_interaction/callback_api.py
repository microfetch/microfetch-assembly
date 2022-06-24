#! /usr/bin/env python
import sys

import requests
import json
import csv

with open("${api_response}", "rb") as f:
    response = json.load(f)

with open("${upload_url_file}", "r") as f:
    upload_url = f.read()

with open("${upload_sha_file}", "r") as f:
    upload_sha = f.read()

data = {
    'assembly_result': 'success',
    'assembled_genome_url': upload_url,
    'assembled_genome_sha1': upload_sha
}

try:
    with open("${qualifyr_report[1]}", "r") as q_report:
        reader = csv.DictReader(q_report, delimiter="\t")
        data['qualifyr_report'] = json.dumps(next(reader))
except StopIteration:
    pass

url = response['action_links']['report_assembly_result']

print(f"POST {url}")
print(data)

r = requests.request('POST', url, data=data)

if r.status_code != 200:
    print(r.text, file=sys.stderr)
    raise ConnectionError(f"API call failed (Status code {r.status_code})")
