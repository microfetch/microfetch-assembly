#! /usr/bin/env python
import sys

import requests
import json
import csv
import re
import os

with open("${api_response}", "rb") as f:
    response = json.load(f)

data = {
    'assembly_result': 'success',
    'assembled_genome_url': """${spaces_url}""".replace("""\n""", '')  # avoid EOL issue
}

try:
    with open("${qualifyr_report[1]}", "r") as q_report:
        reader = csv.DictReader(q_report, delimiter="\t")
        data['qualifyr_report'] = json.dumps(next(reader))
except StopIteration:
    pass

upload_url = response['upload_url']
url = os.path.join('${api_url}', re.sub('^/', '', upload_url))

print(f"PUT {url}")
print(data)

r = requests.request('PUT', url, data=data)

if r.status_code != 204:
    print(r.json()['error'], file=sys.stderr)
    raise ConnectionError(f"API call failed (Status code {r.status_code})")
