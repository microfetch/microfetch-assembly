#! /usr/bin/env python

import boto3
import os
import requests
import re
import json
import csv
import sys
import logging

region = "fra1"
root = os.environ.get("SPACES_ROOT_DIR")
key = os.environ.get("SPACES_KEY")
secret = os.environ.get("SPACES_SECRET")
out_dir = sys.argv[1]
api_url = sys.argv[2]
error_report_file = sys.argv[3]

logger = logging.getLogger(__file__)
logger.addHandler(logging.FileHandler(f"{out_dir}/api_interaction/error_upload.log"))
logger.setLevel(logging.DEBUG)

logger.info(f"{sys.argv[0]} {sys.argv[1]} {sys.argv[2]} {sys.argv[3]}")

session = boto3.session.Session()
logger.debug(session)
client = session.client(
    's3',
    region_name=region,
    endpoint_url=f'https://{region}.digitaloceanspaces.com',
    aws_access_key_id=key,
    aws_secret_access_key=secret
)
logger.debug(client)

# Try to discover the API sample id from the output directory
if out_dir:
    file = f"{out_dir}/api_interaction/api_response.json"
    logger.info(f"{file} exists={os.path.exists(file)}")
    if os.path.exists(file):
        with open(file, 'r') as f:
            j = json.load(f)

logger.debug(json.dumps(j))

if j and 'id' in j.keys():
    client.put_object(
        Bucket=root,
        Key=os.path.basename(f'{j["id"]}.txt'),
        Body=error_report_file,
        ACL='public-read',
        Metadata={
            'x-amz-meta-error-report': 'true'
        }
    )

    data = {
        'assembly_result': 'fail',
        'assembly_error_report_url': f"{region}.digitaloceanspaces.com/{root}/{j['id']}.txt"
    }

    try:
        report_file = f"{out_dir}/quality_reports/qualifyr_report.tsv"
        with open(report_file, "r") as q_report:
            reader = csv.DictReader(q_report, delimiter="\t")
            data['qualifyr_report'] = json.dumps(next(reader))
    except StopIteration:
        pass

    logger.debug(data)

    # Upload report to API
    upload_url = j['upload_url']
    url = os.path.join(api_url, re.sub('^/', '', upload_url))
    r = requests.request('PUT', url, data=data)

    if r.status_code != 204:
        logger.error(f"API call failed (Status code {r.status_code})")
        raise ConnectionError(f"API call failed (Status code {r.status_code})")

else:
    logger.error("No json found")
