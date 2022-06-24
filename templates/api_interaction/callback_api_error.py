#! /usr/bin/env python

import boto3
import os
import requests
import re
import json
import csv
import sys
import logging


def first_match(s: str, txt: str) -> str:
    try:
        m = re.findall(s, txt)
        return m[0]
    except IndexError:
        return f"Failed to parse error file for RegEx: {s}"


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

spaces_url = ""
upload_sha = ""

# Try to discover the API sample id from the output directory
if out_dir:
    file = f"{out_dir}/api_interaction/api_response.json"
    logger.info(f"{file} exists={os.path.exists(file)}")
    if os.path.exists(file):
        with open(file, 'r') as f:
            j = json.load(f)

    try:
        with open(f"{out_dir}/api_interaction/assembled_genome_url.txt", "r") as f:
            spaces_url = f.read()
    except FileNotFoundError:
        pass

    try:
        with open(f"{out_dir}/api_interaction/assembled_genome_sha.txt", "r") as f:
            upload_sha = f.read()
    except FileNotFoundError:
        pass

logger.debug(json.dumps(j))

if j and 'id' in j.keys():
    try:
        logger.info(f"Uploading {error_report_file}")
        with open(error_report_file, 'rb') as error_report:
            client.put_object(
                Body=error_report.read(),
                Bucket=root,
                Key=os.path.basename(f'{j["id"]}.txt'),
                ACL='public-read',
                Metadata={
                    'x-amz-meta-error-report': 'true'
                }
            )
        upload_url = f"{region}.digitaloceanspaces.com/{root}/{j['id']}.txt"
    except BaseException as e:
        logger.info(f"boto3 exception: {e}")
        upload_url = ""

    # Parse error_report_file for individual fields to send to API
    logger.info(f"Parsing {error_report_file}")
    with open(error_report_file, 'r') as error_report:
        text = error_report.read()
        data = {
            'assembly_result': 'fail',
            'assembled_genome_url': spaces_url,
            'assembled_genome_sha1': upload_sha,
            'assembly_error_report_url': upload_url,
            'assembly_error_process': first_match("Caused by:\\s+Process `([^`]*)`", text),
            'assembly_error_exit_code': first_match("Command exit status:\\s+(.*)", text),
            'assembly_error_stdout': first_match("Command output:\n([\\s\\S]*?)\n(?:Command error:|Work dir:)", text),
            'assembly_error_stderr': first_match("Command error:\n([\\s\\S]*?)\nWork dir:", text)
        }

    try:
        report_file = f"{out_dir}/quality_reports/qualifyr_report.tsv"
        logger.info(f"Parsing {report_file}")
        with open(report_file, "r") as q_report:
            reader = csv.DictReader(q_report, delimiter="\t")
            data['qualifyr_report'] = json.dumps(next(reader))
    except FileNotFoundError:
        logger.info('No qualifyr report found.')
    except StopIteration:
        logger.info('Qualifyr Report is empty.')

    logger.debug(data)

    # Upload report to API
    url = j['action_links']['report_assembly_result']
    logger.info(f"POST {url}")
    r = requests.request('POST', url, data=data)

    if r.status_code != 204:
        logger.error(f"API call failed (Status code {r.status_code})")
        raise ConnectionError(f"API call failed (Status code {r.status_code})")

else:
    logger.error("No json found")
