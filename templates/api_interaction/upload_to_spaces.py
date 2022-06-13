#!/usr/bin/env python

# This sends the upload URL back to nextflow using stdout,
# so don't include other print() calls.
import boto3
import os
# import re
# import csv

region = "fra1"
root = os.environ.get("SPACES_ROOT_DIR")
if not root:
    raise EnvironmentError("Required envvar SPACES_ROOT_DIR not set.")
key = os.environ.get("SPACES_KEY")
if not key:
    raise EnvironmentError("Required envvar SPACES_KEY not set.")
secret = os.environ.get("SPACES_SECRET")
if not secret:
    raise EnvironmentError("Required envvar SPACES_SECRET not set.")

session = boto3.session.Session()
client = session.client(
    's3',
    region_name=region,
    endpoint_url=f'https://{region}.digitaloceanspaces.com',
    aws_access_key_id=key,
    aws_secret_access_key=secret
)

# Metadata seems to give us an AccessDenied error
# try:
#     with open('${qualifyr_report[1]}', "r") as q_report:
#         reader = csv.DictReader(q_report, delimiter='\\t')
#         report = next(reader)
#         metadata = {re.sub("[^0-9a-zA-Z]", "", k): v for k, v in report.items()}
# except StopIteration:
#     metadata = {}

with open('${assembled_genome}', 'rb') as f:
    client.put_object(
        Body=f,
        ContentLength=os.path.getsize('${assembled_genome}'),
        Bucket=root,
        Key=os.path.basename('${assembled_genome}'),
        ACL='public-read',
        # Metadata=metadata
    )
print(f"{region}.digitaloceanspaces.com/{root}/{os.path.basename('${assembled_genome}')}", end='')
