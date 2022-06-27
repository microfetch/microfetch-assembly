#!/usr/bin/env python

# This sends the upload URL back to nextflow using stdout,
# so don't include other print() calls.
import boto3
import os
import hashlib
import gzip
import shutil

# Calculate and store SHA1 digest
with open('${assembled_genome}', 'rb') as f:
    sha1 = hashlib.sha1(f.read())

with open("assembled_genome_sha1.txt", "w+") as f:
    f.write(sha1.hexdigest())

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

zipped_file = f"{os.path.basename('${assembled_genome}')}.gz"

with open('${assembled_genome}', 'rb') as f_in:
    with gzip.open(zipped_file, 'w+b') as f_out:
        shutil.copyfileobj(f_in, f_out)

with open(zipped_file, 'rb') as f:
    client.put_object(
        Body=f,
        ContentLength=os.path.getsize(zipped_file),
        Bucket=root,
        Key=zipped_file,
        ACL='public-read',
        # Metadata=metadata
        Metadata={
            'sha1': sha1.hexdigest(),
            'content-type': 'application/gzip'
        }
    )

# Save upload file path
with open("assembled_genome_url.txt", "w+") as f:
    f.write(f"{region}.digitaloceanspaces.com/{root}/{zipped_file}")
