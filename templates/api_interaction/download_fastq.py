#! /usr/bin/env python

from json import loads
import os

with open("${json_file}", "r") as f:
    j = loads(f.read())

os.mkdir("api_input_files")

sources = j['fastq_ftp'].split(';')
for source in sources:
    print(f"Downloading {source}")
    os.system(f"wget --tries=5 --wait=2 --output-document api_input_files/{os.path.basename(source)} ftp://{source}")
