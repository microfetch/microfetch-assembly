#!/usr/bin/env python
from requests import request
import os
import re

url = os.path.join('${api_url}', 'api/request_assembly_candidate/')
print(f"GET {url}")
r = request('GET', url)

if r.status_code == 204:
    # Server says there's nothing for us
    pass
if r.status_code != 200:
    raise ConnectionError(f"API call failed (Status code {r.status_code})")

try:
    j = r.json()

    accept_url = j['accept_url']
    url = os.path.join('${api_url}', re.sub('^/', '', accept_url))
    if request('GET', url).status_code != 204:
        raise ConnectionError(f"Unable to accept assembly candidate via GET {url}.")

    # Save sample id for future reference
    os.environ['API_SAMPLE_ID'] = j['id']
    os.environ['API_UPLOAD_URL'] = j['upload_url']

    with open("api_response.json", "w+") as f:
        f.write(r.text)

except ValueError:
    raise ValueError("Unable to interpret API response as JSON.")
