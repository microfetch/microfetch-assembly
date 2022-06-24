#!/usr/bin/env python
import sys
from requests import request
import os
import json

url = os.path.join('${api_url}', 'records/awaiting_assembly/')
print(f"GET {url}")
r = request('GET', url)

if r.status_code == 204:
    # Server says there's nothing for us
    print("Server responded with 204: No content. There are no records to process.")
else:
    if r.status_code != 200:
        raise ConnectionError(f"API call failed (Status code {r.status_code})")

    try:
        j = r.json()
        if len(j['results']) == 0:
            raise ConnectionError(f"API has no results for processing")
        for record in j["results"]:
            print(record)
            try:
                accept_url = record['action_links']['register_assembly_attempt']
                print(f"GET {accept_url}")
                r = request('GET', accept_url)
            except BaseException as e:
                print(f"Error accepting {record['id']} -- {e}", file=sys.stderr)

            if r.status_code == 200:
                print(f"Accepted record {record['id']} for assembly")
                # Save sample id for future reference
                os.environ['API_SAMPLE_ID'] = record['id']
                os.environ['API_UPLOAD_URL'] = record['action_links']['report_assembly_result']

                # Patch record with the taxon's post-assembly-filters
                taxon_response = request('GET', record['taxon'])
                if taxon_response.status_code == 200:
                    taxon_json = taxon_response.json()
                    record['post_assembly_filters'] = taxon_json['post_assembly_filters']

                    with open("api_response.json", "w+") as f:
                        json.dump(record, f)

                    break
                else:
                    print((
                        "Unable to download post assembly filters for record. "
                        "The record will be marked as undergiong assembly but not assembled."
                    ))

    except ValueError:
        raise ValueError("Unable to interpret API response as JSON.")
