#!/usr/bin/env python3
"""Download Synapse file entities through the REST API without synapseclient.

The personal access token is read from ~/.synapseConfig (first `authtoken`
line) or the SYNAPSE_AUTH_TOKEN environment variable and is never printed.
Existing non-empty files are not overwritten.
"""
import argparse
import hashlib
import json
import os
import pathlib
import re
import sys
import urllib.request

API = "https://repo-prod.prod.sagebase.org"


def read_token():
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if token:
        return token.strip()
    cfg = pathlib.Path.home() / ".synapseConfig"
    for line in cfg.read_text().splitlines():
        m = re.match(r"\s*authtoken\s*=\s*(\S+)", line)
        if m:
            return m.group(1)
    sys.exit("No Synapse token found in SYNAPSE_AUTH_TOKEN or ~/.synapseConfig")


def get_json(url, token):
    req = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
    with urllib.request.urlopen(req) as resp:
        return json.load(resp)


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("ids", nargs="+", help="Synapse file entity IDs")
    args = ap.parse_args()
    token = read_token()
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for syn_id in args.ids:
        ent = get_json(f"{API}/repo/v1/entity/{syn_id}", token)
        handle = ent["dataFileHandleId"]
        dest = out / ent["name"]
        expected_md5 = get_json(
            f"{API}/repo/v1/entity/{syn_id}/filehandles", token
        )["list"]
        expected_md5 = next((fh.get("contentMd5") for fh in expected_md5 if fh["id"] == handle), None)
        if not (dest.exists() and dest.stat().st_size > 0):
            url = (
                f"{API}/file/v1/file/{handle}?fileAssociateType=FileEntity"
                f"&fileAssociateId={syn_id}&redirect=false"
            )
            req = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
            with urllib.request.urlopen(req) as resp:
                presigned = resp.read().decode().strip()
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(presigned) as resp, open(tmp, "wb") as fh:
                while chunk := resp.read(1 << 22):
                    fh.write(chunk)
            tmp.rename(dest)
        local_md5 = md5sum(dest)
        status = "OK" if expected_md5 in (None, local_md5) else "MISMATCH"
        print(f"{syn_id}\t{dest.name}\t{dest.stat().st_size}\t{local_md5}\t{expected_md5}\t{status}")
        if status != "OK":
            sys.exit(f"MD5 mismatch for {dest}")


if __name__ == "__main__":
    main()
