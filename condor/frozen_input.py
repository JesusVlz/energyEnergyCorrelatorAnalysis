"""Check frozen source metadata before and after a worker stages an input.

XRootD exposes modification times to whole seconds. Local test paths retain
nanosecond comparisons; neither metadata check is a content checksum.
"""
import datetime as dt
import os
from pathlib import Path
import re
import subprocess
from urllib.parse import urlsplit


def source_metadata(url, local_test=False):
    parsed = urlsplit(url)
    if not parsed.scheme:
        if not local_test or not Path(url).is_absolute():
            raise ValueError('Frozen input requires an absolute XRootD URL')
        stat = Path(url).stat()
        return {'size_bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
    if (parsed.scheme != 'root' or not parsed.netloc or not parsed.path.startswith('/')
            or parsed.query or parsed.fragment):
        raise ValueError('Frozen input requires an absolute XRootD URL')
    result = subprocess.run(['xrdfs', 'root://' + parsed.netloc, 'stat',
                             '/' + parsed.path.lstrip('/')], check=True,
                            capture_output=True, text=True, timeout=60,
                            env={**os.environ, 'TZ': 'UTC'})
    size = re.search(r'^Size:\s*(\d+)\s*$', result.stdout, re.M)
    modified = re.search(r'^MTime:\s*(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\s*$', result.stdout, re.M)
    if size is None or modified is None:
        raise ValueError('XRootD stat lacks readable Size/MTime: ' + url)
    timestamp = dt.datetime.strptime(modified.group(1), '%Y-%m-%d %H:%M:%S').replace(tzinfo=dt.timezone.utc)
    return {'size_bytes': int(size.group(1)), 'mtime_seconds': int(timestamp.timestamp())}


def check_source_identity(url, expected, local_test=False):
    if (type(expected.get('size_bytes')) is not int or expected['size_bytes'] < 0
            or type(expected.get('mtime_ns')) is not int or expected['mtime_ns'] < 0):
        raise ValueError('Frozen source requires size_bytes and mtime_ns; regenerate the inventory')
    actual = source_metadata(url, local_test)
    if actual['size_bytes'] != expected['size_bytes']:
        raise ValueError('Input size changed since snapshot: ' + url)
    same_time = (actual['mtime_ns'] == expected['mtime_ns'] if 'mtime_ns' in actual else
                 actual['mtime_seconds'] == expected['mtime_ns'] // 1_000_000_000)
    if not same_time:
        raise ValueError('Input modification time changed since snapshot: ' + url)
    return actual


def stage_frozen_input(url, destination, expected, stage, local_test=False):
    before = check_source_identity(url, expected, local_test)
    stage()
    after = check_source_identity(url, expected, local_test)
    if before != after or Path(destination).stat().st_size != expected['size_bytes']:
        raise ValueError('Input changed while staging: ' + url)
    return {'before': before, 'after': after}
