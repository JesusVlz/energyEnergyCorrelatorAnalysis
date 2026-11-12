"""Run-3 ppref-data CRAB entry point; use submitRun3.py for the full campaign."""

import os
import sys
from pathlib import Path

# CRAB may evaluate this file with a different sys.path than a normal script.
_config_path = Path(globals().get('__file__') or os.getcwd()).resolve()
_config_dir = _config_path.parent if _config_path.is_file() else _config_path
sys.path.insert(0, str(_config_dir))

from run3_crab import create_config

config = create_config('ppref-data')
