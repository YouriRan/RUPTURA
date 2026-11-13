#!/usr/bin/env python3
"""Collect execution status, CSS status and model metrics from a HAL run."""
import sys
from hal_examples import main

if __name__ == '__main__':
    raise SystemExit(main(['collect'] + sys.argv[1:]))
