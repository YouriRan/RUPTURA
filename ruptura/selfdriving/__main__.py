"""Runs one experiment directory of a campaign: python -m ruptura.selfdriving <directory>."""

import sys

from .psa import run_experiment

run_experiment(sys.argv[1])
