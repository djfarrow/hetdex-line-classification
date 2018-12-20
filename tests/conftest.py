""" Define global fixtures here """

from __future__ import (print_function, absolute_import)

import py
import os
from ConfigParser import RawConfigParser
import pytest

from line_classifier.lfs_ews.equivalent_width import EquivalentWidthAssigner


collect_ignore = ["setup.py"]

@pytest.fixture
def datadir():
    """ Return a py.path.local object for the test data directory"""
    return py.path.local(os.path.dirname(__file__)).join('data')

@pytest.fixture
def config(datadir):
    config_fn = datadir.join("universe.cfg").strpath
    config = RawConfigParser()
    config.read(config_fn)

    return config

@pytest.fixture
def flim_file(datadir):
    return datadir.join("Line_flux_limit_5_sigma_baseline.dat").strpath

@pytest.fixture
def oii_ew_assigner(config):
    return EquivalentWidthAssigner.from_config(config, "OII_EW")

@pytest.fixture
def lae_ew_assigner(config):
    return EquivalentWidthAssigner.from_config(config, "LAE_EW")



