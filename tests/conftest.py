""" Define global fixtures here """

from __future__ import (print_function, absolute_import)

import py
import os
import shutil
from ConfigParser import RawConfigParser
import pytest

from line_classifier.lfs_ews.equivalent_width import EquivalentWidthAssigner


collect_ignore = ["setup.py"]

@pytest.fixture(scope="session")
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
    return datadir.join("config").join("Line_flux_limit_5_sigma_baseline.dat").strpath

@pytest.fixture
def oii_ew_assigner(config):
    return EquivalentWidthAssigner.from_config(config, "OII_EW")

@pytest.fixture
def lae_ew_assigner(config):
    return EquivalentWidthAssigner.from_config(config, "LAE_EW")


@pytest.fixture(scope="session")
def tmpdir_with_config(tmpdir_factory, datadir):
    """
    A temporary directory with the relevant
    configuration files inside
    """
    dirname = tmpdir_factory.mktemp("tmprun")
 
    # XXX todo, replace with glob loop??
    shutil.copy(datadir.join("config").join("Line_flux_limit_5_sigma_baseline.dat").strpath, dirname.strpath)
    shutil.copy(datadir.join("config").join("lae_log_ew_obs_0_600_glim25_redo.fits").strpath, dirname.strpath)
    shutil.copy(datadir.join("config").join("oii_log_ew_obs_0_600_glim25_redo.fits").strpath, dirname.strpath)

    return dirname.strpath

