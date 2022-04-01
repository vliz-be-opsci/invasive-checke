#!/usr/bin/env python

"""Tests for `invasive_checker` package."""

import pytest
from invasive_checker import Aphia_Checker

@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    checker = Aphia_Checker()

@pytest.fixture
def invasive_coral():
    aphia_id = 132818

@pytest.fixture
def native_cod():
    aphia_id = 126436

@pytest.fixture
def belgium_eez():
    lon = 2.5
    lat = 51.5

@pytest.fixture
def null_island():
    #The one MRGID runs along the equator. Move null island slightly south...
    lon = -0.01
    lat = -0.01

def test_invasive(invasive_coral, belgium_eez):
    """Sample pytest test function with the pytest fixture as an argument."""
    aphia_checker = Aphia_Checker()
    aa, bb = aphia_checker.check_aphia(belgium_eez.lon, belgium_eez.lat, invasive_coral.aphia_id)

    assert aa['aphia_id']==132818
    assert aa['distance [km] to nearest introduced location']==0.0
    assert aa['nearest introduced MRGID'] == [21912]

def test_native(native_cod, belgium_eez):
    """Sample pytest test function with the pytest fixture as an argument."""
    aphia_checker = Aphia_Checker()
    aa, bb = aphia_checker.check_aphia(belgium_eez.lon, belgium_eez.lat, native_cod.aphia_id)

    assert aa['aphia_id']==126436
    assert aa['sample location within <buffer> of aphia distribution'] is True
    assert aa['distance [deg] to nearest introduced location'] == "No known 'introduced' locations"
