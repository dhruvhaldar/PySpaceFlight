import numpy as np
import pytest
import sys
import os

# Add src to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from orbit import G, M_EARTH, MU, get_orbital_elements, hohmann_transfer
from rocket import RocketStage, G0

def test_orbital_period():
    # Test for circular LEO at 500km altitude
    r = 6371000 + 500000
    v = np.sqrt(MU / r)
    
    r_vec = np.array([r, 0.0])
    v_vec = np.array([0.0, v])
    
    elements = get_orbital_elements(r_vec, v_vec)
    
    a = elements['a']
    expected_period = 2 * np.pi * np.sqrt(a**3 / MU)
    
    # Check if calculated period matches formula
    assert np.isclose(elements['period'], expected_period, rtol=1e-5)
    
    # Check if period is roughly 94.6 minutes (5677 s)
    assert 5600 < elements['period'] < 5750

def test_hohmann_transfer_leo_geo():
    # LEO: 300km, GEO: 35786km
    r1 = 6371000 + 300000
    r2 = 6371000 + 35786000
    
    dv1, dv2, total_dv, t_transfer = hohmann_transfer(r1, r2)
    
    # Expected values roughly:
    # dv1 ~ 2.4 km/s
    # dv2 ~ 1.5 km/s
    # total ~ 3.9 km/s
    
    assert 2000 < dv1 < 2600
    assert 1200 < dv2 < 1700
    assert 3500 < total_dv < 4500

def test_tsiolkovsky():
    # Test rocket equation
    dry_mass = 1000
    fuel_mass = 9000
    isp = 300
    thrust = 10000
    
    stage = RocketStage(dry_mass, fuel_mass, thrust, isp, 1.0)
    
    initial_mass = dry_mass + fuel_mass
    final_mass = dry_mass
    
    expected_dv = isp * G0 * np.log(initial_mass / final_mass)
    
    # Simulate burn
    # Note: Our simulation uses drag and gravity, so this simple check only works in vacuum/zero-g
    # But we can check the mass flow rate logic
    
    dt = 1.0
    _, mdot = stage.burn(dt)
    
    expected_mdot = thrust / (isp * G0)
    assert np.isclose(mdot, expected_mdot, rtol=1e-5)
