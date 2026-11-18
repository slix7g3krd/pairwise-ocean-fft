import pytest
from render import ocean_fft as mod

CASES = [
    ("wind_base", "m.cascade_wind(10.0, 0) == approx(10.0)"),
    ("wind_higher", "m.cascade_wind(10.0, 2) > m.cascade_wind(10.0, 0)"),
    ("spectrum_positive", "m.phillips(1.0, 0.0, 0.0, 10.0) > 0.0"),
    ("spectrum_cross_wind", "m.phillips(0.0, 1.0, 0.0, 10.0) < m.phillips(1.0, 0.0, 0.0, 10.0)"),
    ("normal_unit", "sum(v * v for v in m.fold_normal(1.0, 0.0, 1.0)) == approx(1.0)"),
    ("normal_orientation", "m.fold_normal(1.0, 0.0, 1.0)[0] > 0.0"),
    ("foam_resolution_free", "m.foam_threshold(0.5, 64) == approx(m.foam_threshold(0.5, 1024))"),
    ("foam_monotonic", "m.foam_threshold(0.2, 64) > m.foam_threshold(0.8, 64)"),
    ("clipmap_ring0", "m.clipmap_uv(5.0, 0, 4.0) == approx(1.0)"),
    ("clipmap_ring1", "m.clipmap_uv(5.0, 1, 4.0) == approx(5.0)"),
    ("buoyancy_physical", "m.buoyancy_height(1.0, 0.1) == approx(1.0)"),
    ("buoyancy_no_visual_bias", "m.buoyancy_height(1.0, 0.1) == approx(m.buoyancy_height(1.0, 0.5))"),
    ("underwater_true", "m.is_underwater(0.0, -1.0) is True"),
    ("underwater_false", "m.is_underwater(1.0, 0.0) is False"),
    ("spectrum_wind_strength", "m.phillips(2.0, 0.0, 0.0, 20.0) > m.phillips(2.0, 0.0, 0.0, 10.0)")
]

@pytest.mark.parametrize("name,expr", CASES, ids=[c[0] for c in CASES])
def test_contract(name, expr):
    assert eval(expr, {"m": mod, "approx": pytest.approx})
