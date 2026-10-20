"""
Verification tests for basketball_bounce_sim.

Run from the repository root with:
    python -m pytest -q

Each test checks one physical or numerical property the model must satisfy.
See docs/DESIGN.md, Section 7 (Verification and validation).
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import basketball_bounce_sim as bbs  # noqa: E402


def params(**kw):
    p = bbs.default_params()
    p.update(kw)
    return p


def test_version_string():
    major, minor, patch = (int(x) for x in bbs.__version__.split("."))
    assert major >= 1


def test_unit_conversions():
    assert bbs.KPA_PER_PSI == pytest.approx(6.894757, rel=1e-6)
    assert bbs.FT_TO_M == 0.3048
    assert bbs.IN_TO_M == 0.0254


def test_geometry_limits():
    R = 0.12
    assert bbs.ball_volume(0.0, R) == pytest.approx(4 / 3 * np.pi * R**3)
    assert bbs.contact_area(0.0, R) == 0.0
    # dV/dx = -A(x): volume lost equals the flattened patch area swept
    x, dx = 0.02, 1e-7
    dVdx = (bbs.ball_volume(x + dx, R) - bbs.ball_volume(x - dx, R)) / (2 * dx)
    assert dVdx == pytest.approx(-bbs.contact_area(x, R), rel=1e-6)


def test_gas_energy_is_integral_of_gas_force():
    R, Pg = 0.12, 8 * bbs.PSI_TO_PA
    xs = np.linspace(0, 0.03, 4001)
    work = np.trapezoid(bbs.gas_force(xs, R, Pg), xs) if hasattr(np, "trapezoid") \
        else np.trapz(bbs.gas_force(xs, R, Pg), xs)
    assert work == pytest.approx(bbs.gas_energy(0.03, R, Pg), rel=1e-5)


def test_gent_relation_monotonic_and_plausible():
    E = [bbs.shore_A_to_youngs_modulus(s) for s in (45, 60, 75, 85)]
    assert all(np.diff(E) > 0)
    assert 1e6 < E[0] < 3e6          # ~45 Shore A is a soft rubber, ~2 MPa
    assert 9e6 < E[-1] < 15e6        # ~85 Shore A ~ 10-13 MPa


def test_energy_balance():
    res = bbs.simulate_drop(params())
    assert res["max_energy_error_J"] < 1e-6


def test_elastic_ball_has_cor_one():
    res = bbs.simulate_drop(params(damping_lambda=0.0, air_drag=False))
    assert res["COR"] == pytest.approx(1.0, abs=1e-5)


def test_small_drop_matches_linear_spring_theory():
    # For small deflections the ball is a linear spring k = 2 pi R P + k_shell
    # [Tadrist & Darbois Texier 2017, eq. 2.12].  Including gravity, the exact
    # contact time of a linear spring-mass hitting at speed v0 is
    #     t_c = (pi + 2 atan(g / (v0 w))) / w,   w = sqrt(k/m)
    # (gravity's static sag mg/k keeps the ball on the floor a little longer).
    res = bbs.simulate_drop(params(drop_height_ft=0.5, damping_lambda=0.0, air_drag=False))
    b = res["ball"]
    w = np.sqrt((b["k_gas0"] + b["k_shell"]) / b["m"])
    t_exact = (np.pi + 2 * np.arctan(bbs.G / (res["v0"] * w))) / w * 1e3
    assert res["contact_time_ms"] == pytest.approx(t_exact, rel=0.01)


def test_contact_time_in_published_range():
    # Measured basketball contact times: ~11 ms (Tadrist 2017), ~15 ms (Cross)
    res = bbs.simulate_drop(params())
    assert 9.0 < res["contact_time_ms"] < 16.0


def test_higher_pressure_shorter_contact_and_higher_force():
    lo = bbs.simulate_drop(params(pressure_psi=7.5, damping_lambda=2000))
    hi = bbs.simulate_drop(params(pressure_psi=8.5, damping_lambda=2000))
    assert hi["contact_time_ms"] < lo["contact_time_ms"]
    assert hi["peak_force_N"] > lo["peak_force_N"]
    assert hi["COR"] > lo["COR"]     # matches Njock Libii (2011) experiments


def test_calibration_hits_rule_band():
    p = params()
    lam = bbs.calibrate_damping(p)
    h = bbs.rule_test_rebound(p, lam, bbs.DEFAULT_TEST)
    lo, hi = bbs.RULE_TESTS[bbs.DEFAULT_TEST]["rebound_m"]
    assert h == pytest.approx((lo + hi) / 2, abs=0.002)
    assert all(ok for _, ok in bbs.rule_check(p, lam).values())


def test_floor_never_pulls():
    res = bbs.simulate_drop(params(damping_lambda=6000))
    assert (res["timeseries"]["F_floor_N"] >= 0).all()


def test_csv_round_trip(tmp_path):
    res = bbs.simulate_drop(params())
    files = bbs.export_csv(res, str(tmp_path / "out"))
    df = pd.read_csv(files[0], comment="#")
    assert len(df) == len(res["timeseries"])
    assert df["F_floor_N"].max() == pytest.approx(res["peak_force_N"], rel=1e-6)
    with open(files[0]) as f:
        assert f"# tool_version, {bbs.__version__}" in f.read()


# ---------------------------------------------------------------------------
# Initial-velocity option (added in v1.2.0)
# ---------------------------------------------------------------------------
def _ball(drag):
    return bbs.build_ball(params(damping_lambda=0.0, air_drag=drag))


@pytest.mark.parametrize("v_init", [-6.0, -1.0, 0.0, 1.0, 6.0])
def test_no_drag_impact_speed_is_energy_conservation(v_init):
    # Without drag, v0 = sqrt(v_init^2 + 2 g h) whatever the direction.
    h = 2.0
    assert bbs.impact_speed(h, _ball(False), v_init) == pytest.approx(
        np.sqrt(v_init**2 + 2 * bbs.G * h), rel=1e-12)


def test_drag_from_rest_matches_closed_form():
    b, h = _ball(True), 3.0
    vt = bbs.terminal_velocity(b)
    assert bbs.impact_speed(h, b, 0.0) == pytest.approx(
        vt * np.sqrt(1 - np.exp(-2 * bbs.G * h / vt**2)), rel=1e-12)


def test_drag_upward_throw_equals_drop_from_apex():
    b, h, v = _ball(True), 1.0, 5.0
    apex = h + bbs.rise_height(v, b)
    assert bbs.impact_speed(h, b, v) == pytest.approx(bbs.impact_speed(apex, b, 0.0), rel=1e-12)


def test_drag_down_throw_faster_than_drop_and_bounded_by_terminal_speed():
    b, h = _ball(True), 3.0
    vt = bbs.terminal_velocity(b)
    assert bbs.impact_speed(h, b, -5.0) > bbs.impact_speed(h, b, 0.0)
    # thrown down faster than terminal speed: drag slows it toward vt
    v_fast = bbs.impact_speed(h, b, -(vt + 5.0))
    assert vt < v_fast < vt + 5.0


def test_drag_makes_up_throw_slower_than_down_throw():
    b = _ball(True)
    assert bbs.impact_speed(2.0, b, +5.0) < bbs.impact_speed(2.0, b, -5.0)


def test_simulate_drop_with_initial_velocity_energy_and_metadata(tmp_path):
    res = bbs.simulate_drop(params(v_init_m_s=-8.0))
    assert res["max_energy_error_J"] < 1e-6
    assert res["equiv_drop_height_m"] > res["drop_height_m"]
    files = bbs.export_csv(res, str(tmp_path / "throw"))
    with open(files[0]) as f:
        assert "# v_init_m_s_positive_up, -8.0" in f.read()
