# =============================================================================
#  BASKETBALL DROP-IMPACT SIMULATOR            version 1.2.0  (2026-10-20)
#  Author: Tico Hannan   |   License: BSD-3-Clause (see LICENSE)
#  https://github.com/ticohannan/basketball-bounce-sim
#  Floor contact force, contact duration and energy conversion for a basketball
#  dropped onto a rigid floor, with adjustable inflation and viscoelastic
#  (rubber shell) properties.
#
#  HOW TO RUN (Jupyter on Windows / WinPython 3.12)
#  ------------------------------------------------
#    Option A: paste this whole file into ONE notebook cell and run it.
#    Option B: put this file next to your notebook and run, in a cell:
#                  %run basketball_bounce_sim.py
#    Either way the interactive panel (sliders + plots + CSV button) appears
#    below the cell.
#
#    Required packages (all ship with WinPython): numpy, scipy, matplotlib,
#    pandas, ipywidgets.  No internet access is needed.
#
#    You can also use the physics without the GUI, e.g. in another cell:
#        p   = default_params()          # a dict of ball/drop settings
#        p["drop_height_ft"] = 6
#        p["v_init_m_s"] = -3.0          # optional: thrown DOWN at 3 m/s (+ = up)
#        res = simulate_drop(p)          # run one bounce
#        print(res["contact_time_ms"], res["peak_force_N"])
#        res["timeseries"].head()        # pandas DataFrame of the bounce
#
#  -----------------------------------------------------------------------------
#  THE PHYSICAL PICTURE (read this first)
#  -----------------------------------------------------------------------------
#  1. FLIGHT TO THE FLOOR.  The ball starts with its BOTTOM at height h
#     (same convention the FIBA / NFHS rebound tests use), either
#       * released from rest (a simple drop), or
#       * moving with a starting velocity v_init (POSITIVE = UPWARD, negative =
#         thrown downward), e.g. a dribble or a toss.
#     It reaches the floor at speed v0.  Without air drag, energy conservation
#     gives v0 = sqrt(v_init^2 + 2 g h) whatever the direction of v_init.
#     With drag we use exact quadratic-drag solutions (Section 5).
#
#  2. CONTACT.  We follow ONE variable, x(t) = how far the ball's centre has
#     moved down since first touching the floor (the "deflection" or
#     "indentation").  The bottom of the ball flattens into a circular patch.
#     Newton's 2nd law for the ball (down = positive):
#
#             m * d^2x/dt^2  =  m*g  -  F_floor(x, dx/dt)
#
#     F_floor is the push of the floor on the ball (and, by Newton's 3rd law,
#     exactly the force the FLOOR FEELS - the quantity we want).  It has
#     three parts:
#
#       (a) GAS SPRING   F_gas   = (P_inside - P_atm) * A_contact(x)
#           The pressurised air is the main spring of a basketball.  Squashing
#           the ball shrinks its volume, the air pressure rises (adiabatically,
#           because the bounce is too fast for heat to escape), and that
#           pressure pushes on the flat contact patch.  [Tadrist & Darbois
#           Texier 2017; Stronge & Ashcroft 2007]
#
#       (b) SHELL SPRING F_shell = k_shell * x
#           The rubber carcass also resists bending.  For a thin spherical
#           shell  k_shell = 4 E t^2 / (R sqrt(3(1-nu^2)))  where E is the
#           rubber's Young's modulus (we get it from Shore A hardness with
#           Gent's formula), t is wall thickness, R ball radius.
#           [Lazarus, Florijn & Reis 2012, as used by Tadrist 2017]
#           SURPRISE FOR STUDENTS: for a basketball this is only a few % of the
#           gas stiffness, so rubber hardness barely changes contact time.
#           The air does almost all the work.  Try the slider and see!
#
#       (c) VISCOELASTIC DAMPING  F_damp = lambda * x * dx/dt
#           Real rubber (and the ball's layers, seams, air motion...) turn some
#           energy into heat, so the ball rebounds lower than it was dropped.
#           We use the Hunt-Crossley form (damping proportional to deflection
#           times velocity).  Unlike a plain dashpot (c*dx/dt), it is zero at
#           the instant of touch-down and lift-off, so the floor never "pulls"
#           on the ball.  [Hunt & Crossley 1975]
#           lambda is tuned so the ball passes the official FIBA rebound test.
#
#     Contact ends when x returns to 0 (ball's undeformed shape leaves the
#     floor).  Contact time = that instant.
#
#  3. ENERGY BOOKKEEPING.  At every instant
#        KE_impact + m*g*x  =  KE(t) + E_gas(x) + E_shell(x) + E_dissipated(t)
#     At maximum deflection KE = 0: all the incoming kinetic energy is either
#     stored (potential energy in the compressed air + bent rubber) or already
#     lost as heat.  The code checks this balance numerically.
#
#  -----------------------------------------------------------------------------
#  SOURCES (data and equations)
#  -----------------------------------------------------------------------------
#  [NBA]   NBA Official Rules, Rule No. 1, Sec. II: "The ball shall be an
#          officially approved NBA ball between 7 1/2 and 8 1/2 pounds
#          pressure."  https://official.nba.com/rule-no-1-court-dimensions/
#  [FIBA26] FIBA Official Basketball Rules 2026 - Basketball Equipment v1.0
#          (valid from 1 Oct 2026), Sec. 2: size 7 circumference 750-770 mm,
#          mass 580-620 g; dropped from ~1800 mm (underside of ball) it must
#          rebound 1035-1085 mm (measured to the UNDERSIDE of the ball).
#  [FIBA18] FIBA Official Basketball Rules 2018 - Basketball Equipment,
#          Sec. 7: size 7 749-780 mm, 567-650 g; 1800 mm drop -> rebound
#          1200-1400 mm measured to the TOP of the ball.  (Previous edition;
#          kept as a legacy check and for the wider slider ranges.)
#  [NFHS]  NFHS Basketball Rule 1-12-2: dropped from 6 ft (bottom of ball) it
#          must rebound to 49-54 in (top of ball).
#  [TDT17] L. Tadrist & B. Darbois Texier, "Open questions on the impact of an
#          inflated ball", arXiv:1708.01282 (2017).
#          doi:10.48550/arXiv.1708.01282  Contact-patch geometry
#          A(x), volume V(x), adiabatic gas law, linear-spring contact time
#          t_c = sqrt(pi m / (2 R P0)), shell-vs-gas storage number, and
#          measured basketball data: diameter 24.3 cm, mass 608 g, contact
#          time ~11 ms, COR ~0.78; shell E ~ 10 MPa, thickness ~2.5 mm.
#  [SA07]  W.J. Stronge & A.D.C. Ashcroft, "Oblique impact of inflated balls
#          at large deflections", Int. J. Impact Eng. 34 (2007) 1003-1019.
#          doi:10.1016/j.ijimpeng.2006.04.006
#  [LFR12] A. Lazarus, H.C.B. Florijn & P.M. Reis, "Geometry-induced rigidity
#          in nonspherical pressurized elastic shells", Phys. Rev. Lett. 109
#          (2012) 144301.  doi:10.1103/PhysRevLett.109.144301
#          Thin-shell indentation stiffness.
#  [HC75]  K.H. Hunt & F.R.E. Crossley, "Coefficient of restitution
#          interpreted as damping in vibroimpact", J. Appl. Mech. 42 (1975)
#          440-445.  doi:10.1115/1.3423596
#  [Gent]  A.N. Gent, "On the relation between indentation hardness and
#          Young's modulus", Rubber Chem. Technol. 31 (1958) 896-906.
#          doi:10.5254/1.3542351
#          E[MPa] = 0.0981(56+7.62336 S)/(0.137505(254-2.54 S)), S = Shore A.
#  [Cross] R. Cross, "The bounce of a ball", Am. J. Phys. 67 (1999) 222-227,
#          doi:10.1119/1.19229; and "Bouncing Balls" slides (U. Sydney, 2006,
#          https://www.physics.usyd.edu.au/~cross/BallBounce.htm): basketball
#          contact time on a wood floor is about 15 ms; contact time falls
#          slightly as impact speed rises.
#  [Pat1]  US 2008/0188334 A1 (game ball with urethane gum rubber surface):
#          basketball cover rubber Shore A ~45-85; natural-rubber cover
#          stocks 70-78 Shore A.
#  [Pat2]  US 5,310,178 (basketball with polyurethane cover): bladder ~0.8 mm,
#          windings 0.3-0.7 mm, carcass 1.8-2.0 mm, cover panels ~1.8 mm.
#  [ASEE]  J. Njock Libii, "Applying Dynamics to the Engineering of the
#          Perfect Bounce...", ASEE Annual Conf. 2011, doi:10.18260/1-2--17499
#          Higher inflation -> higher rebound; nice undergrad lab companion.
#
#  Full ASME-format reference list: docs/REFERENCES.md in the repository.
#  Tag -> number there:  NBA [1], FIBA26 [2], FIBA18 [3], NFHS [4], TDT17 [5],
#  SA07 [6], LFR12 [7], Gent [8], HC75 [9], Cross [10]-[12], ASEE [13],
#  Pat1 [14], Pat2 [15], drag Cd [16], SciPy [17], RK45 [18], Brent [19].
#
#  MODEL LIMITATIONS (good discussion topics)
#    * Floor is perfectly rigid. A sprung hardwood floor adds compliance
#      (longer contact, lower peak force).  See "HOW TO EXTEND" at the bottom.
#    * Ball is a lumped mass; shell wave/vibration and the "momentum flux" of
#      shell mass landing on the floor [SA07] are ignored.
#    * Shell stiffness is simply ADDED to the gas stiffness (an approximation;
#      [LFR12] gives the fully coupled pressurised-shell result).
#    * One damping constant lambda lumps together all losses (rubber
#      hysteresis, seams, sound, air motion).  [TDT17] explains why the true
#      dominant loss mechanism is still debated.
#    * Vertical drop only (no spin, no friction).
# =============================================================================


# =============================================================================
# SECTION 1 - IMPORTS
# =============================================================================
import os                     # building file paths for CSV export
import datetime               # time-stamping exported files
import itertools              # making every combination of min/max settings

import numpy as np            # arrays and math
import pandas as pd           # tables -> CSV
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp   # ODE solver for the contact phase
from scipy.optimize import brentq       # 1-D root finder (used to calibrate damping)

__version__ = "1.2.0"   # semantic versioning: MAJOR.MINOR.PATCH (see CHANGELOG.md)


# =============================================================================
# SECTION 2 - PHYSICAL CONSTANTS AND UNIT CONVERSIONS
# Everything inside the physics functions is SI (m, kg, s, N, Pa, J).
# Conversions happen only at the "edges" (sliders, readouts, CSV).
# =============================================================================
G = 9.80665                   # gravitational acceleration [m/s^2]
P_ATM = 101_325.0             # atmospheric pressure [Pa]
GAMMA_AIR = 1.4               # ratio of specific heats for air (adiabatic exponent)
RHO_AIR = 1.20                # air density at ~20 C [kg/m^3]  (used for drag)
CD_BALL = 0.47                # drag coefficient of a sphere (textbook value)
NU_RUBBER = 0.5               # Poisson's ratio of rubber (nearly incompressible)

PSI_TO_PA = 6894.757          # 1 psi in pascals
FT_TO_M = 0.3048              # 1 foot in metres (also 1 ft/s in m/s)
IN_TO_M = 0.0254              # 1 inch in metres
N_TO_LBF = 0.224809           # 1 newton in pounds-force
KPA_PER_PSI = PSI_TO_PA / 1000.0   # 1 psi = 6.894757 kPa
G_PER_OZ = 28.349523          # 1 ounce in grams


# =============================================================================
# SECTION 3 - PUBLISHED PARAMETER RANGES
# Each entry: (minimum, default, maximum, slider step, note/source).
# CHANGE THESE to explore other balls (e.g., size 6, youth, volleyball...).
# =============================================================================
RANGES = {
    # Drop height of the BOTTOM of the ball above the floor [ft].
    # 6 ft = NFHS test, 1.8 m (5.9 ft) = FIBA test, 12 ft = your case.
    # (The GUI can show this slider in metres; values are converted.)
    "drop_height_ft":   (0.5, 12.0, 20.0, 0.5,  "user choice; rule tests use 6 ft / 1.8 m"),

    # Inflation pressure, gauge [psi].  NBA rule: 7.5-8.5 psi.  Slider is a bit
    # wider (5-12) so students can see what an under/over-inflated ball does.
    # (The GUI can show this slider in kPa; values are converted.)
    "pressure_psi":     (5.0, 8.0, 12.0, 0.1,   "NBA Rule 1 II: 7.5-8.5 psi"),

    # Mass [g].  FIBA 2018 size 7: 567-650 g (FIBA 2026 tightened to 580-620 g);
    # NBA ball is 22 oz = 624 g.  Slider uses the wider 2018 range.
    "mass_g":           (567.0, 624.0, 650.0, 1.0, "FIBA size 7: 567-650 g (2018), 580-620 g (2026)"),

    # Circumference [mm].  FIBA 2018 size 7: 749-780 mm (2026: 750-770 mm);
    # NBA 29.5 in = 749 mm.
    "circumference_mm": (749.0, 755.0, 780.0, 1.0, "FIBA size 7: 749-780 mm (2018), 750-770 mm (2026)"),

    # Rubber hardness [Shore A].  Basketball cover rubbers ~45-85 [Pat1].
    "shore_A":          (45.0, 75.0, 85.0, 1.0, "cover rubber 45-85 Shore A [Pat1]"),

    # Total wall (shell) thickness [mm]: bladder + windings + carcass + cover.
    # [TDT17] uses 2.5 mm; [Pat2] layer stack adds to ~4-5 mm.
    "shell_thick_mm":   (2.0, 4.0, 6.0, 0.1, "2.5 mm [TDT17] to ~5 mm [Pat2]"),

    # Starting velocity [m/s], POSITIVE = UPWARD, negative = thrown downward.
    # 0 = simple drop.  +/-15 m/s covers hard dribbles and throws; the impact
    # speed then stays below ~20 m/s, where the spherical-cap model still holds
    # (max squash < 2R/3 [TDT17]).  (GUI can show ft/s; values are converted.)
    "v_init_m_s":       (-15.0, 0.0, 15.0, 0.1, "0 = drop from rest; + up, - down"),
}

# Rule-test definitions used for the "legal ball?" check.
#   drop_m       : release height of the ball's BOTTOM (underside)
#   rebound_m    : allowed (min, max) rebound height
#   measured_to  : "bottom" or "top" - which point of the ball the rule
#                  measures the rebound to (the rules differ!)
# The FIRST entry is the default used to calibrate damping.
RULE_TESTS = {
    "FIBA 2026": {"drop_m": 1.800, "rebound_m": (1.035, 1.085), "measured_to": "bottom"},
    "FIBA 2018": {"drop_m": 1.800, "rebound_m": (1.200, 1.400), "measured_to": "top"},
    "NFHS":      {"drop_m": 6 * FT_TO_M, "rebound_m": (49 * IN_TO_M, 54 * IN_TO_M),
                  "measured_to": "top"},
}
DEFAULT_TEST = "FIBA 2026"

# Heights [ft] used for the min/max sweep plots.
SWEEP_HEIGHTS_FT = np.array([1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float)


def default_params():
    """Return a fresh dictionary of default settings (SI conversion happens later)."""
    p = {key: rng[1] for key, rng in RANGES.items()}
    p["air_drag"] = True          # include air resistance during fall/rebound
    # (p["v_init_m_s"] = 0.0 from RANGES -> a simple drop from rest)
    p["damping_lambda"] = None    # None -> auto-calibrate to the FIBA 2026 mid-band
    return p


# =============================================================================
# SECTION 4 - BALL PHYSICS: GEOMETRY, GAS, SHELL, DAMPING
# =============================================================================
def shore_A_to_youngs_modulus(shore_A):
    """
    Gent (1958) relation between Shore A hardness and Young's modulus.
    Returns E in Pa.  Reasonable for roughly 40 < Shore A < 95.
       60 Shore A -> ~3.6 MPa,  75 -> ~7 MPa,  85 -> ~12 MPa
    """
    S = shore_A
    E_MPa = 0.0981 * (56.0 + 7.62336 * S) / (0.137505 * (254.0 - 2.54 * S))
    return E_MPa * 1e6


def ball_volume(x, R):
    """
    Volume [m^3] of a sphere of radius R whose bottom has been flattened by a
    depth x (a sphere with a cap sliced off - the "indented spherical cap"
    shape seen in high-speed video [TDT17 eq. 2.1]).
        x = 0  ->  (4/3) pi R^3   (undeformed)
    """
    return (np.pi / 3.0) * (4.0 * R**3 + x**3 - 3.0 * R * x**2)


def contact_area(x, R):
    """
    Area [m^2] of the flat circular patch touching the floor [TDT17 eq. 2.1]:
        A = 2 pi R x (1 - x/(2R)) = pi (2 R x - x^2)
    For small x this is ~2 pi R x, i.e. the patch grows linearly with
    deflection - that is why an inflated ball behaves like a linear spring.
    """
    return np.pi * (2.0 * R * x - x**2)


def gas_pressure_abs(x, R, p_gauge0):
    """
    Absolute air pressure [Pa] inside the ball at deflection x.
    Adiabatic (no time for heat transfer during ~10 ms):  P V^gamma = const.
    """
    V0 = ball_volume(0.0, R)
    V = ball_volume(x, R)
    return (p_gauge0 + P_ATM) * (V0 / V) ** GAMMA_AIR


def gas_force(x, R, p_gauge0):
    """Force [N] from the pressurised air acting on the contact patch.
    Only the GAUGE pressure (inside minus outside) produces a net force."""
    return (gas_pressure_abs(x, R, p_gauge0) - P_ATM) * contact_area(x, R)


def gas_energy(x, R, p_gauge0):
    """
    Elastic energy [J] stored in the air at deflection x = work done by the
    floor on the gas = integral of F_gas dx.
      = (internal-energy increase of the air)  -  (work done by the atmosphere
        pushing on the shrinking ball)
      = P0abs V0/(gamma-1) * [ (V0/V)^(gamma-1) - 1 ]  -  P_atm (V0 - V)
    """
    V0 = ball_volume(0.0, R)
    V = ball_volume(x, R)
    P0abs = p_gauge0 + P_ATM
    dU = P0abs * V0 / (GAMMA_AIR - 1.0) * ((V0 / V) ** (GAMMA_AIR - 1.0) - 1.0)
    return dU - P_ATM * (V0 - V)


def shell_stiffness(E, thickness, R, nu=NU_RUBBER):
    """
    Indentation stiffness [N/m] of a thin elastic spherical shell
    (Reissner / [LFR12], as written in [TDT17] eq. 2.13):
        k_shell = 4 E t^2 / ( R * sqrt(3 (1 - nu^2)) )
    """
    return 4.0 * E * thickness**2 / (R * np.sqrt(3.0 * (1.0 - nu**2)))


def floor_force_parts(x, v, ball):
    """
    Compute every piece of the floor force at deflection x [m] and deflection
    rate v = dx/dt [m/s] (v > 0 while compressing, v < 0 while rebounding).

    Returns (F_total, F_gas, F_shell, F_damp)  all in N.
    F_total is clipped at zero: a floor can push but can never pull.
    """
    R = ball["R"]
    F_gas = gas_force(x, R, ball["p_gauge"])
    F_shell = ball["k_shell"] * x
    F_damp = ball["lambda"] * x * v            # Hunt-Crossley damping [HC75]
    F_total = F_gas + F_shell + F_damp
    if np.ndim(F_total) == 0:
        F_total = max(F_total, 0.0)
    else:
        F_total = np.maximum(F_total, 0.0)
    return F_total, F_gas, F_shell, F_damp


def build_ball(params):
    """
    Turn the user-facing settings (psi, grams, mm, Shore A...) into the SI
    quantities the physics functions need.  Returns a dict called 'ball'.
    """
    m = params["mass_g"] / 1000.0
    C = params["circumference_mm"] / 1000.0
    R = C / (2.0 * np.pi)
    p_gauge = params["pressure_psi"] * PSI_TO_PA
    E = shore_A_to_youngs_modulus(params["shore_A"])
    t = params["shell_thick_mm"] / 1000.0
    ball = {
        "m": m, "R": R, "p_gauge": p_gauge, "E": E, "t": t,
        "k_shell": shell_stiffness(E, t, R),
        # small-deflection gas stiffness: dF/dx at x=0 = 2 pi R P_gauge
        "k_gas0": 2.0 * np.pi * R * p_gauge,
        "lambda": params["damping_lambda"] if params["damping_lambda"] is not None else 0.0,
        "air_drag": params["air_drag"],
    }
    return ball


# =============================================================================
# SECTION 5 - FLIGHT PHASES (fall before impact, rise after)
# Quadratic air drag  F = 1/2 rho Cd A v^2  has exact closed-form solutions.
# =============================================================================
def terminal_velocity(ball):
    area = np.pi * ball["R"] ** 2
    return np.sqrt(2.0 * ball["m"] * G / (RHO_AIR * CD_BALL * area))


def rise_height(v_up, ball):
    """
    Extra height [m] gained by a ball moving UP at v_up before it stops.
    No drag:   v^2 / (2g).
    Drag:      vt^2/(2g) * ln(1 + v^2/vt^2)   (exact for quadratic drag)
    """
    if not ball["air_drag"]:
        return v_up**2 / (2.0 * G)
    vt = terminal_velocity(ball)
    return vt**2 / (2.0 * G) * np.log(1.0 + v_up**2 / vt**2)


def impact_speed(h, ball, v_init=0.0):
    """
    Speed [m/s] at the floor for a ball whose bottom starts at height h [m]
    with vertical velocity v_init [m/s] (POSITIVE = UPWARD).

    Method:
      * If v_init > 0 the ball first rises by rise_height(v_init) and stops,
        so it then falls FROM REST through h + rise.
      * If v_init <= 0 it falls through h starting with downward speed
        |v_init|.  For quadratic drag, d(v^2)/dy = 2g(1 - v^2/vt^2), which
        integrates exactly to
            v^2 = vt^2 + (v_start^2 - vt^2) * exp(-2 g y / vt^2)
        (a ball thrown down faster than vt slows toward vt).
      * Without drag both cases reduce to  v0 = sqrt(v_init^2 + 2 g h).
    """
    if v_init > 0.0:
        h_total, v_start = h + rise_height(v_init, ball), 0.0
    else:
        h_total, v_start = h, -v_init
    if not ball["air_drag"]:
        return np.sqrt(v_start**2 + 2.0 * G * h_total)
    vt = terminal_velocity(ball)
    v2 = vt**2 + (v_start**2 - vt**2) * np.exp(-2.0 * G * h_total / vt**2)
    return np.sqrt(v2)


def rebound_height(v_up, ball):
    """Height [m] the ball's bottom rises to after leaving the floor at v_up."""
    return rise_height(v_up, ball)


# =============================================================================
# SECTION 6 - THE CONTACT SIMULATION
# =============================================================================
def _contact_ode(t, y, ball):
    """
    Right-hand side of the ODE system, in the form solve_ivp wants.
    State vector y = [x, v, E_diss]
        x      deflection (centre displacement since touch-down)   [m]
        v      dx/dt                                                [m/s]
        E_diss energy dissipated so far (running integral)          [J]
    """
    x, v, _ = y
    F_total, F_gas, F_shell, _ = floor_force_parts(x, v, ball)
    a = G - F_total / ball["m"]                 # Newton's 2nd law (down = +)
    # Power lost = (actual force - purely elastic force) * velocity.
    # This is the damping power, correctly reduced if the clip at F=0 is active.
    P_diss = (F_total - F_gas - F_shell) * v
    return [v, a, P_diss]


def _liftoff_event(t, y, ball):
    """solve_ivp stops when this crosses zero going downward: x back to 0."""
    return y[0]
_liftoff_event.terminal = True
_liftoff_event.direction = -1


def simulate_contact(v0, ball, n_points=600):
    """
    Integrate one contact from touch-down (x=0, v=v0) to lift-off (x=0 again).
    Returns a dict with a pandas DataFrame time series and summary numbers.
    """
    m = ball["m"]
    # Rough contact-time estimate to set the integration window (x3 safety).
    k_est = ball["k_gas0"] + ball["k_shell"]
    t_window = 3.0 * np.pi * np.sqrt(m / k_est) + 0.01

    sol = solve_ivp(_contact_ode, (0.0, t_window), [0.0, v0, 0.0],
                    args=(ball,), events=_liftoff_event,
                    rtol=1e-9, atol=[1e-12, 1e-10, 1e-12],
                    dense_output=True)
    if sol.t_events[0].size == 0:
        raise RuntimeError("Ball never left the floor - damping too large?")

    t_c = sol.t_events[0][0]                    # contact duration [s]
    v_out = -sol.y_events[0][0][1]              # rebound speed (positive up)

    # Sample the smooth (dense) solution on an even time grid for plots/CSV.
    t = np.linspace(0.0, t_c, n_points)
    x, v, E_diss = sol.sol(t)
    x = np.clip(x, 0.0, None)                   # remove -1e-15 round-off at ends
    F_tot, F_gas, F_sh, F_d = floor_force_parts(x, v, ball)

    KE = 0.5 * m * v**2
    E_gas = gas_energy(x, ball["R"], ball["p_gauge"])
    E_sh = 0.5 * ball["k_shell"] * x**2
    W_grav = m * G * x                          # energy gravity adds during contact
    KE0 = 0.5 * m * v0**2
    balance_err = (KE + E_gas + E_sh + E_diss) - (KE0 + W_grav)

    P_in_gauge = gas_pressure_abs(x, ball["R"], ball["p_gauge"]) - P_ATM
    A = contact_area(x, ball["R"])

    df = pd.DataFrame({
        "time_ms": t * 1e3,
        "deflection_mm": x * 1e3,
        "velocity_m_s": v,
        "contact_area_cm2": A * 1e4,
        "contact_patch_diam_mm": 2.0 * np.sqrt(A / np.pi) * 1e3,
        "pressure_gauge_psi": P_in_gauge / PSI_TO_PA,
        "F_gas_N": F_gas,
        "F_shell_N": F_sh,
        "F_damping_N": F_d,
        "F_floor_N": F_tot,
        "F_floor_lbf": F_tot * N_TO_LBF,
        "KE_J": KE,
        "E_gas_stored_J": E_gas,
        "E_shell_stored_J": E_sh,
        "E_dissipated_J": E_diss,
        "W_gravity_J": W_grav,
        "energy_balance_error_J": balance_err,
    })

    i_pk = int(np.argmax(F_tot))
    i_mx = int(np.argmax(x))
    return {
        "timeseries": df,
        "v0": v0, "v_out": v_out, "KE0": KE0,
        "COR": v_out / v0,                       # coefficient of restitution
        "contact_time_ms": t_c * 1e3,
        "peak_force_N": F_tot[i_pk],
        "time_of_peak_ms": t[i_pk] * 1e3,
        "max_deflection_mm": x[i_mx] * 1e3,
        "i_maxdefl": i_mx,
        "max_area_cm2": A[i_mx] * 1e4,
        "peak_pressure_psi": P_in_gauge[i_mx] / PSI_TO_PA,
        "max_energy_error_J": np.max(np.abs(balance_err)),
    }


def simulate_drop(params):
    """
    Full experiment: start at params['drop_height_ft'] with velocity
    params['v_init_m_s'] (+ up; 0 = drop from rest), contact, rebound.
    Adds the flight results to the dict returned by simulate_contact().
    """
    params = dict(params)
    if params.get("damping_lambda") is None:
        params["damping_lambda"] = calibrate_damping(params)
    ball = build_ball(params)
    h = params["drop_height_ft"] * FT_TO_M
    v_init = params.get("v_init_m_s", 0.0)
    v0 = impact_speed(h, ball, v_init)
    if v0 < 1e-6:
        raise ValueError("The ball reaches the floor with zero speed - "
                         "raise the height or give it a starting velocity.")
    res = simulate_contact(v0, ball)
    res["ball"] = ball
    res["params"] = params
    res["drop_height_m"] = h
    res["v_init_m_s"] = v_init
    # highest point before impact (only above h if thrown upward)
    res["apex_height_m"] = h + (rise_height(v_init, ball) if v_init > 0 else 0.0)
    # drop height (from rest, same drag setting) that would give the same impact
    # speed - a handy way to compare a throw with a simple drop
    res["equiv_drop_height_m"] = (brentq(lambda hh: impact_speed(hh, ball) - v0, 0.0, 1e4)
                                  if (not ball["air_drag"] or v0 < 0.999 * terminal_velocity(ball))
                                  else np.inf)
    res["rebound_height_m"] = rebound_height(res["v_out"], ball)
    # simple linear-spring theory for comparison [TDT17 eq. 2.12]:
    #   t_c = pi sqrt(m / k),  k = 2 pi R P_gauge (+ shell)
    res["t_linear_theory_ms"] = np.pi * np.sqrt(ball["m"] / (ball["k_gas0"] + ball["k_shell"])) * 1e3
    return res


# =============================================================================
# SECTION 7 - RULE-TEST CHECK AND DAMPING CALIBRATION
# The rules don't give a damping value; they give a REBOUND TEST.  So we
# choose lambda so that the simulated ball passes that test.
# =============================================================================
def rule_test_rebound(params, lam, test=DEFAULT_TEST):
    """
    Simulated rebound height [m] in an official drop test, measured the way
    that rule measures it (to the underside OR to the top of the ball).
    """
    p = dict(params)
    p["damping_lambda"] = lam
    ball = build_ball(p)
    spec = RULE_TESTS[test]
    v0 = impact_speed(spec["drop_m"], ball)
    res = simulate_contact(v0, ball, n_points=50)
    h_bottom = rebound_height(res["v_out"], ball)
    return h_bottom + (2.0 * ball["R"] if spec["measured_to"] == "top" else 0.0)


def calibrate_damping(params, target_fraction=0.5, test=DEFAULT_TEST):
    """
    Find lambda so the ball rebounds to a chosen point inside the legal band
    of the rule test.  target_fraction = 0 -> lowest legal rebound (most
    damping), 1 -> highest legal rebound (least damping), 0.5 -> middle.
    Uses Brent's root-finding method on lambda.
    """
    lo, hi = RULE_TESTS[test]["rebound_m"]
    target = lo + target_fraction * (hi - lo)
    f = lambda lam: rule_test_rebound(params, lam, test) - target
    return brentq(f, 0.0, 50_000.0, xtol=1.0)


def rule_check(params, lam):
    """Return {test_name: (rebound_m, is_legal)} for every rule test."""
    out = {}
    for name, spec in RULE_TESTS.items():
        h = rule_test_rebound(params, lam, name)
        lo, hi = spec["rebound_m"]
        out[name] = (h, lo <= h <= hi)
    return out


# =============================================================================
# SECTION 8 - MIN / MAX RANGE SWEEP
# "What is the shortest and longest contact time (and lowest/highest peak
#  force) for any regulation-plausible ball?"  We run every corner of the
# parameter box (each parameter at its min or max) for several drop heights
# and take the min and max of each output.  Response is monotonic in each
# parameter, so the corners bracket the full range.
# =============================================================================
SWEEP_BOX = {
    "pressure_psi":   (7.5, 8.5),     # NBA legal range
    "shore_A":        (45.0, 85.0),   # [Pat1]
    "shell_thick_mm": (2.5, 5.0),     # [TDT17] .. [Pat2]
    # damping: from "deadest legal" to "liveliest legal" ball in the
    # DEFAULT_TEST (FIBA 2026) rebound test
    "damping":        ("rule_low", "rule_high"),
}


def sweep_ranges(params, heights_ft=SWEEP_HEIGHTS_FT):
    """
    Returns a DataFrame with, for each drop height, the min and max contact
    time and peak force over all corner combinations of SWEEP_BOX.
    Mass, size, air-drag setting and starting velocity are taken from params.
    """
    base = dict(params)
    v_init = base.get("v_init_m_s", 0.0)
    # Damping values that put an 8.0 psi mid-hardness ball at the bottom /
    # top of the DEFAULT_TEST legal rebound band.
    ref = dict(base, pressure_psi=8.0, shore_A=65.0, shell_thick_mm=3.75)
    lam_dead = calibrate_damping(ref, target_fraction=0.0)
    lam_live = calibrate_damping(ref, target_fraction=1.0)

    rows = []
    corners = itertools.product(SWEEP_BOX["pressure_psi"], SWEEP_BOX["shore_A"],
                                SWEEP_BOX["shell_thick_mm"], (lam_live, lam_dead))
    corner_balls = []
    for P, S, th, lam in corners:
        p = dict(base, pressure_psi=P, shore_A=S, shell_thick_mm=th, damping_lambda=lam)
        corner_balls.append(build_ball(p))

    for h_ft in heights_ft:
        tc, fp, cor = [], [], []
        for ball in corner_balls:
            r = simulate_contact(impact_speed(h_ft * FT_TO_M, ball, v_init), ball, n_points=200)
            tc.append(r["contact_time_ms"])
            fp.append(r["peak_force_N"])
            cor.append(r["COR"])
        rows.append({"drop_height_ft": h_ft,
                     "contact_time_min_ms": min(tc), "contact_time_max_ms": max(tc),
                     "peak_force_min_N": min(fp), "peak_force_max_N": max(fp),
                     "COR_min": min(cor), "COR_max": max(cor)})
    df = pd.DataFrame(rows)
    df.attrs["lambda_range"] = (lam_live, lam_dead)
    return df


def sweep_current_ball(params, heights_ft=SWEEP_HEIGHTS_FT):
    """Contact time & peak force vs height for the ball currently on the sliders."""
    ball = build_ball(params)
    v_init = params.get("v_init_m_s", 0.0)
    rows = []
    for h_ft in heights_ft:
        r = simulate_contact(impact_speed(h_ft * FT_TO_M, ball, v_init), ball, n_points=200)
        rows.append({"drop_height_ft": h_ft, "contact_time_ms": r["contact_time_ms"],
                     "peak_force_N": r["peak_force_N"], "COR": r["COR"]})
    return pd.DataFrame(rows)


# =============================================================================
# SECTION 9 - PLOTTING
# =============================================================================
COLORS = {"floor": "#1f4e79", "gas": "#2e86c1", "shell": "#7dcea0",
          "damp": "#e67e22", "ke": "#8e44ad", "grav": "#7f8c8d", "band": "#aed6f1"}


def plot_bounce(res):
    """Four-panel figure for one bounce. Returns the matplotlib Figure."""
    h_ft = res["params"]["drop_height_ft"]
    df = res["timeseries"]
    t = df["time_ms"]
    i_mx = res["i_maxdefl"]
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
    v_i = res.get("v_init_m_s", 0.0)
    if v_i == 0.0:
        start = f"Drop {h_ft:.1f} ft ({h_ft*FT_TO_M:.2f} m)"
    else:
        start = (f"Start {h_ft:.1f} ft ({h_ft*FT_TO_M:.2f} m), "
                 f"v_init {v_i:+.1f} m/s ({v_i/FT_TO_M:+.1f} ft/s) {'up' if v_i > 0 else 'down'}")
    fig.suptitle(
        f"{start}  |  "
        f"{res['params']['pressure_psi']:.2f} psi ({res['params']['pressure_psi']*KPA_PER_PSI:.1f} kPa)  |  "
        f"impact {res['v0']:.2f} m/s\n"
        f"contact {res['contact_time_ms']:.2f} ms  |  peak floor force "
        f"{res['peak_force_N']:.0f} N ({res['peak_force_N']*N_TO_LBF:.0f} lbf)",
        fontsize=12, fontweight="bold")

    # (a) Floor force vs time, broken into its parts --------------------------
    ax = axs[0, 0]
    ax.plot(t, df["F_floor_N"], color=COLORS["floor"], lw=2.5, label="Total floor force")
    ax.plot(t, df["F_gas_N"], color=COLORS["gas"], lw=1.2, ls="--", label="Gas (air pressure)")
    ax.plot(t, df["F_shell_N"], color=COLORS["shell"], lw=1.2, ls="--", label="Rubber shell")
    ax.plot(t, df["F_damping_N"], color=COLORS["damp"], lw=1.2, ls=":", label="Damping (loss)")
    ax.axhline(0, color="k", lw=0.6)
    ax.plot(res["time_of_peak_ms"], res["peak_force_N"], "o", color="crimson")
    ax.annotate(f"peak {res['peak_force_N']:.0f} N", (res["time_of_peak_ms"], res["peak_force_N"]),
                xytext=(8, -4), textcoords="offset points", color="crimson")
    ax.set_xlabel("time since touch-down [ms]")
    ax.set_ylabel("force [N]")
    ax.set_title("(a) Force on the floor during contact")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

    # (b) Energy vs time -------------------------------------------------------
    ax = axs[0, 1]
    ax.stackplot(t, df["KE_J"], df["E_gas_stored_J"], df["E_shell_stored_J"], df["E_dissipated_J"],
                 colors=[COLORS["ke"], COLORS["gas"], COLORS["shell"], COLORS["damp"]],
                 labels=["Kinetic", "Stored in gas (PE)", "Stored in shell (PE)", "Dissipated (heat)"],
                 alpha=0.85)
    ax.plot(t, res["KE0"] + df["W_gravity_J"], "k--", lw=1,
            label="Impact KE + gravity work (should = top of stack)")
    ax.axvline(t.iloc[i_mx], color="k", lw=0.8, ls=":")
    ax.text(t.iloc[i_mx], 0.02 * res["KE0"], " max deflection", fontsize=8, rotation=90, va="bottom")
    ax.set_xlabel("time since touch-down [ms]")
    ax.set_ylabel("energy [J]")
    ax.set_title("(b) Where the kinetic energy goes")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)

    # (c) Force vs deflection (hysteresis loop) --------------------------------
    ax = axs[1, 0]
    comp = df["velocity_m_s"] >= 0
    ax.plot(df["deflection_mm"][comp], df["F_floor_N"][comp], color="firebrick", lw=2,
            label="compression (going down)")
    ax.plot(df["deflection_mm"][~comp], df["F_floor_N"][~comp], color="seagreen", lw=2,
            label="restitution (coming back)")
    ax.plot(df["deflection_mm"], df["F_gas_N"] + df["F_shell_N"], "k--", lw=1,
            label="elastic only (no loss)")
    ax.fill(df["deflection_mm"], df["F_floor_N"], color="orange", alpha=0.25,
            label=f"loop area = energy lost = {df['E_dissipated_J'].iloc[-1]:.2f} J")
    ax.set_xlabel("deflection x [mm]")
    ax.set_ylabel("floor force [N]")
    ax.set_title("(c) Force-deflection loop")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    # (d) Energy budget bars at max deflection and at lift-off ------------------
    ax = axs[1, 1]
    rows = {"At max\ndeflection": i_mx, "At lift-off": len(df) - 1}
    keys = [("KE_J", "Kinetic", COLORS["ke"]),
            ("E_gas_stored_J", "Gas PE", COLORS["gas"]),
            ("E_shell_stored_J", "Shell PE", COLORS["shell"]),
            ("E_dissipated_J", "Dissipated", COLORS["damp"])]
    total_in = res["KE0"] + df["W_gravity_J"]
    for j, (lab, i) in enumerate(rows.items()):
        left = 0.0
        for col, name, c in keys:
            val = df[col].iloc[i]
            pct = 100.0 * val / total_in.iloc[i]
            ax.barh(j, val, left=left, color=c, edgecolor="white", label=name if j == 0 else None)
            if pct > 4:
                ax.text(left + val / 2, j, f"{pct:.0f}%", ha="center", va="center",
                        color="white", fontsize=9, fontweight="bold")
            left += val
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(list(rows.keys()))
    ax.invert_yaxis()
    ax.set_xlabel("energy [J]  (impact KE = %.2f J)" % res["KE0"])
    ax.set_title("(d) Energy budget (KE -> potential energy -> KE + heat)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=4, frameon=False)
    ax.grid(alpha=0.3, axis="x")

    fig.tight_layout()
    return fig


def plot_sweep(sweep_df, current_df, current_height_ft, height_unit="ft", v_init=0.0):
    """Min/max envelope of contact time and peak force vs drop height.
    height_unit = "ft" or "m" sets the x-axis units.  v_init (m/s, + up) is
    only used for the figure title - the sweep data already includes it."""
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
    if v_init != 0.0:
        fig.suptitle(f"Every point starts with v_init = {v_init:+.1f} m/s "
                     f"({v_init/FT_TO_M:+.1f} ft/s, + = up)", fontsize=10)
    k = 1.0 if height_unit == "ft" else FT_TO_M      # ft -> display unit factor
    h = sweep_df["drop_height_ft"] * k
    for ax, lo, hi, cur, ylab, title in [
        (axs[0], "contact_time_min_ms", "contact_time_max_ms", "contact_time_ms",
         "contact duration [ms]", "(e) Contact duration: min-max range"),
        (axs[1], "peak_force_min_N", "peak_force_max_N", "peak_force_N",
         "peak floor force [N]", "(f) Peak floor force: min-max range")]:
        ax.fill_between(h, sweep_df[lo], sweep_df[hi], color=COLORS["band"], alpha=0.8,
                        label="regulation-range envelope")
        ax.plot(h, sweep_df[lo], color="steelblue", lw=1)
        ax.plot(h, sweep_df[hi], color="steelblue", lw=1)
        ax.plot(current_df["drop_height_ft"] * k, current_df[cur], "o-", color="crimson",
                ms=4, lw=1.5, label="current slider ball")
        ax.axvline(current_height_ft * k, color="k", ls=":", lw=1)
        ax.set_xlabel(f"{'drop' if v_init == 0 else 'start'} height (bottom of ball) [{height_unit}]")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return fig


# =============================================================================
# SECTION 10 - CSV EXPORT
# =============================================================================
def export_csv(res, filename, sweep_df=None):
    """
    Write the bounce time series to CSV.  Lines starting with '#' at the top
    record every setting and summary result so the file is self-describing.
    Re-load later with:  pd.read_csv(filename, comment="#")
    Returns the list of files written.
    """
    if not filename.lower().endswith(".csv"):
        filename += ".csv"
    path = os.path.abspath(filename)
    b = res["ball"]
    meta = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "tool_version": __version__,
        "drop_height_ft": res["params"]["drop_height_ft"],
        "drop_height_m": res["drop_height_m"],
        "v_init_m_s_positive_up": res.get("v_init_m_s", 0.0),
        "v_init_ft_s_positive_up": res.get("v_init_m_s", 0.0) / FT_TO_M,
        "apex_height_m": res.get("apex_height_m", res["drop_height_m"]),
        "equivalent_drop_height_m": res.get("equiv_drop_height_m", res["drop_height_m"]),
        "pressure_gauge_psi": res["params"]["pressure_psi"],
        "pressure_gauge_kPa": res["params"]["pressure_psi"] * KPA_PER_PSI,
        "mass_g": res["params"]["mass_g"],
        "circumference_mm": res["params"]["circumference_mm"],
        "radius_m": b["R"],
        "shore_A": res["params"]["shore_A"],
        "youngs_modulus_MPa": b["E"] / 1e6,
        "shell_thickness_mm": res["params"]["shell_thick_mm"],
        "k_gas_small_defl_N_per_m": b["k_gas0"],
        "k_shell_N_per_m": b["k_shell"],
        "damping_lambda_Ns_per_m2": b["lambda"],
        "damping_calibration_rule": DEFAULT_TEST,
        "air_drag": b["air_drag"],
        "impact_speed_m_s": res["v0"],
        "rebound_speed_m_s": res["v_out"],
        "COR": res["COR"],
        "rebound_height_m": res["rebound_height_m"],
        "contact_time_ms": res["contact_time_ms"],
        "peak_force_N": res["peak_force_N"],
        "max_deflection_mm": res["max_deflection_mm"],
        "max_energy_balance_error_J": res["max_energy_error_J"],
    }
    # newline="" stops Windows from inserting blank lines between rows
    with open(path, "w", newline="", encoding="utf-8") as f:
        for k, v in meta.items():
            f.write(f"# {k}, {v}\n")
        res["timeseries"].to_csv(f, index=False)
    written = [path]
    if sweep_df is not None:
        sweep_path = path[:-4] + "_range_sweep.csv"
        sw = sweep_df.copy()
        sw.insert(1, "drop_height_m", sw["drop_height_ft"] * FT_TO_M)
        sw.to_csv(sweep_path, index=False)
        written.append(sweep_path)
    return written


# =============================================================================
# SECTION 11 - INTERACTIVE PANEL (ipywidgets)
#
# Layout ideas (so it still works in a narrow browser window):
#   * Each slider is a small vertical "block":   label  /  slider  /  value.
#     The current value is printed UNDER the slider (not to its right), in
#     both unit systems, so it never gets cut off.
#   * Slider widths are 100% of their column, and the columns sit in a
#     flex box that WRAPS: wide window -> side by side, narrow -> stacked.
#
# Units: the physics always receives drop height in ft, pressure in psi and
# starting velocity in m/s (the keys of RANGES).  The GUI may DISPLAY metres,
# kPa or ft/s; height_ft(), pressure_psi() and v_init_m_s() below convert
# back before simulating.
#
# Initial condition: a toggle chooses "Drop from rest" (v_init = 0, velocity
# slider hidden) or "Initial velocity" (slider shown; + = upward).
#
# Sliders call update() only when released (continuous_update=False) so the
# plots don't try to redraw dozens of times while you drag.
# =============================================================================

# Slider settings for each display unit: (min, max, step, number format)
UNIT_SLIDER = {
    "ft":  (0.5, 20.0, 0.5, ".1f"),
    "m":   (0.15, 6.10, 0.01, ".2f"),
    "psi": (5.0, 12.0, 0.05, ".2f"),
    "kPa": (34.5, 82.7, 0.1, ".1f"),
    "m/s": (-15.0, 15.0, 0.1, ".1f"),
    "ft/s": (-49.0, 49.0, 0.25, ".1f"),
}
# Conversion factor FROM the canonical unit (ft, psi, m/s) TO each display unit
UNIT_FACTOR = {"ft": 1.0, "m": FT_TO_M, "psi": 1.0, "kPa": KPA_PER_PSI,
               "m/s": 1.0, "ft/s": 1.0 / FT_TO_M}


def launch(height_unit="ft", pressure_unit="psi", velocity_unit="m/s"):
    """Build and display the interactive panel.  Returns the top widget."""
    import ipywidgets as w
    from IPython.display import display, FileLink

    # ---------------------------------------------------------------- helpers
    def slider(lo, hi, step, value):
        """A slider with its built-in readout turned OFF (we print our own)."""
        return w.FloatSlider(value=value, min=lo, max=hi, step=step, readout=False,
                             continuous_update=False,
                             layout=w.Layout(width="100%", min_width="160px"))

    def block(title, sl, value_html, extra=None):
        """Stack: bold title, (optional unit toggle), slider, value text."""
        top = [w.HTML(f"<b>{title}</b>")]
        if extra is not None:
            top.append(extra)
        head = w.HBox(top, layout=w.Layout(justify_content="space-between",
                                           align_items="center", width="100%"))
        return w.VBox([head, sl, value_html],
                      layout=w.Layout(width="100%", margin="0 0 6px 0"))

    def column(children):
        """A column that grows/shrinks and wraps under its neighbour when narrow."""
        return w.VBox(children, layout=w.Layout(flex="1 1 320px", min_width="260px",
                                                max_width="520px", padding="0 12px 0 0"))

    def small_html():
        return w.HTML(layout=w.Layout(margin="-4px 0 0 0"))

    # ---------------------------------------------------------- unit toggles
    tg_h = w.ToggleButtons(options=["ft", "m"], value=height_unit,
                           style={"button_width": "42px"}, tooltips=["feet", "metres"])
    tg_p = w.ToggleButtons(options=["psi", "kPa"], value=pressure_unit,
                           style={"button_width": "46px"}, tooltips=["psi (gauge)", "kPa (gauge)"])
    tg_v = w.ToggleButtons(options=["m/s", "ft/s"], value=velocity_unit,
                           style={"button_width": "46px"}, tooltips=["metres per second",
                                                                     "feet per second"])
    # initial-condition mode: simple drop, or start with a velocity
    tg_mode = w.ToggleButtons(options=["Drop from rest", "Initial velocity"],
                              value="Drop from rest", style={"button_width": "auto"},
                              tooltips=["release from rest at the start height",
                                        "start with a vertical velocity (+ up, - down)"])

    # ---------------------------------------------------------------- sliders
    def r(key):                          # shorthand for the RANGES tuple
        return RANGES[key]

    lo, hi, st, _ = UNIT_SLIDER[height_unit]
    s_h = slider(lo, hi, st, r("drop_height_ft")[1] * UNIT_FACTOR[height_unit])
    lo, hi, st, _ = UNIT_SLIDER[pressure_unit]
    s_p = slider(lo, hi, st, round(r("pressure_psi")[1] * UNIT_FACTOR[pressure_unit] / st) * st)
    s_m = slider(r("mass_g")[0], r("mass_g")[2], r("mass_g")[3], r("mass_g")[1])
    s_c = slider(r("circumference_mm")[0], r("circumference_mm")[2],
                 r("circumference_mm")[3], r("circumference_mm")[1])
    s_S = slider(r("shore_A")[0], r("shore_A")[2], r("shore_A")[3], r("shore_A")[1])
    s_t = slider(r("shell_thick_mm")[0], r("shell_thick_mm")[2],
                 r("shell_thick_mm")[3], r("shell_thick_mm")[1])
    s_lam = slider(0.0, 6000.0, 25.0, 1000.0)
    lo, hi, st, _ = UNIT_SLIDER[velocity_unit]
    s_v = slider(lo, hi, st, 0.0)

    v_h, v_p, v_m, v_c, v_S, v_t, v_lam, v_v = (small_html() for _ in range(8))
    title_h = w.HTML("<b>Drop height (bottom of ball)</b>")   # changes with mode

    cb_drag = w.Checkbox(value=True, description="Include air drag", indent=False)
    btn_cal = w.Button(description=f"Calibrate λ to {DEFAULT_TEST}",
                       tooltip="Set damping so the ball rebounds to the middle of the legal band",
                       layout=w.Layout(width="auto"))
    btn_reset = w.Button(description="Reset defaults", layout=w.Layout(width="auto"))
    txt_file = w.Text(value="basketball_bounce", description="CSV name:",
                      layout=w.Layout(width="100%", max_width="320px"))
    cb_sweep_csv = w.Checkbox(value=True, description="also save range sweep", indent=False)
    btn_csv = w.Button(description="Export CSV", button_style="success",
                       layout=w.Layout(width="auto"))

    out_info = w.HTML()
    out_plot = w.Output()
    out_sweep = w.Output()
    out_csv = w.Output()

    state = {"res": None, "sweep": None, "sweep_key": None, "busy": False}

    # ------------------------------------------------- unit-aware conversions
    def height_ft():
        return s_h.value / UNIT_FACTOR[tg_h.value]

    def pressure_psi():
        return s_p.value / UNIT_FACTOR[tg_p.value]

    def v_init_m_s():
        """Starting velocity in m/s (+ up). Zero in 'Drop from rest' mode."""
        if tg_mode.value == "Drop from rest":
            return 0.0
        return s_v.value / UNIT_FACTOR[tg_v.value]

    def current_params():
        return {"drop_height_ft": height_ft(), "pressure_psi": pressure_psi(),
                "mass_g": s_m.value, "circumference_mm": s_c.value,
                "shore_A": s_S.value, "shell_thick_mm": s_t.value,
                "damping_lambda": s_lam.value, "air_drag": cb_drag.value,
                "v_init_m_s": v_init_m_s()}

    def fmt_h(h_ft):
        """Height in the chosen unit first, the other unit in brackets."""
        if tg_h.value == "ft":
            return f"{h_ft:.2f} ft ({h_ft*FT_TO_M:.2f} m)"
        return f"{h_ft*FT_TO_M:.2f} m ({h_ft:.2f} ft)"

    def fmt_v(v_ms):
        """Velocity in the chosen unit first, signed, with a direction word."""
        word = "upward" if v_ms > 0 else ("downward" if v_ms < 0 else "at rest")
        if tg_v.value == "m/s":
            return f"{v_ms:+.1f} m/s ({v_ms/FT_TO_M:+.1f} ft/s) {word}"
        return f"{v_ms/FT_TO_M:+.1f} ft/s ({v_ms:+.1f} m/s) {word}"

    def fmt_p(p_psi):
        if tg_p.value == "psi":
            return f"{p_psi:.2f} psi ({p_psi*KPA_PER_PSI:.1f} kPa)"
        return f"{p_psi*KPA_PER_PSI:.1f} kPa ({p_psi:.2f} psi)"

    # ---------------------------------------------- value text under sliders
    def refresh_value_labels():
        span = "<span style='font-family:monospace;color:#1f4e79'>{}</span>"
        note = "<span style='color:#666;font-size:12px'> {}</span>"
        p_psi = pressure_psi()
        legal = "NBA-legal" if 7.5 <= p_psi + 1e-9 and p_psi - 1e-9 <= 8.5 else "outside NBA 7.5–8.5 psi"
        v_h.value = span.format(fmt_h(height_ft()))
        v_p.value = span.format(fmt_p(p_psi)) + note.format(legal)
        v_m.value = span.format(f"{s_m.value:.0f} g ({s_m.value/G_PER_OZ:.1f} oz)")
        v_c.value = span.format(f"{s_c.value:.0f} mm ({s_c.value/25.4:.2f} in)") + \
            note.format(f"diameter {s_c.value/np.pi:.0f} mm")
        v_S.value = span.format(f"{s_S.value:.0f} Shore A") + \
            note.format(f"→ E ≈ {shore_A_to_youngs_modulus(s_S.value)/1e6:.1f} MPa (Gent)")
        v_t.value = span.format(f"{s_t.value:.1f} mm ({s_t.value/25.4:.3f} in)")
        v_lam.value = span.format(f"{s_lam.value:.0f} N·s/m²") + \
            note.format("0 = perfectly elastic")
        v_v.value = span.format(fmt_v(s_v.value / UNIT_FACTOR[tg_v.value])) + \
            note.format("+ = up, − = down")

    def info_html(res, rules, sweep_df):
        b = res["ball"]
        i = res["i_maxdefl"]
        df = res["timeseries"]
        E0 = res["KE0"] + df["W_gravity_J"].iloc[i]
        pct = lambda col: 100 * df[col].iloc[i] / E0
        rule_parts = []
        for n, (h, ok) in rules.items():
            where = RULE_TESTS[n]["measured_to"]
            lo, hi = RULE_TESTS[n]["rebound_m"]
            verdict = ("<span style='color:green'>LEGAL</span>" if ok
                       else "<span style='color:red'>not legal</span>")
            rule_parts.append(f"<b>{n}</b>: {h*1000:.0f} mm to {where} "
                              f"(legal {lo*1000:.0f}–{hi*1000:.0f}) → {verdict}")
        row = sweep_df.iloc[(sweep_df["drop_height_ft"] - height_ft()).abs().argmin()]
        weight_N = b["m"] * G
        if res["v_init_m_s"] == 0.0:
            start_line = ""
        else:
            eq = res["equiv_drop_height_m"]
            eq_txt = (fmt_h(eq / FT_TO_M) if np.isfinite(eq)
                      else "none (faster than a drop from any height can reach)")
            start_line = (f"<b>Start</b> {fmt_h(height_ft())} with {fmt_v(res['v_init_m_s'])}"
                          + (f" · <b>apex</b> {fmt_h(res['apex_height_m']/FT_TO_M)}"
                             if res["v_init_m_s"] > 0 else "")
                          + f" · <b>equivalent drop height</b> {eq_txt}<br>")
        return f"""
        <div style="font-family:monospace;font-size:13px;line-height:1.55;
                    overflow-wrap:anywhere;max-width:100%">
        {start_line}<b>Impact speed</b> {res['v0']:.2f} m/s · <b>rebound speed</b> {res['v_out']:.2f} m/s
        · <b>COR</b> {res['COR']:.3f} · <b>rebound height</b> {fmt_h(res['rebound_height_m']/FT_TO_M)}<br>
        <b>Contact time</b> {res['contact_time_ms']:.2f} ms
        (linear-spring theory π√(m/k): {res['t_linear_theory_ms']:.2f} ms) ·
        <b>Peak floor force</b> {res['peak_force_N']:.0f} N = {res['peak_force_N']*N_TO_LBF:.0f} lbf
        = {res['peak_force_N']/weight_N:.0f}× ball weight<br>
        <b>Max deflection</b> {res['max_deflection_mm']:.1f} mm ·
        <b>max contact patch</b> {res['max_area_cm2']:.0f} cm² ·
        <b>peak internal pressure</b> {fmt_p(res['peak_pressure_psi'])} gauge<br>
        <b>Stiffness</b>: gas {b['k_gas0']/1000:.1f} kN/m, shell {b['k_shell']/1000:.2f} kN/m
        (shell = {100*b['k_shell']/(b['k_gas0']+b['k_shell']):.1f}% of total)<br>
        <b>Energy at max deflection</b>: gas {pct('E_gas_stored_J'):.1f}% · shell {pct('E_shell_stored_J'):.1f}%
        · already dissipated {pct('E_dissipated_J'):.1f}% · KE {pct('KE_J'):.1f}%
        (energy-balance error {res['max_energy_error_J']:.1e} J)<br>
        <b>Rule checks</b>: {' · '.join(rule_parts)}<br>
        <b>Regulation-range envelope at {fmt_h(row['drop_height_ft'])}</b>:
        contact {row['contact_time_min_ms']:.2f}–{row['contact_time_max_ms']:.2f} ms,
        peak force {row['peak_force_min_N']:.0f}–{row['peak_force_max_N']:.0f} N,
        COR {row['COR_min']:.2f}–{row['COR_max']:.2f}
        </div>"""

    # ------------------------------------------------------------ main update
    def update(*_):
        if state["busy"]:
            return
        refresh_value_labels()
        p = current_params()
        try:
            # The min/max sweep only depends on mass, size, drag and the
            # starting velocity -> cache it on those.
            key = (p["mass_g"], p["circumference_mm"], p["air_drag"], p["v_init_m_s"])
            if state["sweep_key"] != key:
                out_info.value = "<i>Computing regulation-range sweep (a few seconds)…</i>"
                state["sweep"] = sweep_ranges(p)
                state["sweep_key"] = key
            res = simulate_drop(p)
            rules = rule_check(p, p["damping_lambda"])
            cur = sweep_current_ball(p)
        except Exception as e:                       # show errors in the panel
            out_info.value = f"<b style='color:red'>Error:</b> {e}"
            return
        state["res"] = res
        out_info.value = info_html(res, rules, state["sweep"])
        with out_plot:
            out_plot.clear_output(wait=True)
            fig = plot_bounce(res)
            display(fig)
            plt.close(fig)
        with out_sweep:
            out_sweep.clear_output(wait=True)
            fig = plot_sweep(state["sweep"], cur, height_ft(), height_unit=tg_h.value,
                             v_init=p["v_init_m_s"])
            display(fig)
            plt.close(fig)

    # --------------------------------------------------- unit toggle handling
    def set_slider_range(sl, unit, canonical_value):
        """Re-scale a slider to a new unit while keeping the physical value."""
        lo, hi, st, _ = UNIT_SLIDER[unit]
        new_val = round(canonical_value * UNIT_FACTOR[unit] / st) * st
        new_val = min(max(new_val, lo), hi)
        # widen first so the value is never clipped mid-change, then tighten
        sl.min, sl.max = min(sl.min, lo), max(sl.max, hi)
        sl.step = st
        sl.value = new_val
        sl.min, sl.max = lo, hi

    def on_height_unit(change):
        canonical = s_h.value / UNIT_FACTOR[change["old"]]   # value in ft
        state["busy"] = True
        set_slider_range(s_h, change["new"], canonical)
        state["busy"] = False
        update()

    def on_pressure_unit(change):
        canonical = s_p.value / UNIT_FACTOR[change["old"]]   # value in psi
        state["busy"] = True
        set_slider_range(s_p, change["new"], canonical)
        state["busy"] = False
        update()

    def on_velocity_unit(change):
        canonical = s_v.value / UNIT_FACTOR[change["old"]]   # value in m/s
        state["busy"] = True
        set_slider_range(s_v, change["new"], canonical)
        state["busy"] = False
        if tg_mode.value == "Initial velocity":
            update()
        else:
            refresh_value_labels()

    def on_mode(change):
        """Show the velocity controls only in 'Initial velocity' mode."""
        is_vel = change["new"] == "Initial velocity"
        vel_block.layout.display = "flex" if is_vel else "none"
        title_h.value = ("<b>Start height (bottom of ball)</b>" if is_vel
                         else "<b>Drop height (bottom of ball)</b>")
        update()

    def on_velocity(change):
        if tg_mode.value == "Initial velocity":   # slider is hidden otherwise
            update()

    # ----------------------------------------------------------- buttons
    def on_calibrate(_):
        p = current_params()
        lam = round(calibrate_damping(p) / 25.0) * 25.0
        if s_lam.value == lam:
            update()                                 # value unchanged -> no event fires
        else:
            s_lam.value = lam                        # triggers update() via observe

    def on_reset(_):
        state["busy"] = True                         # avoid many redraws
        set_slider_range(s_h, tg_h.value, RANGES["drop_height_ft"][1])
        set_slider_range(s_p, tg_p.value, RANGES["pressure_psi"][1])
        for s, key in [(s_m, "mass_g"), (s_c, "circumference_mm"),
                       (s_S, "shore_A"), (s_t, "shell_thick_mm")]:
            s.value = RANGES[key][1]
        cb_drag.value = True
        set_slider_range(s_v, tg_v.value, RANGES["v_init_m_s"][1])
        tg_mode.value = "Drop from rest"
        state["busy"] = False
        on_calibrate(None)

    def on_export(_):
        with out_csv:
            out_csv.clear_output()
            if state["res"] is None:
                print("Run a simulation first.")
                return
            files = export_csv(state["res"], txt_file.value,
                               state["sweep"] if cb_sweep_csv.value else None)
            for f in files:
                print("Saved:", f)
                try:                                  # clickable link in Jupyter
                    display(FileLink(os.path.relpath(f)))
                except ValueError:                    # e.g. file on another Windows drive
                    pass

    # --------------------------------------------------------- wire it all up
    for s in (s_h, s_p, s_m, s_c, s_S, s_t, s_lam, cb_drag):
        s.observe(update, names="value")
    tg_h.observe(on_height_unit, names="value")
    tg_p.observe(on_pressure_unit, names="value")
    tg_v.observe(on_velocity_unit, names="value")
    tg_mode.observe(on_mode, names="value")
    s_v.observe(on_velocity, names="value")
    btn_cal.on_click(on_calibrate)
    btn_reset.on_click(on_reset)
    btn_csv.on_click(on_export)

    sec = lambda txt: w.HTML(f"<div style='font-size:14px;font-weight:600;"
                             f"border-bottom:1px solid #ccc;margin:6px 0 4px'>{txt}</div>")
    # height block built by hand so its title can change with the mode
    h_block = w.VBox([w.HBox([title_h, tg_h],
                             layout=w.Layout(justify_content="space-between",
                                             align_items="center", width="100%")),
                      s_h, v_h], layout=w.Layout(width="100%", margin="0 0 6px 0"))
    vel_block = block("Starting velocity (+ up / − down)", s_v, v_v, tg_v)
    vel_block.layout.display = "none"            # hidden until the mode is chosen
    left = column([
        sec("Initial conditions"),
        tg_mode,
        h_block, vel_block, cb_drag,
        sec("Inflation"),
        block("Inflation pressure (gauge)", s_p, v_p, tg_p),
        sec("Ball size &amp; mass (FIBA size 7 range)"),
        block("Ball mass", s_m, v_m),
        block("Circumference", s_c, v_c),
    ])
    right = column([
        sec("Viscoelastic rubber shell"),
        block("Rubber hardness", s_S, v_S),
        block("Shell (wall) thickness", s_t, v_t),
        block("Damping λ (Hunt–Crossley)", s_lam, v_lam),
        w.HBox([btn_cal, btn_reset], layout=w.Layout(flex_flow="row wrap")),
        sec("Export"),
        w.HBox([txt_file, btn_csv], layout=w.Layout(flex_flow="row wrap")),
        cb_sweep_csv, out_csv,
    ])
    controls = w.Box([left, right], layout=w.Layout(display="flex", flex_flow="row wrap",
                                                     width="100%"))
    panel = w.VBox([w.HTML(f"<h3 style='margin-bottom:2px'>Basketball drop-impact simulator "
                           f"<span style='font-size:12px;color:#888'>v{__version__}</span></h3>"),
                    controls, out_info, out_plot, out_sweep],
                   layout=w.Layout(width="100%"))
    display(panel)

    on_calibrate(None)   # start with damping tuned to a legal ball (DEFAULT_TEST)
    return panel


# =============================================================================
# SECTION 12 - ENTRY POINT
# In Jupyter: launches the panel.  From a plain terminal: prints a summary and
# saves the figures as PNG files instead.
# =============================================================================
def _in_notebook():
    try:
        from IPython import get_ipython
        return get_ipython() is not None and "IPKernelApp" in get_ipython().config
    except Exception:
        return False


if __name__ == "__main__":
    if _in_notebook():
        app = launch()
    else:
        p = default_params()
        r = simulate_drop(p)
        print(f"lambda = {r['ball']['lambda']:.0f} N s/m^2,  contact {r['contact_time_ms']:.2f} ms,"
              f"  peak {r['peak_force_N']:.0f} N,  COR {r['COR']:.3f}")
        plot_bounce(r).savefig("bounce.png", dpi=110)
        sw = sweep_ranges(p)
        print(sw.round(3).to_string(index=False))
        plot_sweep(sw, sweep_current_ball(r["params"]), p["drop_height_ft"]).savefig("sweep.png", dpi=110)


# =============================================================================
# HOW TO EXTEND THIS TOOL (ideas for the student)
# -----------------------------------------------------------------------------
# * Different ball:   edit RANGES (Section 3), e.g. size 6 = 724-737 mm,
#                     510-567 g.  Volleyball/soccer: change pressure range too.
# * Compliant floor:  add a floor spring k_f (and mass) in series.  Simplest:
#                     make the floor a second mass-spring; add its x_f, v_f to
#                     the state vector in _contact_ode and use the relative
#                     deflection (x - x_f) in floor_force_parts().
# * Different damping law: change the F_damp line in floor_force_parts().
#                     e.g. plain dashpot  F_damp = c * v  (Kelvin-Voigt), or
#                     Hunt-Crossley with exponent n:  lambda * x**n * v.
# * Isothermal gas:   set GAMMA_AIR = 1.0 to see how much the adiabatic
#                     heating stiffens the ball (gas_energy needs the log form
#                     then: P0 V0 ln(V0/V) - P_atm (V0-V)).
# * Compare with data: drop a ball on a force plate or piezo disk [Cross],
#                     export the CSV, and overlay your measured force curve on
#                     F_floor_N vs time_ms.
# =============================================================================
