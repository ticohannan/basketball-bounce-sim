# References

Reference style: ASME (numbered). The numbers match the bracketed citations in
[`DESIGN.md`](DESIGN.md), [`README.md`](../README.md) and the source code.
DOIs are given as resolvable `https://doi.org/` links. Web pages show the date they were
accessed.

## Rules and standards (ball specifications)

[1] National Basketball Association, 2026, "Rule No. 1: Court Dimensions – Equipment,"
*Official NBA Rulebook 2026–27*, Sec. II(f)(1), accessed Oct. 20, 2026,
https://official.nba.com/rule-no-1-court-dimensions/

[2] Fédération Internationale de Basketball (FIBA), 2026, *Official Basketball Rules 2026 –
Basketball Equipment*, v1.0, FIBA, Mies, Switzerland, Sec. 2 (Basketballs), valid from
Oct. 1, 2026, accessed Oct. 20, 2026,
https://assets.fiba.basketball/image/upload/documents-corporate-fiba-official-rules-2026-equipment-v1-0.pdf

[3] Fédération Internationale de Basketball (FIBA), 2018, *Official Basketball Rules 2018 –
Basketball Equipment*, FIBA, Mies, Switzerland, Sec. 7 (Balls), valid from Oct. 1, 2018,
accessed Oct. 20, 2026,
https://refereeing.fiba.basketball/images/docs/doc/5ebfbbb5a9f17da0f62ab2174cebaa85/FIBABasketballEquipment2018.pdf

[4] National Federation of State High School Associations (NFHS), 2026, *2026–27 NFHS
Basketball Rules Book*, NFHS, Indianapolis, IN, Rule 1-12-2.

## Inflated-ball impact mechanics

[5] Tadrist, L., and Darbois Texier, B., 2017, "Open Questions on the Impact of an Inflated
Ball," arXiv:1708.01282 [cond-mat.soft]. https://doi.org/10.48550/arXiv.1708.01282

[6] Stronge, W. J., and Ashcroft, A. D. C., 2007, "Oblique Impact of Inflated Balls at Large
Deflections," Int. J. Impact Eng., 34(6), pp. 1003–1019.
https://doi.org/10.1016/j.ijimpeng.2006.04.006

[7] Lazarus, A., Florijn, H. C. B., and Reis, P. M., 2012, "Geometry-Induced Rigidity in
Nonspherical Pressurized Elastic Shells," Phys. Rev. Lett., 109(14), p. 144301.
https://doi.org/10.1103/PhysRevLett.109.144301

[8] Gent, A. N., 1958, "On the Relation Between Indentation Hardness and Young's Modulus,"
Rubber Chem. Technol., 31(4), pp. 896–906. https://doi.org/10.5254/1.3542351

[9] Hunt, K. H., and Crossley, F. R. E., 1975, "Coefficient of Restitution Interpreted as
Damping in Vibroimpact," ASME J. Appl. Mech., 42(2), pp. 440–445.
https://doi.org/10.1115/1.3423596

[10] Cross, R., 1999, "The Bounce of a Ball," Am. J. Phys., 67(3), pp. 222–227.
https://doi.org/10.1119/1.19229

[11] Cross, R., 2006, "Bouncing Balls," presentation slides, Physics Department, University
of Sydney, Sydney, Australia, slides 5 and 8, accessed Oct. 20, 2026,
https://www.physics.usyd.edu.au/~cross/BallBounce.htm

[12] Cross, R., 2014, "Physics of Bounce," School of Physics, University of Sydney, Sydney,
Australia, accessed Oct. 20, 2026, https://www.physics.usyd.edu.au/~cross/BOUNCE.htm

[13] Njock Libii, J., 2011, "Applying Dynamics to the Engineering of the Perfect Bounce:
Experimental Investigation of Why the NBA Requires a Specific Inflation Pressure for
Basketballs Used in Professional Games," *Proc. 2011 ASEE Annual Conference & Exposition*,
Vancouver, BC, Canada, June 26–29, pp. 22.218.1–22.218.13.
https://doi.org/10.18260/1-2--17499

## Ball construction (material and thickness ranges)

[14] Feeney, B., Kennedy, T. J., and Nesbitt, R. D., 2008, "Game Ball With Urethane Gum Rubber
Surface," U.S. Patent Application Publication No. US 2008/0188334 A1.
https://patents.google.com/patent/US20080188334A1/en

[15] Walker, A. D., and Baltronis, J. F., 1994, "Basketball With Polyurethane Cover,"
U.S. Patent No. 5,310,178. https://patents.google.com/patent/US5310178A/en

## Fluid mechanics and numerical methods

[16] White, F. M., 2016, *Fluid Mechanics*, 8th ed., McGraw-Hill Education, New York.
ISBN 978-0-07-339827-3.

[17] Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D.,
Burovski, E., et al., 2020, "SciPy 1.0: Fundamental Algorithms for Scientific Computing in
Python," Nat. Methods, 17(3), pp. 261–272. https://doi.org/10.1038/s41592-019-0686-2

[18] Dormand, J. R., and Prince, P. J., 1980, "A Family of Embedded Runge-Kutta Formulae,"
J. Comput. Appl. Math., 6(1), pp. 19–26. https://doi.org/10.1016/0771-050X(80)90013-3

[19] Brent, R. P., 1973, *Algorithms for Minimization Without Derivatives*, Prentice-Hall,
Englewood Cliffs, NJ. ISBN 0-13-022335-2.

## Further reading (not used directly in the model)

[20] Haron, M. A., Jailani, A., Nik Abdullah, N. A. F., Ismail, R. S., Abd Rahim, S. Z., and
Ghazali, M. F., 2017, "Review of the Dynamic Behaviour of Sports Balls During Normal and
Oblique Impacts," AIP Conf. Proc., 1885, p. 020035. https://doi.org/10.1063/1.5002229

[21] Darbois Texier, B., and Tadrist, L., 2023, "Impact of a Pressurized Membrane:
Coefficient of Restitution," Phys. Rev. E, 107(5), p. 055001.
https://doi.org/10.1103/PhysRevE.107.055001

---

## What each reference is used for

| Ref. | Used for |
|---|---|
| [1] | Inflation pressure range 7.5–8.5 psi (gauge) |
| [2] | Default rebound test used to calibrate damping: 1800 mm drop → 1035–1085 mm rebound (to the underside); size 7 = 750–770 mm, 580–620 g |
| [3] | Legacy rebound test (1200–1400 mm to the top of the ball); wider slider ranges 749–780 mm, 567–650 g |
| [4] | High-school rebound test: 6 ft drop → 49–54 in (to the top of the ball) |
| [5] | Contact-patch area *A(x)*, cap volume *V(x)*, adiabatic gas spring, linear-spring contact time, shell vs. gas energy storage; measured basketball contact time ≈ 11 ms, COR ≈ 0.78; shell *E* ≈ 10 MPa, *t* ≈ 2.5 mm |
| [6] | Spherical-cap deformation geometry; momentum-flux loss mechanism (listed as a model limitation) |
| [7] | Thin-shell indentation stiffness |
| [8] | Shore A hardness → Young's modulus conversion |
| [9] | Nonlinear (Hunt–Crossley) damping law |
| [10]–[12] | Sinusoid-like force pulse; basketball contact time on wood ≈ 15 ms; contact time falls slightly with speed; force-plate / piezo measurement method |
| [13] | Experiments: higher inflation → higher rebound and less energy dissipated |
| [14] | Cover rubber hardness 45–85 Shore A |
| [15] | Layer thicknesses: bladder ≈ 0.8 mm, windings 0.3–0.7 mm, carcass 1.8–2.0 mm, cover panels ≈ 1.8 mm |
| [16] | Quadratic drag, sphere drag coefficient *C*<sub>D</sub> ≈ 0.47 |
| [17], [18] | ODE solver (`scipy.integrate.solve_ivp`, RK45 = Dormand–Prince) |
| [19] | Root finder for damping calibration (`scipy.optimize.brentq`) |
