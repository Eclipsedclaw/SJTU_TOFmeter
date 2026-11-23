#!/usr/bin/env python3
"""Muonic X-ray line energies for TOFmeterEmCaptureCascade's kMeasuredLines (/tof/physics/muonXrayEnergies measured).

2026-11-10 for the steel study (C, N, Si, P, S, K, Cr, Mn, Fe, Ni, Mo, Th, U; then in JC/tmp/muonic_refs), extended
2026-11-23 for the muonic X-ray library to every element from H to Mo (the user: measured energies for all 42).

Point-nucleus Dirac energies (reduced mass) of 1s, 2p1/2, 2p3/2, 3d3/2, 3d5/2, 4f5/2, 4f7/2, then per line a factor
R = measured / point Dirac taken from measured lines (anchors) and interpolated in Z for the elements without a
measurement. R holds the finite nuclear size (lowers the lines) and the vacuum polarisation (raises them).
  K lines: Ka1 = 2p3/2 -> 1s, Ka2 = 2p1/2 -> 1s; L lines: La1 = 3d5/2 -> 2p3/2, La2 = 3d3/2 -> 2p1/2,
  La3 = 3d3/2 -> 2p3/2 (R of La1); M lines (4f -> 3d): the 3d and 4f levels barely feel the nucleus, R of La1 at Fe.
Anchors: D.F. Measday, Phys. Rep. 354 (2001) 243, Table 3.4 (natural isotopes' centre of gravity unless an isotope is
named; since 2026-11-23 also O, Na and Cl, and O's 3d-2p); 98Mo from ENSDF (NDS 164, 1 (2020)); 64Ni from ENSDF
(1976Sh21); 238U read off the muX collaboration's spectrum (A. Doinaki, CHIPP Winter School 2025).
Interpolation of R in Z:
  linear (2026-11-10, the rows of C N Si P S K Cr Mn Fe Ni Mo, kept as they were so that the steel runs stay
    comparable): log(1 - R) linear between the two neighbouring anchors (R itself where R > 1), the trend of the last
    two carried on beyond the anchors;
  cubic (2026-11-23, the 31 added elements): the same through the four nearest anchors (two on each side where there
    are two; Lagrange). Leaving each anchor out, its K lines come back within 0.5 keV up to Ca, +1.3 Fe, -3.5 Ni,
    +5.1 Y, -3.1 Mo (linear: +6.6, +4.5, +28.6, +23.6), the L lines within 0.13 keV up to Y.
Not anchors: Measday's values in curly brackets (estimated from systematics: N, Na, Al, Cl 3d-2p) and K's 3d-2p
143.8 keV, which is ~2 keV off the Z^2 trend of its neighbours (Si, Ca; the interpolation gives 141.8).
Fe's L lines in TOFmeter are ENSDF's 56Fe values (265.30, 268.90), not Measday's 265.7, 269.4.

  python3 scripts/muonic_lines.py              the leave-one-out check and every element's lines
  python3 scripts/muonic_lines.py --cpp Z...   kMeasuredLines rows (C++) for these Z
"""
import math
import sys

ALPHA = 1 / 137.035999
MMU = 105.6583755  # MeV
AMU = 931.494102  # MeV
# standard atomic weights; C N Al Si P S K Ca Cr Mn Fe Ni Y Mo and the heavy anchors as in 2026-11-10
MASS = {1: 1.008, 2: 4.0026, 3: 6.94, 4: 9.0122, 5: 10.81, 6: 12.011, 7: 14.007, 8: 15.999, 9: 18.998, 10: 20.180,
        11: 22.990, 12: 24.305, 13: 26.982, 14: 28.086, 15: 30.974, 16: 32.06, 17: 35.45, 18: 39.948, 19: 39.098,
        20: 40.078, 21: 44.956, 22: 47.867, 23: 50.942, 24: 51.996, 25: 54.938, 26: 55.845, 27: 58.933, 28: 58.693,
        29: 63.546, 30: 65.38, 31: 69.723, 32: 72.630, 33: 74.922, 34: 78.971, 35: 79.904, 36: 83.798, 37: 85.468,
        38: 87.62, 39: 88.906, 40: 91.224, 41: 92.906, 42: 95.95, 53: 126.904, 79: 196.967, 82: 207.2,
        83: 208.980, 90: 232.038, 92: 238.029}
NAME = {1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg",
        13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V",
        24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
        35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 53: "I", 79: "Au", 82: "Pb",
        83: "Bi", 90: "Th", 92: "U"}


def binding(Z, n, j):
    """Point-nucleus Dirac binding energy (keV) of level n, j, with the muon's reduced mass."""
    mu = MMU * MASS[Z] * AMU / (MMU + MASS[Z] * AMU)
    za = Z * ALPHA
    k = j + 0.5
    delta = k - math.sqrt(k * k - za * za)
    return (1.0 - 1.0 / math.sqrt(1.0 + (za / (n - delta)) ** 2)) * mu * 1000.0


def point_lines(Z):
    b = {(n, j): binding(Z, n, j) for n, j in [(1, .5), (2, .5), (2, 1.5), (3, 1.5), (3, 2.5), (4, 2.5), (4, 3.5)]}
    return {"Ka1": b[1, .5] - b[2, 1.5], "Ka2": b[1, .5] - b[2, .5],
            "La1": b[2, 1.5] - b[3, 2.5], "La2": b[2, .5] - b[3, 1.5], "La3": b[2, 1.5] - b[3, 1.5],
            "Ma1": b[3, 2.5] - b[4, 3.5], "Ma2": b[3, 1.5] - b[4, 2.5], "Ma3": b[3, 2.5] - b[4, 2.5]}


# measured (keV): (Ka1, Ka2) or (centre of gravity, None); (La1, La2) or (centre of gravity, None)
MEAS_K = {6: (75.2588, None), 7: (102.403, None), 8: (133.535, None), 11: (250.229, None), 13: (346.828, None),
          14: (400.177, None), 17: (578.6, None), 19: (712.69, None), 20: (784.15, 782.7), 26: (1257.19, 1253.06),
          28: (1425.23 + 6.6, 1419.71 + 6.6), 39: (2439.4, 2420.1), 42: (2708.45 + 11.6, 2683.37 + 11.6),
          53: (3723.74, 3667.36), 79: (5760.79, 5591.71), 82: (5962.9, 5778.1), 83: (6034.0, 5843.2),
          92: (6420.0, 6133.0)}
MEAS_L = {6: (13.966, None), 8: (24.915, None), 14: (76.723, None), 20: (156.83, 158.17), 26: (265.7, 269.4),
          39: (599.4, 616.4), 53: (1101.8, 1150.4), 79: (2341.2, 2474.2), 82: (2500.59, 2642.33), 83: (2549.9, 2700.5)}
NOTE = {28: "64Ni + 6.6 keV for natural Ni (rms radius 3.788 against 3.859 fm)",
        42: "98Mo + 11.6 keV for natural Mo (rms radius 4.381 against 4.409 fm)",
        92: "238U read off a plot: 2p1/2 group 6125-6145, 2p3/2 group 6370-6470 keV (deformed nucleus), +-15/+-40 keV"}


def ratios(meas, k1, k2, w1, w2):
    out = {}
    for Z, (m1, m2) in meas.items():
        p = point_lines(Z)
        if m2 is None:
            r = m1 / ((w1 * p[k1] + w2 * p[k2]) / (w1 + w2))
            out[Z] = (r, r)
        else:
            out[Z] = (m1 / p[k1], m2 / p[k2])
    return out


def interp(table, Z, i):
    if Z in table:
        return table[Z][i]
    zs = sorted(table)
    lo = max([z for z in zs if z < Z], default=None)
    hi = min([z for z in zs if z > Z], default=None)
    if lo is None or hi is None:   # beyond the anchors: carry the trend of the last two on
        a, b = (zs[0], zs[1]) if lo is None else (zs[-2], zs[-1])
        lo, hi = a, b
    f = (Z - lo) / (hi - lo)
    a, b = table[lo][i], table[hi][i]
    if a < 1 and b < 1:
        return 1 - math.exp((1 - f) * math.log(1 - a) + f * math.log(1 - b))
    return a + f * (b - a)


def interp4(table, Z, i):
    """R at Z through the four nearest anchors (two on each side where there are two), Lagrange in log(1 - R)."""
    if Z in table:
        return table[Z][i]
    below = [z for z in sorted(table) if z < Z][-2:]
    above = [z for z in sorted(table) if z > Z][:2]
    pts = below + above
    if len(pts) < 2:
        return interp(table, Z, i)
    use_log = all(table[z][i] < 1 for z in pts)
    total = 0.0
    for a in pts:
        w = 1.0
        for b in pts:
            if b != a:
                w *= (Z - b) / (a - b)
        total += w * (math.log(1 - table[a][i]) if use_log else table[a][i])
    return 1 - math.exp(total) if use_log else total


RK = ratios(MEAS_K, "Ka1", "Ka2", 2, 1)
RL = ratios(MEAS_L, "La1", "La2", 9, 5)


def lines(Z, cubic=True):
    p = point_lines(Z)
    f = interp4 if cubic else interp
    k1, k2, l1, l2 = f(RK, Z, 0), f(RK, Z, 1), f(RL, Z, 0), f(RL, Z, 1)
    m = RL[26][0]
    return {"Ka1": p["Ka1"] * k1, "Ka2": p["Ka2"] * k2, "La1": p["La1"] * l1, "La2": p["La2"] * l2,
            "La3": p["La3"] * l1, "Ma1": p["Ma1"] * m, "Ma2": p["Ma2"] * m, "Ma3": p["Ma3"] * m}


def cpp_row(Z):
    """A kMeasuredLines row in the table's layout, with where its K and L lines come from."""
    e = lines(Z)
    if Z <= 5:
        how = ", point Dirac"   # C and N's ratio carried on: within 0.02% of the point-nucleus lines
    elif Z in MEAS_K:
        how = ", K, L measured" if Z in MEAS_L else ", K measured"
    else:
        how = ""                # interpolated (cubic)
    row = f"  {{{Z:2d}, {{{e['Ka1']:7.2f}, {e['Ka2']:7.2f}}}, {{{e['La1']:7.2f}, {e['La2']:7.2f}, {e['La3']:7.2f}}}, {{}}}},"
    return f"{row:<85s}// {NAME[Z]} (2026-11-23{how})"


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--cpp":
        for z in sys.argv[2:]:
            print(cpp_row(int(z)))
        sys.exit(0)
    print("check: anchors left out one at a time (interpolated from the others, cubic | linear) against their measurement")
    for Z in [8, 11, 13, 14, 17, 19, 20, 26, 28, 39, 42, 53, 82]:
        for table, meas, keys, w in [(RK, MEAS_K, ("Ka1", "Ka2"), (2, 1)), (RL, MEAS_L, ("La1", "La2"), (9, 5))]:
            if Z not in meas:
                continue
            saved = table.pop(Z)
            p = point_lines(Z)
            out = []
            for f in (interp4, interp):
                e1, e2 = p[keys[0]] * f(table, Z, 0), p[keys[1]] * f(table, Z, 1)
                m1, m2 = meas[Z]
                if m2 is None:   # a centre of gravity
                    cog = (w[0] * e1 + w[1] * e2) / (w[0] + w[1])
                    out.append(f"centre {cog:8.2f} against {m1:8.2f} ({cog - m1:+7.2f})")
                else:
                    out.append(f"{e1:8.2f} {e2:8.2f} against {m1:8.2f} {m2:8.2f} ({e1 - m1:+7.2f}, {e2 - m2:+7.2f})")
            table[Z] = saved
            print(f"  {NAME[Z]:2s} {keys[0][:2]}: " + " | ".join(out))
    print("\nlines (keV), cubic interpolation (TOFmeter keeps the 2026-11-10 linear rows of C N Si P S K Cr Mn Fe Ni Mo Th U)")
    for Z in sorted(z for z in NAME if z <= 42 or z in (90, 92)):
        e = lines(Z)
        print(f"  {NAME[Z]:2s} Z={Z:2d}  Ka1 {e['Ka1']:8.2f}  Ka2 {e['Ka2']:8.2f} | La1 {e['La1']:8.2f}  La2 {e['La2']:8.2f}"
              f"  La3 {e['La3']:8.2f} | Ma1 {e['Ma1']:7.2f}  Ma2 {e['Ma2']:7.2f}  Ma3 {e['Ma3']:7.2f}")
