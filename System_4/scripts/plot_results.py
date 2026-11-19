#!/usr/bin/env python3
"""Render the System 4 result plots as SVG from results/csv (stdlib only).

usage: plot_results.py [results_dir=results]
Writes results/plots/*.svg.  Reads only CSV produced by build/bin/run_analysis;
it never invents data.  A missing CSV is reported and skipped.
"""
import csv, math, os, sys

RES = sys.argv[1] if len(sys.argv) > 1 else "results"
CSV = os.path.join(RES, "csv")
OUT = os.path.join(RES, "plots")
os.makedirs(OUT, exist_ok=True)

COL = {"float": "#222222", "W16": "#d95f02", "W24": "#1b9e77", "W32": "#7570b3",
       "truth": "#999999", "a": "#e7298a", "b": "#66a61e", "c": "#e6ab02"}
W, H, ML, MR, MT, MB = 720, 420, 78, 20, 40, 56


def read(name):
    p = os.path.join(CSV, name)
    if not os.path.exists(p):
        print("skip (missing):", p)
        return None
    with open(p) as f:
        return list(csv.DictReader(f))


def nice(lo, hi, log=False):
    if log:
        lo = 10 ** math.floor(math.log10(lo))
        hi = 10 ** math.ceil(math.log10(hi))
        ticks = [10 ** e for e in range(int(math.log10(lo)), int(math.log10(hi)) + 1)]
        return lo, hi, ticks
    if hi == lo:
        hi = lo + 1.0
    step = 10 ** math.floor(math.log10((hi - lo) / 5.0))
    for m in (1, 2, 5, 10):
        if (hi - lo) / (m * step) <= 6:
            step *= m
            break
    lo = math.floor(lo / step) * step
    hi = math.ceil(hi / step) * step
    n = int(round((hi - lo) / step))
    return lo, hi, [lo + i * step for i in range(n + 1)]


def fmt(v):
    if v == 0:
        return "0"
    if abs(v) >= 1e4 or abs(v) < 1e-2:
        return "%.0e" % v
    return ("%.3g" % v)


def svg_open(title):
    return ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
            'font-family="sans-serif" font-size="11">' % (W, H, W, H),
            '<rect width="100%" height="100%" fill="#ffffff"/>',
            '<text x="%d" y="22" font-size="14" font-weight="bold">%s</text>' % (ML, title)]


def axes(s, xt, yt, xlabel, ylabel, X, Y, logy):
    for t in xt:
        x = X(t)
        s.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e5e5e5"/>' % (x, MT, x, H - MB))
        s.append('<text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (x, H - MB + 16, fmt(t)))
    for t in yt:
        y = Y(t)
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e5e5e5"/>' % (ML, y, W - MR, y))
        s.append('<text x="%d" y="%.1f" text-anchor="end">%s</text>' % (ML - 6, y + 4, fmt(t)))
    s.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#444"/>' % (ML, MT, W - ML - MR, H - MT - MB))
    s.append('<text x="%d" y="%d" text-anchor="middle">%s</text>' % ((ML + W - MR) // 2, H - 12, xlabel))
    s.append('<text transform="translate(16,%d) rotate(-90)" text-anchor="middle">%s</text>' % ((MT + H - MB) // 2, ylabel))


def legend(s, names):
    x = ML + 8
    for n, c in names:
        s.append('<rect x="%d" y="%d" width="14" height="3" fill="%s"/>' % (x, MT + 8, c))
        s.append('<text x="%d" y="%d">%s</text>' % (x + 18, MT + 13, n))
        x += 28 + 7 * len(n)


def line_plot(fname, title, xlabel, ylabel, series, logy=False, dashed=()):
    """series: list of (name, color, xs, ys)"""
    xs_all = [x for _, _, xs, _ in series for x in xs]
    ys_all = [y for _, _, _, ys in series for y in ys if (y > 0 or not logy)]
    if not xs_all or not ys_all:
        print("skip (no data):", fname)
        return
    if logy:
        floor = max(min(ys_all), 1e-12)
        yl, yh, yt = nice(floor, max(ys_all), True)
    else:
        yl, yh, yt = nice(min(ys_all), max(ys_all))
    xl, xh, xt = nice(min(xs_all), max(xs_all))
    X = lambda v: ML + (v - xl) / (xh - xl) * (W - ML - MR)
    if logy:
        Y = lambda v: H - MB - (math.log10(max(v, yl)) - math.log10(yl)) / (math.log10(yh) - math.log10(yl)) * (H - MT - MB)
    else:
        Y = lambda v: H - MB - (v - yl) / (yh - yl) * (H - MT - MB)
    s = svg_open(title)
    axes(s, xt, yt, xlabel, ylabel, X, Y, logy)
    for name, color, xs, ys in series:
        pts = " ".join("%.1f,%.1f" % (X(x), Y(y)) for x, y in zip(xs, ys) if (y > 0 or not logy))
        dash = ' stroke-dasharray="5,3"' if name in dashed else ""
        s.append('<polyline fill="none" stroke="%s" stroke-width="1.4"%s points="%s"/>' % (color, dash, pts))
    legend(s, [(n, c) for n, c, _, _ in series])
    s.append("</svg>")
    open(os.path.join(OUT, fname), "w").write("\n".join(s))
    print("wrote", os.path.join(OUT, fname))


def bar_plot(fname, title, ylabel, cats, series, logy=True):
    """cats: list of labels; series: list of (name,color,values)"""
    vals = [v for _, _, vs in series for v in vs if v > 0]
    if not vals:
        print("skip (no data):", fname)
        return
    yl, yh, yt = nice(min(vals), max(vals), True) if logy else nice(0, max(vals))
    Y = (lambda v: H - MB - (math.log10(max(v, yl)) - math.log10(yl)) / (math.log10(yh) - math.log10(yl)) * (H - MT - MB)) if logy \
        else (lambda v: H - MB - (v - yl) / (yh - yl) * (H - MT - MB))
    s = svg_open(title)
    for t in yt:
        y = Y(t)
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e5e5e5"/>' % (ML, y, W - MR, y))
        s.append('<text x="%d" y="%.1f" text-anchor="end">%s</text>' % (ML - 6, y + 4, fmt(t)))
    s.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#444"/>' % (ML, MT, W - ML - MR, H - MT - MB))
    s.append('<text transform="translate(16,%d) rotate(-90)" text-anchor="middle">%s</text>' % ((MT + H - MB) // 2, ylabel))
    gw = (W - ML - MR) / len(cats)
    bw = gw * 0.8 / len(series)
    for i, c in enumerate(cats):
        gx = ML + i * gw + gw * 0.1
        s.append('<text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (gx + gw * 0.4, H - MB + 16, c))
        for j, (n, col, vs) in enumerate(series):
            v = vs[i]
            if v <= 0:
                s.append('<text x="%.1f" y="%d" text-anchor="middle" font-size="9">0</text>' % (gx + j * bw + bw / 2, H - MB - 3))
                continue
            y = Y(v)
            s.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>' % (gx + j * bw, y, bw - 2, H - MB - y, col))
    legend(s, [(n, c) for n, c, _ in series])
    s.append("</svg>")
    open(os.path.join(OUT, fname), "w").write("\n".join(s))
    print("wrote", os.path.join(OUT, fname))


def col(rows, key):
    return [float(r[key]) for r in rows]


# ---- 1-4: time series of the dynamic scenario (all configs) -------------------
SCEN = "dynamic"
ts = {c: read("%s_%s_timeseries.csv" % (SCEN, c)) for c in ("W16", "W24", "W32")}
if all(ts.values()):
    t = col(ts["W32"], "t")
    for key, ttl, yl, fn in (("pN", "Position North", "m", "plot1_position"),
                             ("vN", "Velocity North", "m/s", "plot2_velocity"),
                             ("roll", "Roll", "rad", "plot3_attitude")):
        fl = "fl_" + key
        series = [("float", COL["float"], t, col(ts["W32"], fl))]
        for c in ("W16", "W24", "W32"):
            series.append((c, COL[c], t, col(ts[c], "fx_" + key)))
        if key in ("pN", "vN"):
            series.insert(0, ("truth", COL["truth"], t, col(ts["W32"], "truth_" + key)))
        line_plot(fn + ".svg", "Float vs Fixed -- %s (%s scenario)" % (ttl, SCEN), "time [s]", yl, series, dashed=("truth",))
    ser = []
    for c in ("W16", "W24", "W32"):
        ser.append((c, COL[c], t, [max(v, 1e-12) for v in col(ts[c], "err_pos")]))
    line_plot("plot4a_abs_error_position.svg", "|position error| fixed vs float (%s)" % SCEN, "time [s]", "m", ser, logy=True)
    ser = []
    for c in ("W16", "W24", "W32"):
        ser.append((c, COL[c], t, [max(v, 1e-12) for v in col(ts[c], "err_att")]))
    line_plot("plot4b_abs_error_attitude.svg", "attitude error fixed vs float (%s)" % SCEN, "time [s]", "rad", ser, logy=True)
    ser = []
    for c in ("W16", "W24", "W32"):
        ser.append((c, COL[c], t, [max(v, 1e-12) for v in col(ts[c], "err_vel")]))
    line_plot("plot4c_abs_error_velocity.svg", "|velocity error| fixed vs float (%s)" % SCEN, "time [s]", "m/s", ser, logy=True)

# ---- 5: RMSE comparison 16/24/32 bit -----------------------------------------
rm = read("summary_rmse.csv")
if rm:
    for sc in ("constvel", "dynamic", "highrange"):
        groups = ["position", "velocity", "attitude_angle", "accel_bias", "gyro_bias", "P_norm_diag"]
        ser = []
        for c in ("W16", "W24", "W32"):
            vals = []
            for g in groups:
                v = [float(r["rmse"]) for r in rm if r["scenario"] == sc and r["config"] == c
                     and r["kind"] == "fixed_vs_float" and r["group"] == g]
                vals.append(v[0] if v else 0.0)
            ser.append((c, COL[c], vals))
        bar_plot("plot5_rmse_%s.svg" % sc, "RMSE fixed vs float, 16/24/32-bit error path (%s)" % sc,
                 "RMSE (log scale; units per group)", groups, ser)

# ---- 6: quantization error histogram -----------------------------------------
qh = read("quant_hist.csv")
if qh:
    for mode in ("round", "truncate"):
        fmts = []
        for r in qh:
            if r["format"] not in fmts:
                fmts.append(r["format"])
        series = []
        cols = [COL["W16"], COL["W24"], COL["W32"]]
        for k, f in enumerate(fmts):
            rows = [r for r in qh if r["format"] == f and r["mode"] == mode]
            xs = [(float(r["bin_lo_lsb"]) + float(r["bin_hi_lsb"])) / 2 for r in rows]
            n = float(rows[0]["samples"]) if rows else 1.0
            ys = [float(r["count"]) / n * 100.0 for r in rows]
            series.append((f, cols[k % 3], xs, ys))
        line_plot("plot6_quant_hist_%s.svg" % mode,
                  "Quantization error histogram (%s) - error in LSB" % mode, "error [LSB]", "% of samples", series)

# ---- 7: overflow / saturation statistics -------------------------------------
cn = read("summary_counters.csv")
if cn:
    cats, sat, ovf, und = [], [], [], []
    for sc in ("static", "constvel", "dynamic", "highrange"):
        for c in ("W16", "W24", "W32"):
            rows = [r for r in cn if r["scenario"] == sc and r["config"] == c]
            cats.append("%s/%s" % (sc[:5], c))
            sat.append(sum(int(r["saturation"]) for r in rows))
            ovf.append(sum(int(r["overflow"]) for r in rows))
            und.append(sum(int(r["underflow"]) for r in rows))
    bar_plot("plot7_counters.svg", "Saturation / overflow / underflow counts per run", "count (log scale; 0 shown as '0')",
             cats, [("saturation", COL["a"], sat), ("overflow", COL["b"], ovf), ("underflow", COL["c"], und)])
