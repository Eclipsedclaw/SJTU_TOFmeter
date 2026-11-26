#!/usr/bin/env python3
"""
setup_from_xlsx.py  -  a setup file (setups/<name>.toml) made from the setup spreadsheet

The hardware of a setup is kept in a spreadsheet, one per setup (e.g. setups/two_layer_2026-11-26.xlsx).
Its sheet "help" explains every column; in short:

  sheet "channels"  one row per channel, the mapping sheet as it is:
      DAQ CARD ID, DAQ CARD IP    the board (merged cells are fine: a row without one keeps the one above)
      channel ID                  DAQ_YY id = 10000 + 10 x board + channel
      paddle ID                   the paddle's number; in the analysis it is called P<number>
      side                        the end on the cable label, A or B
      geo_layer                   its layer (a row of the "layers" sheet)
      geo_paddle                  its place across the layer: 1 = at first_cm (the simulation's bar 1),
                                  2 = one pitch further, ...
      geo_side                    the end it reads, A or B (A = at the layer's A_end); empty: as side
      x_cm, y_cm, z_cm            optional: the paddle's measured centre, instead of the layer pattern
      time_offset_ns              optional: ns subtracted from this channel's times (a cable delay)
    Rows without a paddle ID (unused channels) are skipped.
  sheet "layers"    one row per layer, the top layer first: layer, z_cm, along (x or y), A_end (-x, +x,
                    -y or +y), pitch_cm, first_cm (centre of geo_paddle 1 across the layer), length_cm,
                    width_cm, thickness_cm
  sheet "setup"     key | value: name, valid_from, description, geometry_source, notes (comment lines for
                    the top of the file), wall and wall_label (the room side marked in the pictures),
                    top_up (the direction drawn upwards in the TOP views: +y, -x, -y or +x; a turn, never a
                    mirror), board_offset_ns_<board>, settings_from

Writes setups/<name>.toml: geometry and channels from the spreadsheet, and the analysis settings ([pulse],
[events], [calibration], [tracks], [flux], ...) copied word for word, comments included, from settings_from
(empty: the file being replaced, else the newest setup file before valid_from). The new file is read back
with the analysis's own reader before it replaces anything, and its TOP view is drawn next to it
(setups/<name>_top_view.png). An existing file is replaced only with --force.

Usage:
    python3 setup_from_xlsx.py                                      # asks for the spreadsheet
    python3 setup_from_xlsx.py setups/two_layer_2026-11-26.xlsx
    python3 setup_from_xlsx.py setups/my_setup.xlsx --force         # replace setups/<name>.toml
    python3 setup_from_xlsx.py setups/my_setup.xlsx --out check.toml   # write it somewhere else

A new setup: copy the spreadsheet, change its rows, give the "setup" sheet a new name and valid_from,
run this script and look at the TOP view. Runs from valid_from on are analysed with the new file.
Needs openpyxl.
"""

import argparse
import datetime
import glob
import json
import os
import re
import sys

import muon_common as mc

REMADE = ("display", "board_offsets_ns", "channel_offsets_ns")     # written from the spreadsheet
ALIASES = {
    "board": ("daq_card_id", "board", "card"), "ip": ("daq_card_ip", "ip"),
    "channel": ("channel_id", "channel"), "paddle": ("paddle_id", "paddle"), "side": ("side",),
    "layer": ("geo_layer", "layer"), "geo_paddle": ("geo_paddle",), "geo_side": ("geo_side",),
    "x_cm": ("x_cm",), "y_cm": ("y_cm",), "z_cm": ("z_cm",), "time_offset_ns": ("time_offset_ns",),
}
LAYER_NUMBERS = ("z_cm", "pitch_cm", "first_cm", "length_cm", "width_cm", "thickness_cm")
HEADER = re.compile(r"^\s*\[\[?([A-Za-z0-9_]+)\]\]?\s*(#.*)?$")


def fail(text):
    sys.exit(f"ERROR: {text}")


def key(text):
    return re.sub(r"[^a-z0-9]+", "_", str(text).strip().lower()).strip("_")


def blank(v):
    return v is None or str(v).strip() == ""


def number(v, what):
    try:
        return float(v)
    except (TypeError, ValueError):
        fail(f"{what}: '{v}' is not a number")


def whole(v, what):
    x = number(v, what)
    if x != int(x):
        fail(f"{what}: '{v}' is not a whole number")
    return int(x)


def sheet(wb, name, first=False):
    for ws in wb.worksheets:
        if key(ws.title) == name:
            return ws
    if first:
        return wb.worksheets[0]
    fail(f"the spreadsheet has no sheet '{name}' (it has: {', '.join(ws.title for ws in wb.worksheets)})")


def rows_of(ws):
    """The rows below the first one as dicts by its normalised headers; empty rows left out."""
    rows = list(ws.iter_rows(values_only=True))
    head = [key(h) if not blank(h) else "" for h in (rows[0] if rows else ())]
    out = []
    for n, row in enumerate(rows[1:], start=2):
        if all(blank(v) for v in row):
            continue
        d = {h: v for h, v in zip(head, row) if h}
        d["_row"] = n
        out.append(d)
    return out


def pick(d, field):
    for name in ALIASES[field]:
        if not blank(d.get(name)):
            return d[name]
    return None


def paddle_name(v):
    """P<number> for a paddle number (29 -> P29); any other text as it is."""
    try:
        x = float(v)
        if x == int(x):
            return f"P{int(x)}"
    except (TypeError, ValueError):
        pass
    return str(v).strip()


def read_channels(ws):
    out, board = [], None
    for d in rows_of(ws):
        where = f"sheet {ws.title}, row {d['_row']}"
        if pick(d, "board") is not None:
            board = whole(pick(d, "board"), f"{where}: DAQ card ID")
        if pick(d, "paddle") is None:
            continue                                            # an unused channel
        if pick(d, "channel") is None or pick(d, "layer") is None:
            fail(f"{where}: paddle {pick(d, 'paddle')} needs a channel ID and a geo_layer")
        chid = whole(pick(d, "channel"), f"{where}: channel ID")
        if chid < 10010 or (board is not None and (chid - 10000) // 10 != board):
            fail(f"{where}: channel {chid} is not a channel of board {board} (id = 10000 + 10 x board + channel)")
        side = str(pick(d, "side") or "").strip().upper()
        end = str(pick(d, "geo_side") or side).strip().upper()
        if end not in ("A", "B"):
            fail(f"{where}: geo_side (or side) must be A or B, not '{end}'")
        r = {"where": where, "channel": chid, "board": (chid - 10000) // 10, "name": paddle_name(pick(d, "paddle")),
             "side": side, "end": end, "layer": str(pick(d, "layer")).strip(),
             "geo": None if pick(d, "geo_paddle") is None else whole(pick(d, "geo_paddle"), f"{where}: geo_paddle")}
        for k in ("x_cm", "y_cm", "z_cm", "time_offset_ns"):
            r[k] = None if pick(d, k) is None else number(pick(d, k), f"{where}: {k}")
        out.append(r)
    if not out:
        fail(f"sheet {ws.title}: no channel with a paddle ID")
    return out


def read_layers(ws):
    layers = {}
    for d in rows_of(ws):
        where = f"sheet {ws.title}, row {d['_row']}"
        name = str(d.get("layer") or "").strip()
        if not name:
            fail(f"{where}: no layer name")
        L = {"name": name, "along": str(d.get("along") or "").strip().lower(),
             "a_end": str(d.get("a_end") or "").strip().lower()}
        for k in LAYER_NUMBERS:
            L[k] = number(d.get(k), f"{where}: {k}")
        if L["along"] not in ("x", "y") or L["a_end"] not in (f"+{L['along']}", f"-{L['along']}"):
            fail(f"{where}: along must be x or y and A_end +{L['along'] or 'x'} or -{L['along'] or 'x'}")
        layers[name] = L
    if not layers:
        fail(f"sheet {ws.title}: no layers")
    return layers


def read_setup(ws):
    out = {}
    for row in ws.iter_rows(values_only=True):
        if row and not blank(row[0]) and key(row[0]) != "key":
            out[key(row[0])] = row[1] if len(row) > 1 else None
    for need in ("name", "valid_from"):
        if blank(out.get(need)):
            fail(f"sheet {ws.title}: '{need}' is missing")
    v = out["valid_from"]
    if isinstance(v, datetime.datetime):
        v = v.date()
    if not isinstance(v, datetime.date):
        try:
            v = datetime.date.fromisoformat(str(v).strip())
        except ValueError:
            fail(f"sheet {ws.title}: valid_from '{out['valid_from']}' is not a date (YYYY-MM-DD)")
    out["valid_from"] = v
    out["name"] = str(out["name"]).strip()
    return out


def build_paddles(channels, layers):
    """Paddles in layer order, then across the layer: name, layer, axis, centre, size, A_end, channels."""
    found, notes = {}, []
    for r in channels:
        if r["layer"] not in layers:
            fail(f"{r['where']}: geo_layer '{r['layer']}' is not in the layers sheet ({', '.join(layers)})")
        p = found.setdefault(r["name"], {"name": r["name"], "layer": r["layer"], "geo": r["geo"], "ends": {},
                                         "xyz": {}, "boards": set()})
        if (r["layer"], r["geo"]) != (p["layer"], p["geo"]):
            fail(f"{r['where']}: {r['name']}'s channels disagree on geo_layer / geo_paddle")
        if r["end"] in p["ends"]:
            fail(f"{r['where']}: {r['name']} has two channels on end {r['end']}")
        p["ends"][r["end"]] = r["channel"]
        p["boards"].add(r["board"])
        for k in ("x_cm", "y_cm", "z_cm"):
            if r[k] is not None:
                if k in p["xyz"] and p["xyz"][k] != r[k]:
                    fail(f"{r['where']}: {r['name']}'s channels disagree on {k}")
                p["xyz"][k] = r[k]
        if r["side"] and r["side"] != r["end"]:
            notes.append(f"{r['name']}: ch{r['channel']} is labelled {r['side']} but reads end {r['end']}")
    taken = {}
    for p in found.values():
        L = layers[p["layer"]]
        if p["geo"] is None and not {"x_cm", "y_cm"} <= set(p["xyz"]):
            fail(f"{p['name']}: needs a geo_paddle (or both x_cm and y_cm)")
        if p["geo"] is not None:
            if (p["layer"], p["geo"]) in taken:
                fail(f"{p['name']} and {taken[(p['layer'], p['geo'])]} are both geo_paddle {p['geo']} of {p['layer']}")
            taken[(p["layer"], p["geo"])] = p["name"]
            across = L["first_cm"] + (p["geo"] - 1) * L["pitch_cm"]
            x, y = (0.0, across) if L["along"] == "x" else (across, 0.0)
        else:
            x = y = 0.0
        p.update(x=p["xyz"].get("x_cm", x), y=p["xyz"].get("y_cm", y), z=p["xyz"].get("z_cm"), L=L)
    order = list(layers)
    paddles = sorted(found.values(), key=lambda p: (order.index(p["layer"]), p["geo"] or 0, p["name"]))
    return paddles, notes


def settings_blocks(path):
    """Every [section] of a setup file except the ones made from the spreadsheet, word for word, each with
    the comment lines just above its header."""
    lines = open(path, encoding="utf-8").read().splitlines()
    heads = []
    for i, line in enumerate(lines):
        m = HEADER.match(line)
        if m:
            j = i
            while j > 0 and lines[j - 1].lstrip().startswith("#"):
                j -= 1
            heads.append((j, m.group(1), line.lstrip().startswith("[[")))
    blocks = []
    for k, (j, name, array) in enumerate(heads):
        stop = heads[k + 1][0] if k + 1 < len(heads) else len(lines)
        if not array and name not in REMADE:
            blocks.append("\n".join(lines[j:stop]).rstrip())
    return blocks


def default_settings(out_path, valid_from):
    if os.path.isfile(out_path):
        return out_path
    older = []
    for path in glob.glob(os.path.join(mc.SETUPS, "*.toml")):
        s = mc.load_setup(path)
        if s.valid_from is None or s.valid_from <= valid_from:
            older.append((s.valid_from or datetime.date.min, path))
    return max(older)[1] if older else None


def toml_number(v):
    text = f"{float(v):.6g}"
    return text if any(c in text for c in ".en") else text + ".0"


def toml_text(setup, layers, paddles, channels, blocks, source):
    q = lambda s: json.dumps(str(s), ensure_ascii=False)                     # noqa: E731 (a TOML basic string)
    out = [f"# Setup {setup['name']}, made by setup_from_xlsx.py from {source['xlsx']}",
           f"#   (spreadsheet of {source['mtime']}, file written {source['now']}).",
           "# Geometry and channels come from the spreadsheet: change them there and run the script again.",
           f"# The analysis settings below come from {source['settings']}; edit them here if needed:",
           "# the script keeps them when it replaces this file. Its TOP view: "
           f"{setup['name']}_top_view.png."]
    if not blank(setup.get("notes")):
        out.append("#")
        out += [f"# {line}".rstrip() for line in str(setup["notes"]).splitlines()]
    out += ["", f"name = {q(setup['name'])}"]
    for k in ("description", "geometry_source"):
        if not blank(setup.get(k)):
            out.append(f"{k} = {q(setup[k])}")
    out.append(f"valid_from = {setup['valid_from'].isoformat()}")
    for L in layers.values():
        out += ["", "[[layers]]", f"name = {q(L['name'])}", f"z_cm = {toml_number(L['z_cm'])}"]
    for p in paddles:
        L = p["L"]
        where = f"geo_paddle {p['geo']}" if p["geo"] is not None else "measured centre"
        out += ["", "[[paddles]]", f"name = {q(p['name'])}".ljust(37) + f"# {p['layer']}, {where}",
                f"layer = {q(p['layer'])}", f"axis = {q(L['along'])}", f"x_cm = {toml_number(p['x'])}",
                f"y_cm = {toml_number(p['y'])}"]
        if p["z"] is not None:
            out.append(f"z_cm = {toml_number(p['z'])}")
        out += [f"length_cm = {toml_number(L['length_cm'])}", f"width_cm = {toml_number(L['width_cm'])}",
                f"thickness_cm = {toml_number(L['thickness_cm'])}", f"A_end = {q(L['a_end'])}"]
        ends = ", ".join(f"{e} = {c}" for e, c in sorted(p["ends"].items()))
        cards = ", ".join(f"board {(c - 10000) // 10} ch{c % 10}" for _, c in sorted(p["ends"].items()))
        out.append(f"channels = {{ {ends} }}".ljust(37) + f"# {cards}")
    out += ["", "# ---- analysis settings ----", ""]
    out += [b + "\n" for b in blocks]
    wall = "" if blank(setup.get("wall")) else str(setup["wall"]).strip()
    label = str(setup.get("wall_label") or ("wall" if wall else "")).strip()
    up = "+y" if blank(setup.get("top_up")) else str(setup["top_up"]).strip()
    out += ["[display]", f"side = {q(wall)}".ljust(29) + "# the room side marked in the pictures (+x, -x, +y, -y)",
            f"side_label = {q(label)}",
            f"top_up = {q(up)}".ljust(29) + "# the direction drawn upwards in the TOP views (a turn, not a mirror)", "",
            "# ns subtracted from every pulse time of a board (clock offsets), from the setup sheet", "[board_offsets_ns]"]
    for k in sorted(k for k in setup if k.startswith("board_offset_ns_")):
        if not blank(setup[k]):
            out.append(f"{whole(k.rsplit('_', 1)[1], k)} = {toml_number(number(setup[k], k))}")
    out += ["", "# ns subtracted from single channels (cable delays), from the channels sheet", "[channel_offsets_ns]"]
    out += [f"{r['channel']} = {toml_number(r['time_offset_ns'])}" for r in channels if r["time_offset_ns"] is not None]
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx", nargs="?", default=None, help="Setup spreadsheet (asked for if left out)")
    ap.add_argument("--out", default=None, help="Setup file to write (default: setups/<name>.toml)")
    ap.add_argument("--force", action="store_true", help="Replace an existing setup file")
    args = ap.parse_args()
    try:
        import openpyxl
    except ImportError:
        fail("reading the spreadsheet needs openpyxl (pip install openpyxl)")
    if args.xlsx:
        path = os.path.abspath(os.path.expanduser(args.xlsx))
        if not os.path.isfile(path):
            fail(f"no such file: {path}")
    else:
        mc.enable_path_completion()
        path = mc.ask("Setup spreadsheet (.xlsx)", convert=mc.existing_file)
    wb = openpyxl.load_workbook(path, data_only=True)
    setup = read_setup(sheet(wb, "setup"))
    layers = read_layers(sheet(wb, "layers"))
    channels = read_channels(sheet(wb, "channels", first=True))
    paddles, notes = build_paddles(channels, layers)
    out = os.path.abspath(os.path.expanduser(args.out)) if args.out else os.path.join(mc.SETUPS, f"{setup['name']}.toml")
    if os.path.exists(out) and not args.force:
        fail(f"{out} exists; --force replaces it (its analysis settings are kept)")
    src = setup.get("settings_from")
    settings = (os.path.join(mc.SETUPS, str(src).strip()) if not blank(src) and not os.path.isabs(str(src).strip())
                else (str(src).strip() if not blank(src) else default_settings(out, setup["valid_from"])))
    if settings and not os.path.isfile(settings):
        fail(f"settings_from: no such setup file {settings}")
    blocks = settings_blocks(settings) if settings else []
    for path_other in glob.glob(os.path.join(mc.SETUPS, "*.toml")):
        other = mc.load_setup(path_other)
        if os.path.abspath(path_other) != out and other.valid_from == setup["valid_from"]:
            notes.append(f"{os.path.basename(path_other)} has the same valid_from: the analysis takes the one "
                         "whose path sorts last; give this setup another date")
    source = {"xlsx": path, "settings": os.path.basename(settings) if settings else "nowhere (the defaults)",
              "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M"),
              "now": datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}
    text = toml_text(setup, layers, paddles, channels, blocks, source)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + ".new"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(text)
    try:
        check = mc.load_setup(tmp)                 # the analysis's own reader: stops here if anything is wrong
    except SystemExit:
        os.remove(tmp)
        raise
    os.replace(tmp, out)
    picture = os.path.splitext(out)[0] + "_top_view.png"
    fig = mc.draw_top_view(check)
    mc.footer(fig, f"from {os.path.basename(path)} by setup_from_xlsx.py  ·  {source['now']}")
    mc.save_figure(fig, picture)
    boards = sorted({b for p in paddles for b in p["boards"]})
    print(f"Setup      : {setup['name']}, valid from {setup['valid_from']}")
    for L in layers.values():
        names = [p["name"] for p in paddles if p["layer"] == L["name"]]
        print(f"  layer {L['name']:<8} z = {L['z_cm']:g} cm, along {L['along']}, A ends at {L['a_end']}: "
              f"{len(names)} paddles {names[0]} ... {names[-1]}" if names else f"  layer {L['name']}: no paddles")
    print(f"Channels   : {len(channels)} on boards {', '.join(map(str, boards))}")
    print(f"Settings   : {source['settings']}")
    for n in notes:
        print(f"Note       : {n}")
    print(f"Wrote      : {out}\n             {picture} (look at it before analysing)")


if __name__ == "__main__":
    main()
