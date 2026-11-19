#!/usr/bin/env python3
"""Static check of the System 4 coding rules (ICD section 4 + task brief).

 * no malloc/calloc/realloc/free anywhere in src/ include/ tools/ tests/
 * the fixed-point CORE files contain no float/double type and no libm call,
   except between /* BEGIN-CONVERSION */ and /* END-CONVERSION */ markers
   (explicit float<->fixed I/O helpers and one-time parameter conversion).
Exit code 1 when a rule is violated.
"""
import re, sys, pathlib

root = pathlib.Path(__file__).resolve().parent.parent
core = ["src/core/ekf_fixed_math.c", "src/core/ekf_fixed_divider.c", "src/core/ekf_fixed_trig.c",
        "src/model/ekf_s4_ref_fixed.c"]
forbid_core = re.compile(r"\b(float|double|sinf?|cosf?|tanf?|atanf?|atan2f?|sqrtf?|ldexpf?|floorf?|"
                         r"powf?|logf?|expf?|fabsf?|roundf?)\b")
forbid_mem = re.compile(r"\b(malloc|calloc|realloc|free)\s*\(")

def strip_comments(t):
    t = re.sub(r"/\*.*?\*/", lambda m: "\n" * m.group(0).count("\n"), t, flags=re.S)
    return re.sub(r"//[^\n]*", "", t)

def strip_regions(t):
    # work on the raw text so the markers are still present
    return re.sub(r"/\*\s*BEGIN-CONVERSION\s*\*/.*?/\*\s*END-CONVERSION\s*\*/",
                  lambda m: "\n" * m.group(0).count("\n"), t, flags=re.S)

bad = 0
for d in ("src", "include", "tools", "tests"):
    for f in sorted((root / d).rglob("*.[ch]")):
        txt = strip_comments(f.read_text())
        for n, line in enumerate(txt.splitlines(), 1):
            if forbid_mem.search(line):
                print(f"VIOLATION {f.relative_to(root)}:{n}: dynamic memory: {line.strip()}")
                bad += 1
for rel in core:
    raw = (root / rel).read_text()
    txt = strip_comments(strip_regions(raw))
    for n, line in enumerate(txt.splitlines(), 1):
        if line.lstrip().startswith("#"):
            continue
        m = forbid_core.search(line)
        if m:
            print(f"VIOLATION {rel}:{n}: '{m.group(0)}' in fixed-point core: {line.strip()}")
            bad += 1
print("core rules: %s" % ("FAILED (%d violations)" % bad if bad else "OK"))
sys.exit(1 if bad else 0)
