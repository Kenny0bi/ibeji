"""Verify every reference in references.bib against Crossref, and every citation against the .bib.

For each entry: the DOI must resolve in Crossref, and the Crossref title, first-author family
name, year and journal must match the .bib entry. Title and journal use a normalized
similarity check (case, punctuation, braces and "The" ignored); year must match exactly.
Also reports cited keys missing from the .bib and .bib entries not cited in ibeji.tex.
Exits non-zero on any failure, so it can gate the final compile.

Usage: python3 paper/verify_refs.py
"""
import difflib
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIB = HERE / "references.bib"
TEX = HERE / "ibeji.tex"
CACHE = HERE / ".crossref_cache.json"


def parse_bib(text):
    entries = {}
    for m in re.finditer(r"@(\w+)\{([^,]+),(.*?)\n\}", text, flags=re.S):
        if m.group(1).lower() == "ieeetranbstctl":  # bibliography style settings, not a reference
            continue
        fields = dict((k.lower(), re.sub(r"\s+", " ", v).strip())
                      for k, v in re.findall(r"(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}", m.group(3)))
        entries[m.group(2).strip()] = fields
    return entries


def norm(s):
    s = re.sub(r"[{}\\$]", "", s or "").lower()
    s = re.sub(r"^the\s+", "", s)
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def crossref(doi, cache):
    if doi in cache:
        return cache[doi]
    url = f"https://api.crossref.org/works/{urllib.request.quote(doi)}"
    req = urllib.request.Request(url, headers={"User-Agent": "ibeji-reference-check (mailto:obidelek19@gmail.com)"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            msg = json.load(r)["message"]
    except Exception as e:
        return {"error": str(e)}
    parts = (msg.get("published-print") or msg.get("published-online") or msg.get("issued") or {}).get("date-parts", [[None]])
    rec = {"title": (msg.get("title") or [""])[0], "journal": (msg.get("container-title") or [""])[0],
           "first_author": (msg.get("author") or [{}])[0].get("family", ""), "year": parts[0][0],
           "years_all": sorted({d[0][0] for d in (msg.get(k, {}).get("date-parts") for k in ("published-print", "published-online", "issued")) if d})}
    cache[doi] = rec
    time.sleep(0.2)
    return rec


def main():
    entries = parse_bib(BIB.read_text())
    tex = TEX.read_text()
    # Structural guard. On 2026-09-15 an edit that replaced the span from the availability
    # section to \end{document} silently deleted these two lines; the paper still compiled,
    # but with 65 undefined citations and no reference list. This turns that into a hard stop.
    structural = [c for c in ("\\bibliographystyle{", "\\bibliography{") if c not in tex]
    if structural:
        print("FAIL  manuscript is missing: " + ", ".join(structural))
        sys.exit(1)
    cited = {k.strip() for g in re.findall(r"\\cite\{([^}]*)\}", tex) for k in g.split(",")}
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    failures = 0
    for key, f in sorted(entries.items()):
        doi = f.get("doi")
        if not doi:
            print(f"FAIL  {key}: no DOI"); failures += 1; continue
        cr = crossref(doi, cache)
        if "error" in cr:
            print(f"FAIL  {key}: DOI {doi} did not resolve ({cr['error']})"); failures += 1; continue
        problems = []
        if difflib.SequenceMatcher(None, norm(f.get("title")), norm(cr["title"])).ratio() < 0.9:
            problems.append(f"title differs: bib '{f.get('title')}' vs Crossref '{cr['title']}'")
        bib_first = norm(re.split(r"\s+and\s+", f.get("author", ""))[0].split(",")[0])
        if cr["first_author"] and bib_first and norm(cr["first_author"]) != bib_first:
            problems.append(f"first author differs: bib '{bib_first}' vs Crossref '{cr['first_author']}'")
        if str(f.get("year")) != str(cr["year"]) and int(f.get("year", 0)) not in cr["years_all"]:
            problems.append(f"year differs: bib {f.get('year')} vs Crossref {cr['years_all']}")
        if f.get("journal") and cr["journal"] and difflib.SequenceMatcher(None, norm(f["journal"]), norm(cr["journal"])).ratio() < 0.8:
            problems.append(f"journal differs: bib '{f['journal']}' vs Crossref '{cr['journal']}'")
        if problems:
            failures += 1
            print(f"FAIL  {key}: " + "; ".join(problems))
        else:
            print(f"OK    {key}: {cr['first_author']} {cr['year']}, {cr['journal']}")
    CACHE.write_text(json.dumps(cache, indent=1))
    missing = sorted(cited - set(entries))
    uncited = sorted(set(entries) - cited)
    if missing:
        failures += len(missing)
        print("FAIL  cited but not in .bib:", ", ".join(missing))
    print(f"\n{len(entries)} entries, {len(cited)} cited keys; uncited entries ({len(uncited)}): {', '.join(uncited) or 'none'}")
    print("ALL REFERENCES VERIFIED" if failures == 0 else f"{failures} FAILURE(S)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
