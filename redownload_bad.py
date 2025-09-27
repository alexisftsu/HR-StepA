# redownload_bad.py
import hashlib, os, sys, time, re
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

BASE = "https://beta.lmfdb.org/data/riemann-zeta-zeros/"

def load_md5_map(md5file):
    m = {}
    pat = re.compile(r"^([0-9a-fA-F]{32})\s+\*([^\s]+)$")
    with open(md5file, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            mline = pat.match(line.strip())
            if mline:
                h, fn = mline.group(1).lower(), mline.group(2)
                m[fn] = h
    return m

def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def looks_like_html(b):
    head = b[:2048].lstrip()
    return head.startswith(b"<!DOCTYPE") or head.startswith(b"<html") or b"<script" in head

def download(url, headers, timeout=60):
    req = Request(url, headers=headers, method="GET")
    with urlopen(req, timeout=timeout) as resp:
        data = resp.read()
        return resp.getcode(), resp.headers.get_content_type(), data

def main():
    root = os.getcwd()
    datadir = os.path.join(root, "data")
    md5file = os.path.join(root, "md5_subset.txt")
    badlist = os.path.join(root, "bad_files.txt")

    if not os.path.exists(badlist):
        print("No existe bad_files.txt; ejecuta primero el paso 1.", file=sys.stderr)
        sys.exit(1)

    md5map = load_md5_map(md5file)
    headers = {
        # Evita compresión en tránsito y aparenta navegador
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",
        "Accept": "text/plain, */*;q=0.1",
        "Accept-Encoding": "identity",
        "Referer": "https://beta.lmfdb.org/",
        "Connection": "close",
    }
    ok, fail = 0, 0

    for fn in open(badlist, "rt", encoding="utf-8", errors="ignore"):
        fn = fn.strip()
        if not fn: 
            continue
        url = BASE + fn
        dest = os.path.join(datadir, fn)
        tmp  = dest + ".tmp"

        exp = md5map.get(fn)
        if exp is None:
            print(f"[WARN] {fn}: sin hash esperado en md5_subset.txt; lo descargo igual.")

        for attempt in range(1, 4):
            try:
                code, ctype, blob = download(url, headers)
                if code != 200:
                    raise RuntimeError(f"HTTP {code}")
                if looks_like_html(blob):
                    raise RuntimeError("contenido HTML/JS en vez de texto")
                if len(blob) < 4096:
                    raise RuntimeError(f"demasiado pequeño ({len(blob)} bytes)")

                # MD5 si lo tenemos
                if exp is not None:
                    got = hashlib.md5(blob).hexdigest()
                    if got != exp:
                        raise RuntimeError(f"MD5 mismatch (got {got}, exp {exp})")

                # escribe atomico
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                with open(tmp, "wb") as f:
                    f.write(blob)
                os.replace(tmp, dest)
                ok += 1
                print(f"[OK] {fn}")
                break

            except Exception as e:
                print(f"[RETRY {attempt}/3] {fn}: {e}")
                time.sleep(1.0*attempt)
        else:
            fail += 1
            print(f"[FAIL] {fn}")

    print(f"\nResumen: OK={ok}, FAIL={fail}")
    sys.exit(0 if fail==0 else 2)

if __name__ == "__main__":
    main()
