"""Capture Figure 2: a screenshot of the live web viewer (3Dmol.js / WebGL).

Generates web/data.js from the pilot results (so no 'sample' banner shows),
serves web/ locally, loads it in headless Chromium, waits for the structure to
render, and writes paper/figures/fig2_viewer.png.
"""
import json
import threading
import time
from functools import partial
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / "web"
OUT = ROOT / "paper" / "figures" / "fig2_viewer.png"
OUT.parent.mkdir(parents=True, exist_ok=True)
PORT = 8137

# 1. Clean data.js from the canonical 4-candidate pilot (matches paper Section 3 /
#    Table 1). Using sample_results.json avoids the 'sample' banner (no flag) and
#    keeps the figure consistent with the manuscript.
results = json.loads((ROOT / "results" / "sample_results.json").read_text(encoding="utf-8"))
(WEB / "data.js").write_text(
    "window.RESCUE_DATA = " + json.dumps(results) + ";\n", encoding="utf-8")

# 2. Serve web/ in a background thread.
handler = partial(SimpleHTTPRequestHandler, directory=str(WEB))
httpd = ThreadingHTTPServer(("127.0.0.1", PORT), handler)
threading.Thread(target=httpd.serve_forever, daemon=True).start()
print(f"serving {WEB} at http://127.0.0.1:{PORT}")

# 3. Screenshot with Playwright.
from playwright.sync_api import sync_playwright

with sync_playwright() as p:
    browser = p.chromium.launch(args=["--use-gl=swiftshader", "--enable-webgl",
                                       "--ignore-gpu-blocklist"])
    page = browser.new_page(viewport={"width": 1280, "height": 820},
                            device_scale_factor=2)
    page.goto(f"http://127.0.0.1:{PORT}/index.html", wait_until="networkidle")
    # Wait for the first row to be selected and the 3D canvas to appear.
    try:
        page.wait_for_selector("#viewer canvas", timeout=15000)
    except Exception:
        print("warning: canvas selector not found within timeout")
    time.sleep(6)  # allow the AlphaFold structure fetch + WebGL render to finish
    page.screenshot(path=str(OUT))
    browser.close()

httpd.shutdown()
print(f"wrote {OUT}")
