"""
Download all Virgo logbook search result pages for a keyword.

The logbook requires EGO SSO authentication. You need to provide your
session cookie, obtained from your browser after logging in:

  Browser → Developer Tools → Network → any logbook request
  → Request Headers → Cookie → copy the full value

Usage:
    python fetch_logbook_entries.py --keyword "25-minute" --cookie "PHPSESSID=xxx..."
    python fetch_logbook_entries.py --keyword "25-minute" --cookie-file ~/.virgo_cookie

Output: one HTML file per page in --output-dir (default: data/logbook_25min/)
"""
import argparse
import time
from pathlib import Path

import requests
from bs4 import BeautifulSoup

BASE_URL = "https://logbook.virgo-gw.eu/virgo/"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--keyword", default="25-minute")
    p.add_argument("--cookie", default=None, help="Raw cookie string from browser")
    p.add_argument("--cookie-file", default=None, help="File containing the cookie string")
    p.add_argument("--output-dir", default=None)
    p.add_argument("--delay", type=float, default=1.0, help="Seconds between requests")
    return p.parse_args()


def get_cookie_str(args):
    if args.cookie:
        return args.cookie
    if args.cookie_file:
        return Path(args.cookie_file).read_text().strip()
    raise SystemExit("Provide --cookie or --cookie-file")


def cookie_str_to_dict(cookie_str):
    cookies = {}
    for part in cookie_str.split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            cookies[k.strip()] = v.strip()
    return cookies


def fetch_page(session, keyword, page):
    params = {
        "r": "search",
        "kw": keyword,
        "page": page,
    }
    resp = session.get(BASE_URL, params=params, timeout=30)
    resp.raise_for_status()
    return resp.text


def count_pages(html):
    soup = BeautifulSoup(html, "html.parser")
    # Look for pagination links
    pages = set()
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if "page=" in href:
            try:
                p = int(href.split("page=")[1].split("&")[0])
                pages.add(p)
            except ValueError:
                pass
    return max(pages) if pages else 1


def main():
    args = parse_args()
    cookie_str = get_cookie_str(args)
    cookies = cookie_str_to_dict(cookie_str)

    output_dir = Path(args.output_dir) if args.output_dir else \
        Path(__file__).parent.parent / "data" / f"logbook_{args.keyword.replace(' ', '_')}"
    output_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.cookies.update(cookies)
    session.headers.update({"User-Agent": "Mozilla/5.0 (NOISEHOUND logbook fetcher)"})

    print(f"Fetching page 1 for keyword '{args.keyword}'...")
    html = fetch_page(session, args.keyword, 1)

    # Check authentication
    if "login" in html.lower() and "password" in html.lower():
        raise SystemExit("Cookie rejected — please refresh your browser session cookie.")

    n_pages = count_pages(html)
    print(f"Found {n_pages} page(s) of results.")

    (output_dir / "page_001.html").write_text(html, encoding="utf-8")
    print(f"  Saved page_001.html")

    for page in range(2, n_pages + 1):
        time.sleep(args.delay)
        print(f"Fetching page {page}/{n_pages}...")
        html = fetch_page(session, args.keyword, page)
        fname = output_dir / f"page_{page:03d}.html"
        fname.write_text(html, encoding="utf-8")
        print(f"  Saved {fname.name}")

    print(f"\nDone. {n_pages} page(s) saved to {output_dir}/")


if __name__ == "__main__":
    main()
