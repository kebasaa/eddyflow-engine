#!/usr/bin/env python3
"""Serve a local folder the way Google Drive and Dropbox serve a shared link.

The engine reads raw data, metadata and biomet from shared links (see
src/src_common/remote_source.f90). To test that without the network, this
server answers the same requests the real providers do, over a local folder,
and the engine is pointed at it with EDDYFLOW_REMOTE_BASE=http://127.0.0.1:<port>.

The shapes below are what the providers were found to send on 2026-09-25:

Google Drive
    GET /embeddedfolderview?id=<id>     HTML, one `flip-entry` per item, each
                                        linking /drive/folders/<id> or
                                        /file/d/<id>/view
    GET /download?id=<id>&...           the file's bytes

Dropbox
    GET  /scl/fo/<key>/<hash>[/<sub>]?rlkey=..   sets the CSRF cookie `t`
    POST /list_shared_link_folder_entries        JSON {"entries": [...],
         form: t, link_key, secure_hash,          "has_more_entries",
         sub_path (no leading slash), rlkey,      "next_request_voucher"}
         and voucher for the next page
    GET  <an entry's href with dl=1>             the file's bytes
    Every item has its own hash, and a folder is only listed with ITS hash -
    the real service answers 404 to the parent's hash with a sub path, and so
    does this one, so the engine is held to what the service accepts.

Failures, to test that the engine survives them:
    --html NAME     serve an HTML "quota exceeded" page instead of NAME
    --missing NAME  answer 404 for NAME
    --page N        Dropbox entries per page, to exercise paging (default 1000)
    --lifetime S    stop serving after S seconds (default 1800), so a server
                    left behind by a killed test cannot outlive it for long
    --log FILE      append the path of every file served whole, one per line,
                    so a test can check that nothing was downloaded twice

Usage:
    remote_server.py --root DIR --port-file FILE [--port 0] [...]
Prints nothing; writes the port it bound to into FILE, then serves until
killed. The links to use are written next to it, in FILE.links.
"""

import argparse
import hashlib
import html
import json
import os
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, quote, unquote, urlparse

TOKEN = "TESTTOKEN"
KEY = "fixturekey"
RLKEY = "fixturerlkey"
log_lock = threading.Lock()


def item_id(rel):
    """A stable id for a path relative to the root; '' is the root."""
    return "g" + hashlib.sha1(("id:" + rel).encode()).hexdigest()[:24]


def item_hash(rel):
    return "A" + hashlib.sha1(("hash:" + rel).encode()).hexdigest()[:22]


class Tree:
    def __init__(self, root):
        self.root = Path(root)
        self.by_id = {}
        for dirpath, dirnames, filenames in os.walk(self.root):
            for name in dirnames + filenames:
                rel = Path(dirpath, name).relative_to(self.root).as_posix()
                self.by_id[item_id(rel)] = rel
        self.by_id[item_id("")] = ""

    def path(self, rel):
        return self.root / rel if rel else self.root

    def children(self, rel):
        p = self.path(rel)
        out = []
        for child in sorted(p.iterdir(), key=lambda c: c.name):
            crel = child.relative_to(self.root).as_posix()
            out.append((child.name, crel, child.is_dir()))
        return out


def make_handler(tree, opts):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def send(self, code, body, ctype="text/html; charset=utf-8", headers=()):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            for k, v in headers:
                self.send_header(k, v)
            self.end_headers()
            self.wfile.write(body)

        def not_found(self):
            self.send(404, b"<!DOCTYPE html><html><title>404</title></html>")

        def file_bytes(self, rel):
            name = rel.rsplit("/", 1)[-1]
            if name in opts.missing:
                return self.not_found()
            if name in opts.html:
                return self.send(200, b"<!DOCTYPE html><html><head><title>"
                                 b"Google Drive - Quota exceeded</title></head></html>")
            p = tree.path(rel)
            if not p.is_file():
                return self.not_found()
            self.send(200, p.read_bytes(), "application/octet-stream")
            if opts.log:
                with log_lock, open(opts.log, "a", encoding="utf-8") as f:
                    f.write(rel + "\n")

        # --- Google Drive -------------------------------------------------
        def gdrive_listing(self, rel):
            p = tree.path(rel)
            if not p.is_dir():
                return self.not_found()
            entries = []
            for name, crel, is_dir in tree.children(rel):
                cid = item_id(crel)
                link = (f"https://drive.google.com/drive/folders/{cid}" if is_dir
                        else f"https://drive.google.com/file/d/{cid}/view?usp=drive_web")
                entries.append(
                    f'<div class="flip-entry" id="entry-{cid}" tabindex="0" role="link">'
                    f'<div class="flip-entry-info"><a href="{link}" target="_blank">'
                    f'<div class="flip-entry-title">{html.escape(name)}</div></a></div>'
                    f'<div class="flip-entry-last-modified"><div>3:17 am</div></div></div>')
            body = ('<!DOCTYPE html><html><head><title>' + html.escape(p.name) +
                    '</title></head><body><div class="flip-entries">' +
                    "".join(entries) + "</div></body></html>")
            self.send(200, body.encode("utf-8"))

        # --- Dropbox ------------------------------------------------------
        def dbx_href(self, rel):
            return (f"https://www.dropbox.com/scl/fo/{KEY}/{item_hash(rel)}/"
                    f"{quote(rel)}?rlkey={RLKEY}&dl=0")

        def do_GET(self):
            url = urlparse(self.path)
            q = parse_qs(url.query)
            if url.path == "/embeddedfolderview":
                rel = tree.by_id.get(q.get("id", [""])[0])
                return self.not_found() if rel is None else self.gdrive_listing(rel)
            if url.path == "/download":
                rel = tree.by_id.get(q.get("id", [""])[0])
                return self.not_found() if rel is None else self.file_bytes(rel)
            if url.path.startswith(f"/scl/fo/{KEY}/"):
                rest = url.path[len(f"/scl/fo/{KEY}/"):]
                hsh, _, sub = rest.partition("/")
                sub = unquote(sub)
                if hsh != item_hash(sub) or q.get("rlkey", [""])[0] != RLKEY:
                    return self.not_found()
                if q.get("dl", ["0"])[0] == "1" and tree.path(sub).is_file():
                    return self.file_bytes(sub)
                return self.send(200, b"<!DOCTYPE html><html><body>app</body></html>",
                                 headers=[("Set-Cookie", f"t={TOKEN}; Path=/")])
            return self.not_found()

        def do_POST(self):
            url = urlparse(self.path)
            if url.path != "/list_shared_link_folder_entries":
                return self.not_found()
            n = int(self.headers.get("Content-Length", "0"))
            form = {k: v[0] for k, v in parse_qs(self.rfile.read(n).decode(),
                                                    keep_blank_values=True).items()}
            cookie = self.headers.get("Cookie", "")
            if (form.get("t") != TOKEN or f"t={TOKEN}" not in cookie
                    or form.get("link_key") != KEY or form.get("rlkey") != RLKEY):
                return self.send(403, b'{"error": "forbidden"}', "application/json")
            sub = form.get("sub_path", "")
            if sub.startswith("/") or form.get("secure_hash") != item_hash(sub) \
                    or not tree.path(sub).is_dir():
                return self.not_found()
            kids = tree.children(sub)
            start = int(form.get("voucher") or 0)
            page = kids[start:start + opts.page]
            more = start + opts.page < len(kids)
            entries = []
            for name, crel, is_dir in page:
                e = {"_mount_access_perms": ["can_view"], "filename": name,
                     "href": self.dbx_href(crel), "is_dir": is_dir,
                     "sort_key": ["x"], "open_in_app": None}
                if not is_dir:
                    e["bytes"] = tree.path(crel).stat().st_size
                    e["ts"] = 1788170662
                entries.append(e)
            body = {"entries": entries, "share_tokens": [], "has_more_entries": more,
                    "next_request_voucher": str(start + opts.page) if more else None}
            self.send(200, json.dumps(body).encode(), "application/json")

    return Handler


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--port", type=int, default=0)
    ap.add_argument("--port-file", required=True)
    ap.add_argument("--html", action="append", default=[])
    ap.add_argument("--missing", action="append", default=[])
    ap.add_argument("--page", type=int, default=1000)
    ap.add_argument("--lifetime", type=float, default=1800)
    ap.add_argument("--log", default="")
    opts = ap.parse_args()

    tree = Tree(opts.root)
    server = ThreadingHTTPServer(("127.0.0.1", opts.port), make_handler(tree, opts))
    port = server.server_address[1]

    links = {
        "gdrive_root": f"https://drive.google.com/drive/folders/{item_id('')}?usp=sharing",
        "dropbox_root": (f"https://www.dropbox.com/scl/fo/{KEY}/{item_hash('')}"
                         f"?rlkey={RLKEY}&st=abc&dl=0"),
        "gdrive_files": {rel: f"https://drive.google.com/file/d/{i}/view?usp=sharing"
                         for i, rel in tree.by_id.items() if rel and tree.path(rel).is_file()},
        "gdrive_folders": {rel: f"https://drive.google.com/drive/folders/{i}"
                           for i, rel in tree.by_id.items() if rel and tree.path(rel).is_dir()},
        "dropbox_files": {rel: (f"https://www.dropbox.com/scl/fo/{KEY}/{item_hash(rel)}/"
                                f"{quote(rel)}?rlkey={RLKEY}&dl=0")
                          for rel in tree.by_id.values() if rel and tree.path(rel).is_file()},
    }
    Path(opts.port_file + ".links").write_text(json.dumps(links, indent=1))
    Path(opts.port_file).write_text(str(port))
    threading.Timer(opts.lifetime, server.shutdown).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    sys.exit(main())
