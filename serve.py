#!/usr/bin/env python3
"""Serve this folder for the reviewer page, with HTTP range requests (which igv.js needs for BAM and FASTA files).

    python3 serve.py [port]        # default 8000; then open http://localhost:8000/

Python's own http.server ignores Range headers, so the read panels would not load from it; this one answers them.
Standard library only. Serves the folder it is in.
"""
import os
import re
import sys
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


class RangeHandler(SimpleHTTPRequestHandler):
    extensions_map = dict(SimpleHTTPRequestHandler.extensions_map, **{
        '.js': 'text/javascript', '.css': 'text/css', '.html': 'text/html', '.tsv': 'text/tab-separated-values',
        '.bed': 'text/plain', '.fa': 'text/plain', '.fai': 'text/plain', '.json': 'application/json',
        '.bam': 'application/octet-stream', '.bai': 'application/octet-stream', '.md': 'text/markdown', '': 'application/octet-stream'})

    def send_head(self):
        path = self.translate_path(self.path)
        if os.path.isdir(path) or not os.path.exists(path):
            return super().send_head()
        rng = self.headers.get('Range')
        size = os.path.getsize(path)
        m = re.match(r'bytes=(\d*)-(\d*)$', rng or '')
        if not m:
            self.send_response(200)
            self.send_header('Content-type', self.guess_type(path))
            self.send_header('Content-Length', str(size))
            self.send_header('Accept-Ranges', 'bytes')
            self.end_headers()
            return open(path, 'rb')
        start = int(m.group(1)) if m.group(1) else max(0, size - int(m.group(2)))
        end = int(m.group(2)) if m.group(1) and m.group(2) else size - 1
        end = min(end, size - 1)
        if start > end or start >= size:
            self.send_response(416)
            self.send_header('Content-Range', f'bytes */{size}')
            self.end_headers()
            return None
        self.send_response(206)
        self.send_header('Content-type', self.guess_type(path))
        self.send_header('Content-Range', f'bytes {start}-{end}/{size}')
        self.send_header('Content-Length', str(end - start + 1))
        self.send_header('Accept-Ranges', 'bytes')
        self.end_headers()
        f = open(path, 'rb')
        f.seek(start)
        return _Slice(f, end - start + 1)

    def log_message(self, fmt, *args):
        pass


class _Slice:
    def __init__(self, f, n):
        self.f, self.n = f, n

    def read(self, k=-1):
        if self.n <= 0:
            return b''
        k = self.n if k < 0 else min(k, self.n)
        b = self.f.read(k)
        self.n -= len(b)
        return b

    def close(self):
        self.f.close()


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    srv = ThreadingHTTPServer(('127.0.0.1', port), RangeHandler)
    print(f'serving {os.getcwd()} at http://localhost:{port}/  (Ctrl-C stops)')
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
