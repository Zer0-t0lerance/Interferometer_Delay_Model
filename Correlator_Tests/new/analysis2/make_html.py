#!/usr/bin/env python
# make_html.py - собирает OTCHET.md в один самодостаточный OTCHET.html
# с картинками, вшитыми прямо в файл (base64). Файл можно открыть в любом браузере
# и переслать одним куском, без папки figures.
import base64, io, os, re
import markdown

H = os.path.dirname(os.path.abspath(__file__))
src = io.open(os.path.join(H, "OTCHET.md"), encoding="utf-8").read()

def inline(m):
    alt, path = m.group(1), m.group(2)
    full = os.path.join(H, path.replace("/", os.sep))
    if not os.path.exists(full):
        return m.group(0)
    b64 = base64.b64encode(open(full, "rb").read()).decode("ascii")
    return "![%s](data:image/png;base64,%s)" % (alt, b64)

src = re.sub(r'!\[([^\]]*)\]\(([^)]+\.png)\)', inline, src)
body = markdown.markdown(src, extensions=["tables", "fenced_code"])

CSS = """
body{font-family:'Segoe UI',Arial,sans-serif;max-width:1100px;margin:0 auto;padding:30px 40px;
     line-height:1.55;color:#1c1c1c;background:#fff}
h1{border-bottom:3px solid #24506b;padding-bottom:10px}
h2{margin-top:2.2em;border-bottom:1px solid #ccc;padding-bottom:5px;color:#24506b}
h3{margin-top:1.6em;color:#33607d}
img{max-width:100%;display:block;margin:18px auto;border:1px solid #ddd;border-radius:4px}
table{border-collapse:collapse;margin:16px 0;font-size:.93em}
th,td{border:1px solid #bbb;padding:5px 10px;text-align:left}
th{background:#eef3f7}
tr:nth-child(even) td{background:#fafafa}
code,pre{font-family:Consolas,monospace;background:#f4f4f4}
pre{padding:12px;border-left:4px solid #24506b;overflow-x:auto}
blockquote{border-left:4px solid #d9a441;background:#fdf8ee;margin:16px 0;padding:10px 18px}
strong{color:#0b3d5c}
"""

html = ("<!DOCTYPE html><html lang='ru'><head><meta charset='utf-8'>"
        "<title>Сравнение моделей задержки</title><style>%s</style></head><body>%s</body></html>"
        % (CSS, body))
out = os.path.join(H, "OTCHET.html")
io.open(out, "w", encoding="utf-8").write(html)
print("->", out, "%.1f МБ" % (os.path.getsize(out) / 1048576.))
