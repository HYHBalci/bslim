"""Make page contact sheets and a text sanity record for manual visual QA."""
import argparse
import json
from pathlib import Path

from PIL import Image, ImageDraw
import pypdfium2 as pdfium
from pypdf import PdfReader

p = argparse.ArgumentParser()
p.add_argument("pdf", type=Path)
p.add_argument("--out", type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
doc = pdfium.PdfDocument(str(a.pdf))
reader = PdfReader(str(a.pdf))
pages = []
for i, page in enumerate(doc):
    image = page.render(scale=1.12).to_pil().convert("RGB")
    path = a.out / f"page-{i+1:03d}.png"
    image.save(path)
    extracted = reader.pages[i].extract_text() or ""
    pages.append({"page": i+1, "words": len(extracted.split()), "replacement_characters": extracted.count("\ufffd"),
                  "width_pt": float(reader.pages[i].mediabox.width), "height_pt": float(reader.pages[i].mediabox.height)})
for first in range(0, len(doc), 6):
    thumb_w, thumb_h = 670, 965
    contact = Image.new("RGB", (2*thumb_w, 3*thumb_h), "#d8dde0")
    draw = ImageDraw.Draw(contact)
    for offset in range(min(6, len(doc)-first)):
        i = first+offset
        image = Image.open(a.out / f"page-{i+1:03d}.png")
        image.thumbnail((thumb_w-16, thumb_h-35))
        x, y = (offset % 2)*thumb_w, (offset // 2)*thumb_h
        draw.text((x+10, y+8), f"Page {i+1}", fill="black")
        contact.paste(image, (x+8, y+28))
    contact.save(a.out / f"contact-{first//6+1:02d}.png")
(a.out / "qa.json").write_text(json.dumps({"pdf": str(a.pdf), "page_count": len(doc), "pages": pages}, indent=2), encoding="utf-8")
print(json.dumps({"pages": len(doc), "contact_sheets": (len(doc)+5)//6, "out": str(a.out)}))
