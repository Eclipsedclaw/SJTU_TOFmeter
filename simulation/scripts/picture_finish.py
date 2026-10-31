#!/usr/bin/env python3
"""Finish a picture of the Geant4 offscreen driver (TSG_OFFSCREEN) for slides.

    python3 picture_finish.py in.png out.png [--scale 0.5] [--legend right] "label=#rrggbb" ...

Cuts the white border (keeping a margin), shrinks the picture (the driver draws without smoothing, so a large
render shrunk by half has smooth edges) and adds a legend of the colours: in a white band added on the right or
the left (--legend right, left), so that it hides nothing, or over a corner (tl, tr, bl, br). An entry
"label=#rrggbb/o" gets an outline swatch, for parts drawn as outlines. Needs Pillow and the DejaVu Sans font.
scripts/geometry_czt.sh runs it on the pictures of macro/vis_czt.mac.
"""
import argparse

from PIL import Image, ImageChops, ImageDraw, ImageFont

FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
INK = (34, 34, 34)


def main():
    parser = argparse.ArgumentParser(description="Finish a Geant4 offscreen picture for slides")
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("entries", nargs="*", help='legend entries "label=#rrggbb", "label=#rrggbb/o" for an outline')
    parser.add_argument("--scale", type=float, default=0.5, help="size of the output relative to the input")
    parser.add_argument("--legend", choices=["right", "left", "tl", "tr", "bl", "br"], default="right",
                        help="where the legend goes")
    parser.add_argument("--margin", type=float, default=0.04, help="white margin kept, fraction of the content")
    args = parser.parse_args()

    image = Image.open(args.input).convert("RGB")
    box = ImageChops.difference(image, Image.new("RGB", image.size, (255, 255, 255))).getbbox()
    if box:
        pad = int(args.margin * max(box[2] - box[0], box[3] - box[1]))
        image = image.crop((max(box[0] - pad, 0), max(box[1] - pad, 0),
                            min(box[2] + pad, image.width), min(box[3] + pad, image.height)))
    image = image.resize((round(image.width * args.scale), round(image.height * args.scale)),
                         Image.Resampling.LANCZOS)

    if args.entries:
        size = max(14, round(image.height / 40))
        font = ImageFont.truetype(FONT, size)
        rows = []
        for entry in args.entries:
            label, _, colour = entry.rpartition("=")
            outline = colour.endswith("/o")
            rows.append((label, colour.removesuffix("/o"), outline))
        step, gap = round(1.6 * size), round(0.6 * size)
        text_width = max(ImageDraw.Draw(image).textlength(label, font=font) for label, _, _ in rows)
        width = round(size + gap + text_width + 2 * gap)
        height = len(rows) * step + gap
        edge = round(0.03 * image.height)
        if args.legend in ("right", "left"):
            band = Image.new("RGB", (image.width + width + 2 * edge, image.height), (255, 255, 255))
            band.paste(image, (0 if args.legend == "right" else width + 2 * edge, 0))
            image = band
            x0 = image.width - edge - width if args.legend == "right" else edge
            y0 = edge
        else:
            x0 = edge if args.legend[1] == "l" else image.width - edge - width
            y0 = edge if args.legend[0] == "t" else image.height - edge - height
        draw = ImageDraw.Draw(image)
        draw.rectangle((x0, y0, x0 + width, y0 + height), fill=(255, 255, 255), outline=(200, 200, 200))
        for k, (label, colour, outline) in enumerate(rows):
            top = y0 + gap + k * step
            swatch = (x0 + gap, top, x0 + gap + size, top + size)
            if outline:
                draw.rectangle(swatch, outline=colour, width=max(2, size // 8))
            else:
                draw.rectangle(swatch, fill=colour, outline=(90, 90, 90))
            draw.text((x0 + 2 * gap + size, top + size / 2), label, font=font, fill=INK, anchor="lm")
    image.save(args.output)
    print(f"{args.output}: {image.width} x {image.height}")


if __name__ == "__main__":
    main()
