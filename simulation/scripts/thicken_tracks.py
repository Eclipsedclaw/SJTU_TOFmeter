#!/usr/bin/env python3
"""Widen the particle tracks of a Geant4 offscreen picture, before scripts/picture_finish.py shrinks it.

    python3 thicken_tracks.py in.png out.png [--width 7] [--mask-from tracks.png] "#rrggbb" ...

The offscreen driver (TSG_OFFSCREEN) draws lines one pixel wide, whatever line width is set. Every pixel with the hue
of one of the given track colours, also where a see-through volume in front of the track has blended it, is widened
to a square of --width pixels painted in that colour; the colours are painted in the order given, so the last one
ends on top where tracks cross. --mask-from takes the tracks from another picture of the same view, drawn without
the volumes that hide them (the driver hides lines even behind see-through faces), and paints them over this one. A pixel has a colour's hue when its chroma (the pixel minus its grey level) points
the same way (cosine above 0.95) and is at least 40% as strong: the shaded volumes of the TOF-meter pictures never
pass this for red or blue. scripts/trident_views_lyso.sh runs it. Needs numpy and Pillow.
"""
import argparse

import numpy as np
from PIL import Image, ImageFilter


def main():
    parser = argparse.ArgumentParser(description="Widen the tracks of a Geant4 offscreen picture")
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("colours", nargs="+", help='track colours "#rrggbb", painted in this order')
    parser.add_argument("--width", type=int, default=7, help="width of the widened tracks in pixels (odd)")
    parser.add_argument("--mask-from", help="find the tracks in this picture of the same view instead")
    args = parser.parse_args()
    if args.width < 1 or args.width % 2 == 0:
        parser.error("--width must be an odd number of pixels")

    image = Image.open(args.input).convert("RGB")
    source = Image.open(args.mask_from).convert("RGB") if args.mask_from else image
    if source.size != image.size:
        parser.error(f"{args.mask_from} is {source.size[0]} x {source.size[1]}, {args.input} {image.size[0]} x {image.size[1]}")
    pixels = np.asarray(source, dtype=np.float32)
    chroma = pixels - pixels.mean(axis=2, keepdims=True)
    strength = np.linalg.norm(chroma, axis=2)
    result = np.array(image)
    for colour in args.colours:
        rgb = np.array([int(colour[i:i + 2], 16) for i in (1, 3, 5)], dtype=np.float32)
        reference = rgb - rgb.mean()
        cosine = (chroma @ reference) / np.maximum(strength * np.linalg.norm(reference), 1e-6)
        mask = (cosine > 0.95) & (strength >= 0.4 * np.linalg.norm(reference))
        wide = Image.fromarray(mask.astype(np.uint8) * 255).filter(ImageFilter.MaxFilter(args.width))
        result[np.asarray(wide) > 0] = rgb.astype(np.uint8)
        print(f"{colour}: {int(mask.sum())} track pixels")
    Image.fromarray(result).save(args.output)


if __name__ == "__main__":
    main()
