"""Real elevation surfaces: Terrarium terrain tiles (zoom 10) at the 16 locations of PLAN.md.

Each tile (256 x 256) is cut into its 16 non-overlapping 61 x 61 blocks (offsets 6 + 61 i). A block is excluded
only if some pixel is below -50 m. Elevation = R * 256 + G + B / 256 - 32768 (Terrarium encoding).
Writes data/elevation_blocks.npz (blocks, location names, block positions, tile ids).
"""
import math
import os
import urllib.request

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
TILES = os.path.join(DATA, "tiles")
Z = 10
LOCATIONS = [
    ("Luxembourg Ardennes", 50.05, 6.03), ("Swiss Alps", 46.02, 7.75), ("Flevoland", 52.55, 5.60),
    ("Iowa farmland", 42.00, -93.50), ("Grand Canyon", 36.10, -112.10), ("Everest", 27.99, 86.93),
    ("Erg Chebbi dunes", 31.15, -3.97), ("Pennsylvania ridges", 40.60, -77.60), ("Cairngorms", 57.08, -3.67),
    ("Ethiopian highlands", 9.00, 38.70), ("Amazon lowland", -3.50, -62.00), ("Tibetan plateau", 33.00, 88.00),
    ("Black Forest", 48.00, 8.00), ("Mount Fuji", 35.36, 138.73), ("Loess Plateau", 36.50, 109.50),
    ("Central Australia", -25.00, 131.00),
]


def tile_of(lat, lon, z):
    n = 2 ** z
    x = int((lon + 180.0) / 360.0 * n)
    y = int((1.0 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2.0 * n)
    return x, y


blocks, names, pos, tiles = [], [], [], []
for name, lat, lon in LOCATIONS:
    x, y = tile_of(lat, lon, Z)
    path = os.path.join(TILES, f"{Z}_{x}_{y}.png")
    if not os.path.exists(path):
        urllib.request.urlretrieve(f"https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{Z}/{x}/{y}.png", path)
    rgb = np.asarray(Image.open(path).convert("RGB"), dtype=np.float64)
    elev = rgb[..., 0] * 256 + rgb[..., 1] + rgb[..., 2] / 256 - 32768
    assert elev.shape == (256, 256)
    for bi in range(4):
        for bj in range(4):
            b = elev[6 + 61 * bi: 67 + 61 * bi, 6 + 61 * bj: 67 + 61 * bj]
            if b.min() < -50:
                continue
            blocks.append(b)
            names.append(name)
            pos.append((bi, bj))
            tiles.append(f"{Z}/{x}/{y}")
    print(f"{name:22s} tile {Z}/{x}/{y}: range {elev.min():8.1f} .. {elev.max():8.1f} m, "
          f"{sum(n == name for n in names)} blocks kept")
np.savez_compressed(os.path.join(DATA, "elevation_blocks.npz"), blocks=np.array(blocks), names=np.array(names),
                    pos=np.array(pos), tiles=np.array(tiles))
print(len(blocks), "blocks")
