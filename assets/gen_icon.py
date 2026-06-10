"""Generate the GenomeForge app icon: a DNA double helix on a rounded blue tile.
Outputs assets/icon.png (256) and assets/icon.ico (multi-size)."""
import math
from PIL import Image, ImageDraw, ImageFilter

SS = 1024  # supersample, downscaled at the end for smooth edges


def lerp(a, b, t):
    return tuple(int(a[i] + (b[i] - a[i]) * t) for i in range(3))


# --- rounded-rect gradient background -------------------------------------
bg = Image.new("RGB", (SS, SS))
top, bot = (22, 28, 44), (12, 46, 92)
px = bg.load()
for y in range(SS):
    col = lerp(top, bot, y / SS)
    for x in range(SS):
        px[x, y] = col

mask = Image.new("L", (SS, SS), 0)
md = ImageDraw.Draw(mask)
m = int(SS * 0.05)
md.rounded_rectangle([m, m, SS - m, SS - m], radius=int(SS * 0.22), fill=255)

tile = Image.new("RGBA", (SS, SS), (0, 0, 0, 0))
tile.paste(bg, (0, 0), mask)

# subtle inner accent border
bd = ImageDraw.Draw(tile)
bd.rounded_rectangle([m, m, SS - m, SS - m], radius=int(SS * 0.22),
                     outline=(90, 150, 220, 160), width=int(SS * 0.012))

# --- DNA double helix ------------------------------------------------------
helix = Image.new("RGBA", (SS, SS), (0, 0, 0, 0))
hd = ImageDraw.Draw(helix)
cx = SS / 2
amp = SS * 0.17
y0, y1 = SS * 0.20, SS * 0.80
turns = 2.1
strand_w = int(SS * 0.035)

STR_A = (130, 210, 235)   # cyan
STR_B = (150, 180, 255)   # periwinkle
BASE = [(96, 200, 120), (220, 120, 120), (235, 190, 90), (120, 160, 235)]


def strand_points(phase):
    pts = []
    n = 220
    for i in range(n + 1):
        t = i / n
        y = y0 + (y1 - y0) * t
        ang = turns * 2 * math.pi * t + phase
        x = cx + amp * math.sin(ang)
        pts.append((x, y))
    return pts


pa = strand_points(0.0)
pb = strand_points(math.pi)

# rungs (base pairs) — behind the strands
rungs = 14
for k in range(rungs):
    t = (k + 0.5) / rungs
    y = y0 + (y1 - y0) * t
    ang = turns * 2 * math.pi * t
    xa = cx + amp * math.sin(ang)
    xb = cx + amp * math.sin(ang + math.pi)
    depth = (math.cos(ang) + 1) / 2  # front/back shading
    col = lerp((70, 90, 120), BASE[k % len(BASE)], depth)
    hd.line([(xa, y), (xb, y)], fill=col + (235,), width=int(SS * 0.022))

# strands on top
hd.line(pa, fill=STR_A + (255,), width=strand_w, joint="curve")
hd.line(pb, fill=STR_B + (255,), width=strand_w, joint="curve")
# nucleotide caps
for (x, y) in pa[::20] + pb[::20]:
    r = SS * 0.022
    hd.ellipse([x - r, y - r, x + r, y + r], fill=(235, 240, 255, 255))

helix = helix.filter(ImageFilter.GaussianBlur(SS * 0.002))
tile = Image.alpha_composite(tile, helix)

# --- export ----------------------------------------------------------------
icon256 = tile.resize((256, 256), Image.LANCZOS)
icon256.save("E:/Repos/genomeforge/assets/icon.png")
tile.resize((256, 256), Image.LANCZOS).save(
    "E:/Repos/genomeforge/assets/icon.ico",
    sizes=[(256, 256), (128, 128), (64, 64), (48, 48), (32, 32), (16, 16)],
)
print("wrote assets/icon.png and assets/icon.ico")
