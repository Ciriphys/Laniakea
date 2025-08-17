import matplotlib.pyplot as plt
import numpy as np

cmap = plt.get_cmap("seismic", 256)
rgb = (cmap(np.linspace(0, 1, 256))[:, :3] * 255).astype(int)

with open("seismic.h", "w") as f:
    f.write("unsigned char seismic[256][3] = {\n")
    for i, (r, g, b) in enumerate(rgb):
        f.write(f"    {{ {r:3d}, {g:3d}, {b:3d} }},")
        if (i + 1) % 8 == 0:  # 8 per riga
            f.write("\n")
    f.write("};\n")