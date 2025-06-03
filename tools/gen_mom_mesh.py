import numpy as np

mesh_number = 100
p_min = 0.0
p_max = 1200.0
mesh_file = "table_momentum_mesh.txt"


def get_mesh_ori(precision, deg):
    np.set_printoptions(precision)
    temp = np.polynomial.legendre.leggauss(deg)
    xlist = []
    wlist = []
    for i in range(0, len(temp[0])):
        xlist.append(temp[0][i])
        wlist.append(temp[1][i])
    return xlist, wlist


def map(x, w, a, b):
    xx = []
    ww = []
    for _, xi in enumerate(x):
        wi = w[_]
        xxi = (b - a) / 2.0 * xi + (b + a) / 2.0
        wwi = (b - a) / 2.0 * wi
        xx.append(xxi)
        ww.append(wwi)
    return xx, ww


xlist, wlist = get_mesh_ori(precision=16, deg=mesh_number)


mesh_mom, weight_mom = map(xlist, wlist, p_min, p_max)


print("mom mesh :\n")
print(mesh_mom, sep=",", end="\n\n", flush=False)
print(weight_mom, sep=",", end="\n\n", flush=False)


with open(mesh_file, "w") as f:
    f.write("# momentum mesh points number;\n")
    f.write("# momentum mesh points and weights;\n")
    f.write(f"{mesh_number}\n")
    for pi, wi in zip(mesh_mom, weight_mom):
        f.write(f"{pi:.17e} {wi:.17e}\n")

print(f"file written in : {mesh_file}")
