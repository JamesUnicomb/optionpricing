import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import parse

c = []
h = []
k = []

n = 0
with open("output.txt", 'r') as f:
    for l in f.readlines():
        l = l.rstrip()

        s,level,ht,kt,_,_,_,_,_= parse.parse("set {:d} level = {:02d}, h = {:02d}, k = {:02d}, idx = {:02d}, i = {:02d}, j = {:02d}, idx = {:02d}, blocksize = {:02d}", l)

        if s == 2 or s == 1:
            h.append(ht)
            k.append(kt)
            c.append(s)
            n += 1

print(h)
print(k)

plt.figure(figsize=(6,6))
plt.scatter([kt for kt, m in zip(k,c) if m == 1], [ht for ht, m in zip(h,c) if m == 1])
plt.scatter([kt for kt, m in zip(k,c) if m == 2], [ht for ht, m in zip(h,c) if m == 2], marker="x")
plt.xlim(-1,49)
plt.ylim(-1,49)
plt.savefig("output.png")
plt.close()


fig = plt.figure(figsize=(6,6))
s = plt.scatter([],[])
plt.xlim(-1,49)
plt.ylim(-1,49)

def f(i):
    x = k[:i]
    y = h[:i]

    s.set_offsets(np.column_stack([x,y]))

anim = FuncAnimation(fig, f, frames=n, interval=25)
anim.save("plot.gif", writer='imagemagick', fps=25)
