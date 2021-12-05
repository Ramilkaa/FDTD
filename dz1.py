import pylab
import numpy as np
import json
import pathlib

def f(x):
    return 100 * (np.sqrt(abs(1 - 0.01 * x ** 2))) + 0.01 *(abs(x + 10))
x_min = -15
x_max = 5
dx = 0.05

x = np.arange(x_min, x_max + 0.05, dx)
y = f(x)

res = {
    "x": x.tolist(),
    "y": y.tolist(),
}

path = pathlib.Path("results")
path.mkdir(exist_ok=True)
file = path / "result_task1.json"
out = file.open("w")
json.dump(res, out, indent=4)
out.close()

pylab.plot(x, y)
pylab.grid()
pylab.show()
pylab.savefig("results/task1.png")
