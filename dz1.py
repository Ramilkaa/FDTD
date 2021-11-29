import pylab
import numpy as np
import json
import pathlib

def f(x):
    return 100 * (np.sqrt(abs(1 - 0.01 * x ** 2))) + 0.01 *(abs(x + 10))
x_min = -15
x_max = 5
dx = 0.1

x = np.arange(x_min, x_max, dx)
y = f(x)

res = {
    "x": x.tolist(),
    "y": y.tolist(),
}

path = pathlib.Path("results")
path.mkdir(exist_ok=True)
file = path / "result_task1.json"

pylab.plot(x, y)
pylab.grid()
pylab.savefig("results/task1.png")
