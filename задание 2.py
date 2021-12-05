#!/bin/env python3
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import requests
from bs4 import BeautifulSoup
from scipy.constants import pi, speed_of_light
from scipy.special import spherical_jn as jn
from scipy.special import spherical_yn as yn


def hn(n, x):
    return jn(n, x) + 1j * yn(n, x)


def bn(n, x):
    return (x * jn(n - 1, x) - n * jn(n, x)) / (x * hn(n - 1, x) - n * hn(n, x))


def an(n, x):
    return jn(n, x) / hn(n, x)


link = "https://jenyay.net/uploads/Student/Modelling/task_02.xml"

result = requests.get(link)
if result.status_code == requests.codes.OK:
    file = result.text
else:
    exit(1)

soup = BeautifulSoup(file, "xml")
var = soup.find(attrs={"number": "3"})
D = float(var["D"])
fmin = float(var["fmin"])
fmax = float(var["fmax"])
step = 1e6

r = D / 2
freq = np.arange(fmin, fmax, step)
lambd = speed_of_light / freq
k = 2 * pi / lambd

arr_sum = [
    ((-1) ** n) * (n + 0.5) * (bn(n, k * r) - an(n, k * r)) for n in range(1, 50)
]
sum_ = np.sum(arr_sum, axis=0)
rcs = (lambd ** 2) / pi * (np.abs(sum_) ** 2)

path = pathlib.Path("results")
path.mkdir(exist_ok=True)

xmlfale="<?xml version=\"1.1\" encoding=\"UTF-8\" ?>\n<data>\n"
xmlfale+="    <frequencydata>\n"
for i in freq:
    xmlfale+=("        <f>"+str(i)+"</f>\n")
xmlfale+="    </frequencydata>\n"
xmlfale+="    <lambdadata>\n"
for i in lambd:
    xmlfale+=("        <lambdada>"+str(i)+"</lambdada>\n")
xmlfale+="    </lambdadata>\n"
xmlfale+="    <rcsdata>\n"
for i in rcs:
    xmlfale+=("        <rcs>"+str(i)+"</rcs>\n")
xmlfale+="    </rcsdata>\n"
xmlfale+="</data>"


file = path / "result_task2.xml"
out = file.open("w")
out.write(xmlfale)
out.close()


plt.plot(freq / 10e6, rcs)
plt.xlabel("$f, МГц$")
plt.ylabel(r"$\sigma, м^2$")
plt.grid()
plt.show()
plt.savefig("results/task2.png")
