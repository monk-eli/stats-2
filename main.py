import csv
import math

from math import *
import numpy as np


def get_sample(file, line, sample_size):
    with open(file) as csvfile:
        reader = csv.reader(csvfile)
        listedData = list(reader)

    return [float(i) for i in listedData[line - 1][:sample_size]]


def count_freq(sample, a, b):
    count = 0
    for i in sample:
        if a < i < b:
            count += 1

    return count


def count_nij(sample, a, b):  # for continuous distribution
    result = 0
    for i in sample:
        if a < i < b:
            result += 1
    return result


def count_nij_disc(sample, a):  # for discrete distribution
    result = 0
    for i in sample:
        if i == a:
            result += 1
    return result


def binomial(n, p, xi):
    return math.comb(n, xi) * (p ** xi) * ((1 - p) ** (n - xi))


print("Task 1 - 19 (Normal)")

sample1 = get_sample("task1.csv", 19, 30)

freq = [count_freq(sample1, -inf, -2), count_freq(sample1, -2, -1), count_freq(sample1, -1, 0),
        count_freq(sample1, 0, 1), count_freq(sample1, 1, 2), count_freq(sample1, 2, inf)]
p = [0.0228, 0.1359, 0.3413, 0.3413, 0.1359, 0.0228]
n_p = [i * 30 for i in p]
xi = [(freq[i] - n_p[i]) ** 2 / n_p[i] for i in range(6)]

print("Frequency:", freq)
print("P:", p)
print("nP:", n_p)
print("xi:", xi)
print("xi observed:", np.sum(xi), '\n')

print("Task 1 - 18 (Uniform)")

sample2 = get_sample("task1.csv", 18, 30)

freq = [count_freq(sample2, -2, -1), count_freq(sample2, -1, 0),
        count_freq(sample2, 0, 1), count_freq(sample2, 1, 2)]
n_p = [7.5] * 4
xi = [(freq[i] - n_p[i]) ** 2 / n_p[i] for i in range(4)]

print("Frequency:", freq)
print("P:", n_p)
print("xi:", xi)
print("xi observed:", np.sum(xi), '\n')

print("Task 1 - 11 (Uniform discrete)")

sample3 = get_sample("task1.csv", 11, 30)

freq = [count_freq(sample3, 0, 2), count_freq(sample3, 1, 3),
        count_freq(sample3, 2, 4), count_freq(sample3, 3, 5)]
n_p = [7.5] * 4
xi = [(freq[i] - n_p[i]) ** 2 / n_p[i] for i in range(4)]

print("Frequency:", freq)
print("P:", n_p)
print("xi:", xi)
print("xi observed:", np.sum(xi), '\n')

print("Task 1 - 24 (Exponential)")

sample4 = get_sample("task1.csv", 24, 30)

freq = [count_freq(sample4, 0, 1), count_freq(sample4, 1, 2),
        count_freq(sample4, 2, inf)]
n_p = [18.964, 6.975, 4.06]
xi = [(freq[i] - n_p[i]) ** 2 / n_p[i] for i in range(3)]

print("Frequency:", freq)
print("P:", p)
print("nP:", n_p)
print("xi:", xi)
print("xi observed:", np.sum(xi), '\n')

print("Task 1 - 22 (Bin)")

sample5 = get_sample("task1.csv", 22, 30)

freq = [count_freq(sample5, -1, 1), count_freq(sample5, 0, 2), count_freq(sample5, 1, 3),
        count_freq(sample5, 2, 4)]
p = [0.173, 0.412, 0.328, 0.261]
n_p = [30 * i for i in p]
xi = [(freq[i] - n_p[i]) ** 2 / n_p[i] for i in range(4)]
print("xi observed:", np.sum(xi), '\n')

print("Frequency:", freq)
print("nP:", n_p)
print("Mean value", np.mean(sample5))
print("xi:", xi)
print("xi observed:", np.sum(xi), '\n')

print("Task 2 - 6 (Bin)")

sample1 = get_sample("task2.csv", 27, 30)
sample2 = get_sample("task2.csv", 28, 30)
sample3 = get_sample("task2.csv", 29, 30)

nij = [[0] * 4, [0] * 4, [0] * 4]
vj = [0] * 4

for j in range(4):
    i = 0
    for sample in sample1, sample2, sample3:
        vj[j] += count_nij_disc(sample, j)
        nij[i][j] = count_nij_disc(sample, j)
        i += 1

xi = 0
for i in range(3):
    for j in range(4):
        xi += (nij[i][j] - vj[j] * 30 / 90) ** 2 / (vj[j] * 30)

xi *= 90

print("vj:", vj)
print("nij:", nij)
print("xi:", xi, '\n')

sample_w = sample1 + sample2 + sample3
p = np.mean(sample_w) / 3

n_p = [0] * 4
for i in range(4):
    n_p[i] = 90 * binomial(3, p, i)

xi = [0] * 4
for i in range(4):
    xi[i] = (vj[i] - n_p[i]) ** 2 / n_p[i]

print("p:", p)
print("n_p:", n_p)
print("xi", xi)
print("xi observed:", np.sum(xi))

print("Task 2 - 21 (Exp)")

sample1 = get_sample("task2.csv", 102, 30)
sample2 = get_sample("task2.csv", 103, 30)
sample3 = get_sample("task2.csv", 104, 30)

nij = [[0] * 3, [0] * 3, [0] * 3]
vj = [0] * 3

for j in range(3):
    if j == 0:
        left, right = 0, 1
    elif j == 1:
        left, right = 1, 2
    elif j == 2:
        left, right = 2, inf

    i = 0
    for smp in sample1, sample2, sample3:
        vj[j] += count_nij(smp, left, right)
        nij[i][j] = count_nij(smp, left, right)
        i += 1

xi = 0
for i in range(3):
    for j in range(3):
        xi += (nij[i][j] - vj[j] * 30 / 90) ** 2 / (vj[j] * 30)

xi *= 90

print("vj:", vj)
print("nij:", nij)
print("xi:", xi, '\n')

sample_w = sample1 + sample2 + sample3
# p = np.mean(sample_w) / 3

n_p = [0.6321 * 90, 0.2325 * 90, 0.1353 * 90]

xi = [0] * 3
for i in range(3):
    xi[i] = (vj[i] - n_p[i]) ** 2 / n_p[i]

# print("p:", p)
print("n_p:", n_p)
print("xi", xi)
print("xi observed:", np.sum(xi))

print("Task 2 - 19 (Normal)")

sample1 = get_sample("task2.csv", 92, 30)
sample2 = get_sample("task2.csv", 93, 30)
sample3 = get_sample("task2.csv", 94, 30)

nij = [[0] * 6, [0] * 6, [0] * 6]
vj = [0] * 6

for j in range(6):
    if j == 0:
        left, right = -inf, -2
    elif j == 1:
        left, right = -2, -1
    elif j == 2:
        left, right = -1, 0
    elif j == 3:
        left, right = 0, 1
    elif j == 4:
        left, right = 1, 2
    elif j == 5:
        left, right = 2, inf

    i = 0
    for smp in sample1, sample2, sample3:
        vj[j] += count_nij(smp, left, right)
        nij[i][j] = count_nij(smp, left, right)
        i += 1

xi = 0
for i in range(3):
    for j in range(6):
        xi += (nij[i][j] - vj[j] * 30 / 90) ** 2 / (vj[j] * 30)

xi *= 90

print("vj:", vj)
print("nij:", nij)
print("xi:", xi, '\n')

sample_w = sample1 + sample2 + sample3
# p = np.mean(sample_w) / 3

n_p = [0.0228 * 90, 0.1359 * 90, 0.3413 * 90, 0.3413 * 90, 0.1359 * 90, 0.0228 * 90]

xi = [0] * 6
for i in range(6):
    xi[i] = (vj[i] - n_p[i]) ** 2 / n_p[i]

# print("p:", p)
print("n_p:", n_p)
print("xi", xi)
print("xi observed:", np.sum(xi))

print("Task 2 - 16 (Uniform discrete)")

sample1 = get_sample("task2.csv", 77, 30)
sample2 = get_sample("task2.csv", 78, 30)
sample3 = get_sample("task2.csv", 79, 30)

nij = [[0] * 4, [0] * 4, [0] * 4]
vj = [0] * 4

for j in range(4):
    if j == 0:
        t = 1
    elif j == 1:
        t = 2
    elif j == 2:
        t = 3
    elif j == 3:
        t = 4

    i = 0
    for smp in sample1, sample2, sample3:
        vj[j] += count_nij_disc(smp, t)
        nij[i][j] = count_nij_disc(smp, t)
        i += 1

xi = 0
for i in range(3):
    for j in range(4):
        xi += (nij[i][j] - vj[j] * 30 / 90) ** 2 / (vj[j] * 30)

xi *= 90

print("vj:", vj)
print("nij:", nij)
print("xi:", xi, '\n')

sample_w = sample1 + sample2 + sample3
# p = np.mean(sample_w) / 3

n_p = [0] * 4
for i in range(4):
    n_p[i] = 90 * 0.25

xi = [0] * 4
for i in range(4):
    xi[i] = (vj[i] - n_p[i]) ** 2 / n_p[i]

# print("p:", p)
print("n_p:", n_p)
print("xi", xi)
print("xi observed:", np.sum(xi), '\n')

print("Task 2 - 3 (Uniform cont)")

sample1 = get_sample("task2.csv", 12, 30)
sample2 = get_sample("task2.csv", 13, 30)
sample3 = get_sample("task2.csv", 14, 30)

nij = [[0] * 4, [0] * 4, [0] * 4]
vj = [0] * 4

for j in range(4):
    if j == 0:
        left, right = -2, -1
    elif j == 1:
        left, right = -1, 0
    elif j == 2:
        left, right = 0, 1
    elif j == 3:
        left, right = 1, 2

    i = 0
    for smp in sample1, sample2, sample3:
        vj[j] += count_nij(smp, left, right)
        nij[i][j] = count_nij(smp, left, right)
        i += 1

xi = 0
for i in range(3):
    for j in range(4):
        xi += (nij[i][j] - vj[j] * 30 / 90) ** 2 / (vj[j] * 30)

xi *= 90

print("vj:", vj)
print("nij:", nij)
print("xi:", xi, '\n')

sample_w = sample1 + sample2 + sample3
# p = np.mean(sample_w) / 3

n_p = [22.5, 22.5, 22.5, 22.5]

xi = [0] * 4
for i in range(4):
    xi[i] = (vj[i] - n_p[i]) ** 2 / n_p[i]

# print("p:", p)
print("n_p:", n_p)
print("xi", xi)
print("xi observed:", np.sum(xi))
