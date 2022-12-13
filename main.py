from functools import reduce
from math import sqrt

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn


def read_table():
    df = pd.read_table("ekspresije.tsv", index_col=0)
    cols_t = df.columns.size
    rows_t = df.transpose().columns.size

    temp = [(cell, gene, value) for cell in df.columns
            for gene, value in df[cell].items()]
    return rows_t, cols_t, temp


def draw():
    embedding = pd.read_table('umap.tsv')
    embedding['cluster'] = 0
    plt.figure(figsize=(8, 6))
    plt.scatter(
        embedding.umap1,
        embedding.umap2,
        c=[sn.color_palette()[x] for x in embedding.cluster]
    )


def avg_func(ans, tmp):
    global avg_val, cnt, cols

    cnt += 1

    if cnt >= rows:
        avg_val[tmp[0]] = ans[2] / rows
        ans = (None, None, 0)
        cnt = 0

    return None, None, ans[2] + tmp[2]


def var_func(ans, tmp):
    global avg_val, cnt, cols, i, flag

    cnt += 1

    if cnt >= rows:
        var_val[f"cell_{i}"] = ans / (rows - 1)
        ans = 0
        cnt = 0
        i += 1

    if flag:
        flag = False
        return ans ** 2 + tmp ** 2

    return ans + tmp ** 2


def standardize_func(val):
    global standard_dev, avg_val, cnt, i, rows

    if cnt >= rows:
        cnt = 0
        i += 1

    cnt += 1

    return f"cell_{i}", val[1], (val[2] - avg_val[f"cell_{i}"])/standard_dev[i]

#def twopointone:

#    for j in

rows, cols, data = read_table()

avg_val = {}
var_val = {}
cnt = 1
flag = True
i = 0

reduce(avg_func, data)

cnt = 1

centered_val = map(lambda a: a[2] - avg_val[a[0]], data)

reduce(var_func, centered_val)

standard_dev = list(map(lambda a: sqrt(a), var_val.values()))

i = 0
cnt = 0

standardized_val = list(map(standardize_func, data))

std_val_by_gene = sorted(standardized_val, key=lambda x: x[1])
data_by_gene = sorted(data, key=lambda x: x[1])

temp = rows
rows = cols
temp_avg = avg_val
temp_var = var_val
flag = True

reduce(avg_func, data_by_gene)
centered_val_by_gene = list(map(lambda a: a[2] - avg_val[a[0]], data_by_gene))
reduce(var_func, centered_val_by_gene)

avg_val_by_gene = avg_val
var_val_by_gene = var_val
avg_val = temp_avg
var_val = temp_var
rows = temp

print(var_val_by_gene)
