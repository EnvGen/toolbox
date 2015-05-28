#!/usr/bin/env python
import os, csv, sys

anns = {}
try: filelist = sys.argv[1]
except IndexError: sys.exit("Usage: collate.annotations.py <filelist>")
d = {}
files = []

for line in open(filelist, 'r'): files.append(line.rstrip())

for f in files:
    cl = (f.split("/")[-1]).split(".")[0]
    d[cl] = {}
    hin = open(f, 'r')
    hincsv = csv.reader(hin, delimiter = '\t')
    for row in hincsv:
        ann = row[1]
        anns[ann] = ""
        try: d[cl][ann] += 1
        except KeyError: d[cl][ann] = 1
    hin.close()

hout = sys.stdout
houtcsv = csv.writer(hout, delimiter = '\t')
houtcsv.writerow(["Family"]+sorted(d.keys()))

for ann in anns.keys():
    out = [ann]
    for cl in sorted(d.keys()):
        try: out.append(d[cl][ann])
        except KeyError: out.append(0)
    houtcsv.writerow(out)
hout.close()
