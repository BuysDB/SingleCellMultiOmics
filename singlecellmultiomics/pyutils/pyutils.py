#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import subprocess
def meanOfCounter(counter):
    sum_of_numbers = sum(number*count for number, count in counter.most_common())
    count = sum(count for n, count in counter.most_common())
    if count==0:
        return 0
    return sum_of_numbers / count

def varianceOfCounter(counter):

    sum_of_numbers = sum(number*count for number, count in counter.most_common())
    if sum_of_numbers==0:
        return 0
    count = sum(count for n, count in counter.most_common())
    total_squares = sum(number*number * count for number, count in counter.most_common())
    mean_of_squares = total_squares / count
    mean = sum_of_numbers / count
    variance = mean_of_squares - mean * mean
    return math.sqrt(variance)


def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

def wccountgz(filename):
    try:
        cmd = f'zcat {filename} | wc -l'
        out = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True
                             ).communicate()[0]
        return int(out.partition(b' ')[0])
    except Exception as e:
        return 0
