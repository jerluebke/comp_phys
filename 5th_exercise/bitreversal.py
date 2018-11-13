# -*- coding: utf-8 -*-

def reverse_stdlib(x):
    return int(bin(x)[:1:-1], 2)

def reverse_bit(x, s):
    r = 0
    for i in range(s):
        r = (r << 1) | (x & 1)
        x >>= 1
    return r
