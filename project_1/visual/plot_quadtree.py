# -*- coding: utf-8 -*-

from visualise import PyQuadtreeEnv

qtenv = PyQuadtreeEnv() # TODO: add data

while 1:
    try:
        qtenv.insert_next()
        # TODO: do plotting here
    except StopIteration:
        break


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
