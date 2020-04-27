#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from copy import copy

class Prefetcher(object):

    def prefetch(self,contig, start, end):

        prefetched_clone = copy(self)
        prefetched_clone._prefetch(contig, start, end)
        return prefetched_clone


class UnitialisedClass(object):

    def __init__(self, class_name, *args, **kwargs):
        self.class_name = class_name
        self.pos_arguments = args if args is not None else ()
        self.keyword_arguments = kwargs if kwargs is not None else dict()

    def initialise(self) -> object:
        return self.class_name(*self.pos_arguments, **self.keyword_arguments)


def initialise(o):
    if getattr(o, 'initialise', None) is not None:
        return o.initialise()
    return o


def initialise_dict(d: dict):
    """ Initialise all unitialised classes present in the dictionary d"""
    new = {}
    for key, value in d.items():
        if type( value ) is UnitialisedClass:
            new[key] = d[key].initialise()
        else:
            new[key] = value
    return new
