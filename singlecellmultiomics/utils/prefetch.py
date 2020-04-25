#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from copy import copy

class Prefetcher(object):

    def prefetch(self,contig, start, end):

        prefetched_clone = copy(self)
        prefetched_clone._prefetch(contig, start, end)
        return prefetched_clone
