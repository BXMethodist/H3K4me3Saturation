import numpy as np


class WigChrom:
    ## This is a class for storing the chromosome information from a wig file

    def __init__(self, chr_name, start, size, step, span, fixed=True):
        self.chr_name = chr_name
        self.start = start
        self.size = size
        self.step = step
        self.span = span
        self.fixed = fixed

        vector_size = self.size / step if self.size % step == 0 else self.size/step + 1
        self.signals = np.zeros(vector_size)

    def get_signals(self, start, end):
        if start < 0:
            print "start position need to be bigger than 0"
            return
        if end > self.size:
            print "end position need to be smaller than genome size"
            return
        if end%self.step == 0:
            return self.signals[start/self.step:end/self.step]
        else:
            return self.signals[start/self.step:end/self.step+1]
