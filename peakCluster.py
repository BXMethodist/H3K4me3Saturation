



class peakCluster:
    def __init__(self, chrNumber, start, end, path="", step=10):
        self.chrNumber = chrNumber
        self.start = start
        self.end = end
        self.step = step
        # this will be the default path where I store the peak information
        self.path = path

        self.clusterNameMap = {}
        self.clusterMap = {}


