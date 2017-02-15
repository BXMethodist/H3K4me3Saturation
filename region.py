class region():
    """
    Store the inforamtion of a reference peak:
        Parameters
        ----------
        chromosome: chromosome name
        start: start position
        end: end position
        step: step or bin size use to store the signal

        Attributes
        ----------
        chromosome name
        start
        end
        signals: np.array that keep the peaks signal
        step
        """

    def __init__(self, chromosome, start, end, signals, peaks, step=10, cutoff=200):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.signals = signals
        self.step = step
        self.peaks = peaks
        self.cutoff = cutoff

