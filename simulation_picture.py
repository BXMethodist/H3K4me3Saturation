from simulation_generator import get_random_peak
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


def simulation_illustration():
    locations =[200, 400, 600, 800]
    width = [180, 160, 80, 120, 150, 100]


    fig, axes = plt.subplots(5, 1, sharex=True)
    for i in range(5):
        signals1 = np.zeros(1000)
        real_peak_number = np.random.randint(2)
        for j in range(real_peak_number+1):
            peak_index = locations[np.random.randint(4)]

            while signals1[peak_index] != 0:
                peak_index = locations[np.random.randint(4)]
            peak_width = width[np.random.randint(6)]
            peak_height = np.random.randint(250,500)

            signals1[peak_index-peak_width/2:peak_index+peak_width/2] = get_random_peak(peak_height, peak_width)

        noise_peak_number = np.random.randint(5)
        signals2 = np.zeros(1000)
        for k in range(noise_peak_number+1):
            aval = np.where(signals1 ==0)[0]

            aval = aval[aval<970]

            index = aval[np.random.randint(len(aval))]

            peak_width = np.random.randint(20, 30)

            peak_height = np.random.randint(5, 200)

            signals2[index :index + peak_width] = get_random_peak(peak_height, peak_width)

        ax = axes[i]
        x = np.arange(1000)
        ax.fill_between(x, 0, signals1, color='red')
        ax.fill_between(x, 0, signals2, color='blue')
        ax.set_yticks([0, 200, 400])
        ax.set_ylim(0, 500)

        n = 0
        cutoffs = [n for n in range(0 ,310, 100)]
        # print cutoffs
        color = cm.rainbow(np.linspace(0, 1, len(cutoffs)))
        for m, c  in zip(range(len(cutoffs)-1), color):
            cutoff1 = np.asarray([cutoffs[m]]*1000)
            cutoff2 = np.asarray([cutoffs[m+1]]*1000)
            # print cutoff1, cutoff2
            ax.fill_between(x, cutoff1, cutoff2, color=c, alpha=0.4)

    plt.show()

    fig.savefig('simulated_illustration.pdf', dpi=600, format='pdf', facecolor='w', edgecolor='w',
                    orientation='portrait', bbox_inches='tight')

def plot_simulation_distribution(distribution, outputpath, color):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.hist(distribution, 20, facecolor=color, alpha=0.75)
    fig.savefig(outputpath, format='pdf', dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', bbox_inches='tight')
    plt.close('all')

Hreal = np.random.negative_binomial(4, 0.01, 100)
Hnoise = np.random.negative_binomial(1, 0.1, 150)

Wreal = np.random.negative_binomial(3, 0.05, 100)
Wnoise = np.random.negative_binomial(1, 0.05, 150)

distributions = [(Hreal, 'simulated_distribution_HeightReal.pdf', 'red'),
                 (Hnoise, 'simulated_distribution_HeightNoise.pdf', 'blue'),
                 (Wreal, 'simulated_distribution_WidthReal.pdf', 'red'),
                 (Wnoise, 'simulated_distribution_WidthNoise.pdf', 'blue')]

for distribution in distributions:
    d, output, color = distribution
    plot_simulation_distribution(d, output, color)



