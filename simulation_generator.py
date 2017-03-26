import numpy as np, pandas as pd, os
from ReferenceMap import refMap
from table import *

#
# s = np.random.negative_binomial(2, 0.2, 100000)
# print np.mean(s)
# plt.hist(s, bins=100)
# plt.show()

def simulation_genome_size_generater(genome_size):
    f = open("simulation_genome_sizes.txt", "w")
    f.write("simulation\t"+str(genome_size))
    f.close()
    path = os.getcwd()
    if not path.endswith("/"):
        path += "/"
    return path + 'simulation_genome_sizes.txt'

def get_negative_binomial_random_samples(mean, p, size):
    """
    random sample generator for negative binomial distribution
    :param mean: distribution mean
    :param p: probability
    :return: 1-d array containing numbers
    """
    n = int(p*mean/(1-p))
    return np.random.negative_binomial(n, p, size)

def get_random_locations(locations, peak_number):
    """
    get random locations for peaks
    :param locations:
    :param peak_number:
    :return: a 1-d array containing the random index
    """
    return locations[np.random.randint(low=0, high=len(locations), size=peak_number)]

def draw_simulated_sample(ndarray, start, end):
    """
    :param ndarray: signals
    :param start: start index
    :param end: end index
    :return: None, draw the picture
    """
    cur_array = ndarray[start:end]
    df = pd.DataFrame(cur_array, columns=['value'])
    df = df.set_index(np.arange(len(cur_array)))
    ax = df.plot(y=['value'], kind='area')
    plt.show()
    return

def generater(H_real, p_real,  H_noise, p_noise, W_real, W_noise, wp_real, wp_noise,
              genome_size=100000, real_peak_number=250,
              noise_peak_number=1000, real_peak_locations=None, sample_number=300,
              cutoffs=[x for x in range(10, 310, 10)],
              simulation_directory='./simulation_results/'):
    """
    :param H_real:
    :param p_real:
    :param H_noise:
    :param p_noise:
    :param W_real:
    :param W_noise:
    :param genome_size:
    :param real_peak_number:
    :param noise_peak_number:
    :param real_peak_locations:
    :param noise_peak_locations:
    :param sample_number:
    :return: simulated_sample_peaks_table, simulated_genome_size_file
    """
    if not simulation_directory.endswith("/"):
        simulation_directory += "/"

    for cutoff in cutoffs:
        if not os.path.isdir(simulation_directory+str(cutoff)):
            os.system("mkdir "+simulation_directory+str(cutoff))

    if real_peak_locations is None:
        real_peak_locations = get_random_locations(np.arange(genome_size), real_peak_number*2)

    for i in range(sample_number):
        # generate simulated real peaks
        cur_real_peak_locations = get_random_locations(real_peak_locations, real_peak_number)
        cur_real_peak_widths = get_negative_binomial_random_samples(W_real, wp_real, real_peak_number)
        cur_real_peak_heights = get_negative_binomial_random_samples(H_real, p_real, real_peak_number)

        # generate simulated noise peaks
        noise_peak_locations = get_random_locations(np.arange(genome_size), noise_peak_number)
        noise_peak_heights = get_negative_binomial_random_samples(W_noise, wp_noise, noise_peak_number)
        noise_peak_widths = get_negative_binomial_random_samples(H_noise, p_noise, noise_peak_number)

        for cutoff in cutoffs:
            results = []
            for j in range(real_peak_number):
                cur_width = int(cur_real_peak_widths[j])
                cur_height = int(cur_real_peak_heights[j])
                cur_location = int(cur_real_peak_locations[j])

                if cur_height >= cutoff:
                    results.append(('simulation', cur_location-cur_width/2, cur_location+cur_width/2,
                                   0,0,0,cur_height))

            for k in range(noise_peak_number):
                cur_width = int(noise_peak_widths[k])
                cur_height = int(noise_peak_heights[k])
                cur_location = int(noise_peak_locations[k])
                if cur_height >= cutoff:
                    results.append(('simulation', cur_location-cur_width/2, cur_location+cur_width/2,
                                   0, 0, 0, cur_height))

            df = pd.DataFrame(results, index=None, columns=None)
            df.to_csv(simulation_directory+ str(cutoff)+'/'+'simulated_sample'+str(i)+'cutoff'+str(cutoff)+".tsv",
                      sep='\t', index=None, header=None)
    return

if __name__ == "__main__":
    # step 1
    H_real = 100
    p_real = 0.01
    H_noise = 30
    p_noise = 0.3
    W_real = 100
    wp_real = 0.05
    W_noise = 10
    wp_noise = 0.1

    generater(H_real=H_real, p_real=p_real,
              H_noise=H_noise, p_noise=p_noise,
              W_real=W_real, W_noise=W_noise,
              wp_real=wp_real,
              wp_noise=wp_noise)
    genome_file_path = simulation_genome_size_generater(100000)
    #
    simulation_path = "/home/tmhbxx3/archive/WigChrSplits/code/simulation_results/"

    for cutoff in range(10, 310, 10):
        cur_refmap = refMap(1, genome_size_path=genome_file_path)
        print cutoff, 'is start'
        # cur_refmap.trainMap("/home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff/pooled", cutoff=cutoff,
        #                 individual=True)

        cur_refmap.trainMap(simulation_path + str(cutoff),
                            outputname='simulation', cutoff=cutoff,
                            individual=False, saveRefMap=False)

    refmap_path = "/home/tmhbxx3/archive/WigChrSplits/code/"
    outputname = 'simulation'+'_H_real_'+str(H_real) + '_p_real_'+str(p_real)+\
                 "_H_noise_"+str(H_noise)+"_p_noise_" + str(p_noise) \
                 + "_W_real_" + str(W_real) + "_W_noise_" + str(W_noise)
    finalpoint_cutoff_vs_stat(cutoffs=[x for x in range(10, 310, 10)],
                              file_addresses=refmap_path,
                              number_sample=300, outputname=outputname, prefix='simulation')