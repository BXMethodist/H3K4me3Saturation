import numpy as np, pandas as pd, os
from ReferenceMap import refMap
from table import *
from CallPeak import callpeak
import matplotlib.pyplot as plt


# #

def simulation_genome_size_generater(genome_size):
    f = open("simulation_genome_sizes.txt", "w")
    f.write("simulation\t"+str(genome_size))
    f.close()
    path = os.getcwd()
    if not path.endswith("/"):
        path += "/"
    return path + 'simulation_genome_sizes.txt'

def get_negative_binomial_random_samples(n, p, size):
    """
    random sample generator for negative binomial distribution
    :param mean: distribution mean
    :param p: probability
    :return: 1-d array containing numbers
    """
    return np.random.negative_binomial(n, p, size)

def get_random_locations(locations, peak_number):
    """
    get random locations for peaks
    :param locations:
    :param peak_number:
    :return: a 1-d array containing the random index
    """
    return locations[np.random.randint(low=0, high=len(locations), size=peak_number)]

def get_random_peak(height, width):
    """
    create a random shape of a peak
    :param width: peak width
    :return:
    """
    x = np.tri(height, width, int(width/height))
    x = np.sum(x, axis=0)
    return x
    # if width <= 80:
    #     return np.sin(np.linspace(0, np.pi, width))
    # else:
    #     reads_number = width/20 * 100
    #     reads_length = 20
    #     signals = np.zeros(width)
    #     for i in range(reads_number):
    #         left_increment = np.random.randint(0, i+1)
    #         right_increment = np.random.randint(0, i+1)
    #         if left_increment >= width-20-right_increment:
    #             break
    #         start_index = np.random.randint(left_increment, width-20-right_increment)
    #         signals[start_index: start_index+reads_length]+=1
    #     return signals

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
              genome_size=1000000, real_peak_number=200,
              noise_peak_number=900, real_peak_locations=None, sample_number=300,
              cutoffs=[x for x in [6]+ range(10, 310, 10)],
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

    samples = []
    for i in range(sample_number):
        # generate simulated real peaks
        cur_real_peak_locations = get_random_locations(real_peak_locations, real_peak_number)
        cur_real_peak_widths = get_negative_binomial_random_samples(W_real, wp_real, real_peak_number)
        cur_real_peak_heights = get_negative_binomial_random_samples(H_real, p_real, real_peak_number)

        # generate simulated noise peaks
        noise_peak_locations = get_random_locations(np.arange(genome_size), noise_peak_number)
        noise_peak_heights = get_negative_binomial_random_samples(H_noise, p_noise, noise_peak_number)
        noise_peak_widths = get_negative_binomial_random_samples(W_noise, wp_noise, noise_peak_number)

        signals = np.zeros(genome_size)

        for j in range(real_peak_number):
            cur_width = int(cur_real_peak_widths[j])
            cur_height = int(cur_real_peak_heights[j])
            cur_location = int(cur_real_peak_locations[j])
            cur_width = cur_width/2*2

            if cur_width == 0 or cur_height==0:
                continue
            cur_signals = get_random_peak(cur_height, cur_width)
            if cur_location - (cur_width / 2) >= 0 and cur_location + (cur_width / 2) < signals.shape[0]:
                signals[cur_location-(cur_width/2):cur_location+(cur_width/2)] = cur_signals*(cur_height/np.max(cur_signals))
            elif cur_location - (cur_width / 2) < 0 and cur_location + cur_width < signals.shape[0]:
                signals[cur_location:cur_location + (cur_width)] = cur_signals * (cur_height/np.max(cur_signals))
            elif cur_location - (cur_width) >= 0 and cur_location + cur_width/2 >= signals.shape[0]:
                signals[cur_location-cur_width:cur_location] = cur_signals * (cur_height/np.max(cur_signals))

        for k in range(noise_peak_number):
            cur_width = int(noise_peak_widths[k])
            cur_height = int(noise_peak_heights[k])
            cur_location = int(noise_peak_locations[k])
            cur_width = cur_width/2*2
            if cur_width == 0 or cur_height==0:
                continue
            cur_signals = get_random_peak(cur_height, cur_width)
            if cur_location - (cur_width / 2) >= 0 and cur_location + (cur_width / 2) < signals.shape[0]:
                signals[cur_location-(cur_width/2):cur_location+(cur_width/2)] = cur_signals*(cur_height/np.max(cur_signals))
            elif cur_location - (cur_width / 2) < 0 and cur_location + cur_width < signals.shape[0]:
                signals[cur_location:cur_location + (cur_width)] = cur_signals * (cur_height/np.max(cur_signals))
            elif cur_location - (cur_width) >= 0 and cur_location + cur_width/2 >= signals.shape[0]:
                signals[cur_location-cur_width:cur_location] = cur_signals * (cur_height/np.max(cur_signals))
        samples.append(signals)

        # signals = signals[1000:3000]
        # plt.plot(np.arange(2000), signals)
        # plt.show()

        for cutoff in cutoffs:
            genome_signals = {'simulation':signals}
            results = callpeak(genome_signals, cutoff, step=1)
            df = pd.DataFrame(results, index=None, columns=None)
            df.to_csv(simulation_directory+ str(cutoff)+'/'+'simulated_sample'+str(i)+'cutoff'+str(cutoff)+".tsv",
                      sep='\t', index=None, header=None)
    return

if __name__ == "__main__":
    H_reals = [2]
    p_reals = [0.003,0.004,0.005]
    # p_reals = [x/1000.0 for x in range(1, 1000)]
    H_noises = [1]
    p_noises = [0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]
    # p_noises = [x/100.0 for x in range(1, 100)]
    W_reals = [1,2,3,4]
    wp_reals = [0.01, 0.02,0.03,0.04,0.05]
    W_noises = [1,2]
    wp_noises = [0.02,0.04,0.06,0.08,0.1,0.15,0.2,0.25,0.3]
    genome_size = 1000000
    real_peak_numbers = [100]
    noise_peak_numbers =[400]
    # real_peak_numbers = [x for x in range(100,310, 100)]
    # noise_peak_numbers = [x for x in range(100,1510,100)]
    # step 1
    parameters = [(H_real, p_real, H_noise, p_noise, W_real, wp_real, W_noise, wp_noise, real_peak_number, noise_peak_number)
                  for H_real in H_reals
                  for p_real in p_reals
                  for H_noise in H_noises
                  for p_noise in p_noises
                  for W_real in W_reals
                  for wp_real in wp_reals
                  for W_noise in W_noises
                  for wp_noise in wp_noises
                  for real_peak_number in real_peak_numbers
                  for noise_peak_number in noise_peak_numbers
                  if real_peak_number <= noise_peak_number
                  if H_real >= H_noise
                  if p_real <= p_noise
                  if W_real >= W_noise
                  if wp_real <= wp_noise]

    print len(parameters)

    # H_real = 2
    # p_real = 0.003
    # # H_real = 2
    # # p_real = 1
    # H_noise = 1
    # p_noise = 0.15
    # W_real = 3
    # wp_real = 0.01
    # W_noise = 1
    # wp_noise = 0.01
    # genome_size = 1200000
    for parameter in parameters[:10000]:
        H_real, p_real, H_noise, p_noise, W_real, wp_real, W_noise, wp_noise, real_peak_number, noise_peak_number = \
        parameter
        generater(H_real=H_real, p_real=p_real,
                  H_noise=H_noise, p_noise=p_noise,
                  W_real=W_real, W_noise=W_noise,
                  wp_real=wp_real,
                  wp_noise=wp_noise,
                  genome_size=genome_size,
                  real_peak_number=real_peak_number,
                  noise_peak_number=noise_peak_number)
        genome_file_path = simulation_genome_size_generater(genome_size)
        #
        simulation_path = "/home/tmhbxx3/archive/WigChrSplits/code/simulation_results/"

        for cutoff in [6]+range(10, 310, 10):
            cur_refmap = refMap(1, genome_size_path=genome_file_path)
            # print cutoff, 'is start'
            # cur_refmap.trainMap("/home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff/pooled", cutoff=cutoff,
            #                 individual=True)

            cur_refmap.trainMap(simulation_path + str(cutoff),
                                outputname='simulation', cutoff=cutoff,
                                individual=False, saveRefMap=False)

        refmap_path = "/home/tmhbxx3/archive/WigChrSplits/code/"
        outputname = 'simulation'+'_H_real_'+str(H_real) + '_p_real_'+str(p_real)+\
                     "_H_noise_"+str(H_noise)+"_p_noise_" + str(p_noise) \
                     + "_W_real_" + str(W_real) + "_W_noise_" + str(W_noise) + \
                     "_wp_real_" + str(wp_real) + "_wp_noise_" + str(wp_noise)
        df = finalpoint_cutoff_vs_stat(cutoffs=[x for x in [6]+range(10, 310, 10)],
                                  file_addresses=refmap_path,
                                  number_sample=300, outputname=outputname, prefix='simulation')

        sub_df = df.ix[1:, :]
        print df['Length of Region'].argmin(), type(df['Length of Region'].argmin())
        print sub_df['Length of Region'].argmax(), type(sub_df['Length of Region'].argmax())
        print df.ix[10, 'Length of Region'], type(df.ix[10, 'Length of Region'])
        print df.ix[sub_df['Length of Region'].argmax(), 'Length of Region']
        print parameter
        if df['Length of Region'].argmin() == 10:
            if 70<= sub_df['Length of Region'].argmax() <= 110:
                if 100 <= df.ix[10, 'Length of Region'] < 300 and \
                    600 < df.ix[100, 'Length of Region'] < 1000:
                    df.to_csv(outputname)