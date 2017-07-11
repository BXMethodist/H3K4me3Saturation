

import os

def danposRecall():
    wigPath = "/home/tmhbxx3/scratch/ENC_H3K4me3/wigs_with_input/"

    wigFiles = [wigPath + path for path in os.listdir(wigPath) if path.endswith("wig")]

    cutoffs = [str(x) for x in [3,6] + range(10, 310, 10)]

    for wig in wigFiles[1:]:
        danpos_no_input(wig, cutoffs)

def danpos_no_input(sample_id, cutoffs):

    sample_name = sample_id[sample_id.rfind('/') + 1:sample_id.find('.')]
    danpos_cmd = 'python /home/tmhbxx3/archive/tools/danposTemp_multi_q/danpos.py dpeak '
    danpos_parameters = " -q "+','.join(cutoffs) +" -f 0 --smooth_width 0 -o /home/tmhbxx3/scratch/ENC_H3K4me3/peaks"

    cmd = danpos_cmd + sample_id + danpos_parameters

    pbs = open(sample_name + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + sample_id +'\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    # pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.2.1\n")
    pbs.write(cmd + '\n')
    pbs.close()
    os.system('qsub ' + sample_name + ".pbs")
    return


if __name__ == "__main__":
    danposRecall()


