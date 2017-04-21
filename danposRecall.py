

import os

def danposRecall():
    wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

    wigFiles = [wigPath + path for path in os.listdir(wigPath) if path.endswith("wig")]

    for cutoff in range(10, 310, 10) + [3, 6]:
        os.system('mkdir '+str(cutoff))


    for wig in wigFiles[:5]:
        for cutoff in [100]:
            danpos_no_input(wig, cutoff)

def danpos_no_input(sample_id, cutoff):
    'python /archive/tmhkxc48/tools/danposTemp/danpos.py dpeak /archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/BI.Adipose_Nuclei.H3K4me3.92.Fnor.wig -q 10 -f 0 -z 0 -o /home/tmhbxx3/archive/test_danpos/10'

    sample_name = sample_id[sample_id.rfind('/') + 1:]
    danpos_cmd = 'python /archive/tmhkxc48/tools/danposTemp/danpos.py dpeak '
    danpos_parameters = " -q "+str(cutoff)+" -f 0 -z 0 -o /home/tmhbxx3/archive/test_danpos/"+str(cutoff) + "/" +sample_name

    cmd = danpos_cmd + sample_id + danpos_parameters

    pbs = open(sample_name + str(cutoff)+ ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + sample_id + str(cutoff)+'\n')
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
    os.system('qsub ' + sample_name + str(cutoff)+ ".pbs")
    return


if __name__ == "__main__":
    danposRecall()


