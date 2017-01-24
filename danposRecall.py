

import os

def danposRecall(k, cutoff):
    wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]

    n = 0
    m = 0

    finishedjob = os.listdir("/home/tmhbxx3/archive/KFH3K4me3/" + str(cutoff) + "cutoff/pooled")
    finishedjob = [x.replace("archive_tmhkxc48_BroadH3K4me3_broadpeak201401_H3K4me3_dregion_pooled_", "") for x in
                   finishedjob]
    finishedjob = [x.replace(".peaks.xls", ".wig") for x in finishedjob]

    finishedjob = set(finishedjob)

    for wig in wigFiles:
        if k*60 <=n <= (k+1)*60:
            if wig not in finishedjob:
                cmd = "python /archive/tmhkxc48/tools/danposTemp/danpos.py dpeak "+wigPath+wig+" -q "+str(cutoff)+" -f 0 -z 0 -o /home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff"
                #os.system(cmd)
                #print cmd
                m +=1
                print wig
        n+=1
    print m

if __name__ == "__main__":
    danposRecall(0, 90)




