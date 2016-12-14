# This is one is to call Danpos to run dpeak with different peak height cutoff


import os


def danposRecall(k, cutoff):
    wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]


    n = 0
    for wig in wigFiles:
        if 180*k <=n< 180*(k+1):
            cmd = "python /archive/tmhkxc48/tools/danposTemp/danpos.py dpeak "+wigPath+wig+" -q "+str(cutoff)+" -f 0 -z 0 -o /home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff"
            os.system(cmd)
        n+=1




danposRecall(1, 500)
