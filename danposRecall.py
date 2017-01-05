# This is one is to call Danpos to run dpeak with different peak height cutoff


import os


# def danposRecall(k, cutoff):
#     wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"
#
#     wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]
#
#
#     n = 0
#     for wig in wigFiles:
#         if 180*k <=n< 180*(k+1):
#             cmd = "python /archive/tmhkxc48/tools/danposTemp/danpos.py dpeak "+wigPath+wig+" -q "+str(cutoff)+" -f 0 -z 0 -o /home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff"
#             os.system(cmd)
#         n+=1




# danposRecall(1, 500)

# bed1 = "/archive/tmhkxc48/BroadH3K4me3/test/ENCFF218WSN.chr21.bed"
# bed2 = "/archive/tmhkxc48/BroadH3K4me3/test/ENCFF621OIP.chr21.bed"
#
#
# cmd1 = "python /archive/tmhkxc48/tools/danposTemp/danpos.py dpeak "
# cmd2 = " -c 400000 --extend 200 --frsz 200 -q "
#
#
# for i in range(10, 310, 10):
#     os.system("mkdir "+str(i))
#
#     cmd = cmd1 + bed1 + cmd2 + str(i) +" -o /home/tmhbxx3/archive/test/"+str(i)
#     os.system(cmd)
#
#     cmd = cmd1 + bed2 + cmd2 + str(i) +" -o /home/tmhbxx3/archive/test/"+str(i)
#     os.system(cmd)






