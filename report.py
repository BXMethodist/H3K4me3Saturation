## This is module is used to check running result from the standard out and standard error from different analysis software.


import os, pandas as pd, numpy as np


def bowtie_results(path, surffix='.error'):
    """
    :param path: the path containing bowtie results in standard error
           surffix: usurally is named as '.error', if not please specify
    :return: a dataframe table that containing error reports, and will print out the list if the error report is empty
    which happened a lot when using sra toolkit
    """
    error_files = [x for x in os.listdir(path) if x.find(surffix) != -1]
    print error_files[0]
    path = path + '/' if not path.endswith("/") else path

    results = {}
    failed = []
    for error_file in error_files:
        f = open(path+error_file, "r")
        info = f.readlines()
        f.close()
        error = True
        filename = error_file[:error_file.find(surffix)]
        filename = filename.replace('bowtie', '')
        filename = filename.replace("_", "")
        for line in info:
            if line.startswith("# reads with at least one reported alignment:"):
                error = False
                percentage = line[line.find("(")+1:line.find(")")]
                if filename in results:
                    results[filename] = max(results[filename], percentage)
                else:
                    results[filename] = percentage
            if line.find("storage exhausted while writing file within file system module") != -1:
                print error_file
        if error:
            failed.append(filename)

    final_results = []
    for key, value in results.items():
        final_results.append((key, value))

    df = pd.DataFrame(final_results, columns=['sample_ID','percentage'])

    re_failed = []
    for fail in failed:
        if fail not in df['sample_ID'].unique():
            re_failed.append(fail)
    print re_failed
    df.to_csv('bowtie_results.csv', index=False)
    return df


bowtie_results('./', surffix='.e')
