# Jan 11, 2022: enabled the use of reliability
# Jan 13, 2022: corrected genotype centering
# Mar 10, 2022: renamed SSGP (ssgp) to SLEMM (slemm)
# Jul 02, 2022: printed the last modified date
# Jul 30, 2022: improved the speed of associations

#!/usr/bin/python3
import pgenlib
import numpy as np
from scipy.stats import chi2
from scipy import linalg
import sys
import os
import time

'''
Important info for PgenReader overhead:

For a SET of variants, continuous or not,
1) repeat read() and convert int8 into float32 2d-array
2) repeat read_dosages()
3) use read_list() and convert int8 2d-array to float32
3 is slightly faster than 1 and 2.
Using int8 is faster than using int32.

For a number of continous variants,
1) use read() and convert int8 into float32 1d-array
2) use read_dosages()
3) use read_range() convert int8 2d-array to float32
1 is faster than 2 and 3, per variant.
'''

def get_parser():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--pfile", required=False,
                        dest="pfile",
                        help="PLINK2 pfile prefix",
                        metavar="<prefix>")
    parser.add_argument("-g", "--pgen", required=False,
                        dest="pgen",
                        help="PLINK2 pgen file",
                        metavar="<filename>")
    parser.add_argument("-v", "--pvar", required=False,
                        dest="pvar",
                        help="PLINK2 pvar file",
                        metavar="<filename>")
    parser.add_argument("-s", "--psam", required=False,
                        dest="psam",
                        help="PLINK2 psam file",
                        metavar="<filename>")
    parser.add_argument("-b", "--slemm", required=True,
                        dest="slemm",
                        help="SLEMM LMM filename prefix",
                        metavar="<prefix>")
    parser.add_argument("-o", "--out", required=True,
                        dest="out",
                        help="output filename",
                        metavar="<filename>")
    parser.add_argument("-f", "--maf",
                        dest="maf",
                        type=float,
                        default=0.0,
                        help="min minor allele frequency",
                        metavar="<freq>")
    parser.add_argument("-m", "--mac",
                        dest="mac",
                        type=int,
                        default=20,
                        help="min minor allele count",
                        metavar="<count>")
    parser.add_argument("-i", "--begin",
                        dest="begin",
                        type=int,
                        default=0,
                        help="begin index",
                        metavar="<index>")
    parser.add_argument("-j", "--end",
                        dest="end",
                        type=int,
                        default=1000000000,
                        help="end index",
                        metavar="<index>")
    return parser

def mean_wo_outliers(x):
    q1 = np.percentile(x, 25)
    q3 = np.percentile(x, 75)
    lower = q1 - (q3 - q1)*1.5
    upper = q3 + (q3 - q1)*1.5
    rsum = 0.0
    rcnt = 0.0
    for i in x:
        if i<upper and i>lower:
            rsum += i
            rcnt += 1
    if int(len(x)-rcnt)>0:
        print(int(len(x)-rcnt),"outliers identified when calculating correction factor")
        print("If not removing outliers, the correction factor is", sum(x)/len(x))
        print("Removed outliers")
    return(rsum/rcnt)


if __name__ == "__main__":
    start = time.time()
    
    print("SLEMM-Gamma by Jicai Jiang")
    print("Last Modified: Sat, 30 Jul 2022\n")
    
    args = get_parser().parse_args()
    # check arguments
    pgen = None
    psam = None
    pvar = None
    if args.pfile is not None:
        pgen = args.pfile + ".pgen"
        psam = args.pfile + ".psam"
        pvar = args.pfile + ".pvar"
        if not (os.path.isfile(pgen) and os.path.isfile(psam) and os.path.isfile(pvar)):
            print("Files do not exit:", pgen, psam, pvar)
            print("Searching for bed/bim/fam files.")
            pgen = args.pfile + ".bed"
            psam = args.pfile + ".fam"
            pvar = args.pfile + ".bim"
    elif (args.pgen is not None) and (args.psam is not None) and (args.pvar is not None):
        pgen = args.pgen
        psam = args.psam
        pvar = args.pvar
    else:
        print("Genotype files are not specified. Use following arguments.")
        print("    --pfile <prefix>")
        print("    --pgen <filename> --pvar <filename> --psam <filename>")
        sys.exit()
    if not (os.path.isfile(pgen) and os.path.isfile(psam) and os.path.isfile(pvar)):
        print("Files do not exit:", pgen, psam, pvar)
        sys.exit()
    
    phenofile = args.slemm + ".reml.py.txt"
    if not os.path.isfile(phenofile):
        print("File path {} does not exist.".format(phenofile))
        sys.exit()
    gsfile = args.slemm + ".gstat.txt"
    if not os.path.isfile(gsfile):
        print("File path {} does not exist.".format(gsfile))
        sys.exit()
    vcfile = args.slemm + ".reml.vc.csv"
    if not os.path.isfile(vcfile):
        print("File path {} does not exist.".format(vcfile))
        sys.exit()
    
    # SLEMM gstat file
    # Compute correction factor
    cfactor = None
    with open(gsfile) as fp:
        next(fp)
        r = []
        for line in fp:
            elem = line.strip().split()
            x = float(elem[7])
            y = float(elem[9])
            r.append(y/x)
        cfactor = mean_wo_outliers(r)
    print("Correction factor is", cfactor)
    # SLEMM vc file
    # Get variances
    snp_var = None
    err_var = None
    with open(vcfile) as fp:
        fp.readline()
        fp.readline()
        snp_var = float(fp.readline().strip().split(',')[1])
        err_var = float(fp.readline().strip().split(',')[1])
    print("SNP var is", snp_var)
    print("ERR var is", err_var)
    tau = err_var/snp_var
    
    '''
    Get phenotypes in SLEMM reml.py file and store 
    them in dict pheno. All the SLEMM samples are 
    in .psam file.
    '''
    pheno = {}
    covar = {}
    diagr = {}
    covar_ct = None
    with open(phenofile) as fp:
        line = fp.readline()
        covar_ct = len(line.strip().split()) - 3
        for line in fp:
            elem = line.strip().split()
            pheno.update({elem[0] : elem[1]})
            covar.update({elem[0] : elem[2:(2+covar_ct)]})
            diagr.update({elem[0] : elem[-1]})
    sample_ct = len(pheno)
    print("\nSample count in SLEMM file is", sample_ct)
    print("Number of covariates is",covar_ct)
    
    if 2*sample_ct*args.maf > args.mac:
        args.mac = int(2*sample_ct*args.maf + 0.5)
    
    '''
    Read .psam file.
    The order of samples in vec_y is the same 
    as in sample_subset.
    '''
    raw_sample_ct = 0
    sample_subset = np.empty(sample_ct, dtype=np.uint32)
    vec_y = np.empty(sample_ct, dtype=np.float32)
    mat_covar = np.empty([sample_ct,covar_ct], dtype=np.float32)
    vec_rr = np.empty(sample_ct, dtype=np.float32)
    with open(psam) as fp:
        header = None
        first_cline = None 
        for line in fp:
            if line[0] == '#':
                header = line
            else:
                first_cline = line
                break
        iid_col = 1
        if (header is not None) and (header[1] == 'I'):
            iid_col = 0
        idx = 0
        iid = first_cline.strip().split()[iid_col]
        if iid in pheno:
            sample_subset[idx] = raw_sample_ct
            vec_y[idx] = pheno[iid]
            mat_covar[idx,] = covar[iid]
            vec_rr[idx] = diagr[iid]
            idx += 1
        for line in fp:
            iid = line.strip().split()[iid_col]
            if iid is not None:
                raw_sample_ct += 1
            if iid in pheno:
                sample_subset[idx] = raw_sample_ct
                vec_y[idx] = pheno[iid]
                mat_covar[idx,] = covar[iid]
                vec_rr[idx] = diagr[iid]
                idx += 1
    raw_sample_ct += 1
    print("Raw sample count in psam is", raw_sample_ct)
    if sample_ct == raw_sample_ct:
        sample_subset = None
        print("Sample count in GWA is", sample_ct)
    else:
        print("Sample count in GWA is", sample_subset.size)
    
    is_one_rr = False
    if(np.min(vec_rr)>0.99999 and np.max(vec_rr)<1.00001):
        is_one_rr = True
    is_const_covar = False
    if(covar_ct==1 and np.max(mat_covar)==np.min(mat_covar)):
        is_const_covar = True
    
    # reciprocal square root of diagonal R
    vec_rr = 1.0/np.sqrt(vec_rr)
    mat_covar = np.transpose(np.transpose(mat_covar) * vec_rr)
    # qr decomposition of covariate matrix
    Q = linalg.qr(mat_covar, mode = "economic")[0]
    del mat_covar
    print("Completed QR decomposition of covariate matrix.")
    
    '''
    Read .pvar file.
    Get the number of variants to be tested.
    '''
    tested_var_ct = 0
    with open(pvar) as fp:
        var_idx = -1
        for line in fp:
            if line[0] == '#':
                continue
            if line.strip() is not None:
                var_idx += 1
                if var_idx > args.end:
                    break
                if var_idx >= args.begin:
                    tested_var_ct += 1
    
    '''
    Read pgen/bed file and do single-variant assoc tests.
    Keep intermediate terms, xy and xx.
    '''
    print("\nStarting assoc tests...")
    int8vec_x = np.empty(sample_ct, dtype=np.int8)
    vec_x = np.empty(sample_ct, dtype=np.float32)
    xx = np.empty(tested_var_ct, dtype=np.float32)
    xy = np.empty(tested_var_ct, dtype=np.float32)
    pf = pgenlib.PgenReader(pgen.encode(), raw_sample_ct = raw_sample_ct, sample_subset = sample_subset)
    for list_idx in range(tested_var_ct):
        var_idx = list_idx + args.begin
        if var_idx % 1000 == 0:
            print("    Scanning variant", var_idx)
        pf.read(var_idx, int8vec_x)
        vec_x = int8vec_x.astype('float32')
        vec_x_sum = vec_x.sum()
        if vec_x_sum < args.mac or (2*sample_ct - vec_x_sum) < args.mac:
            xy[list_idx] = np.nan
            xx[list_idx] = np.nan
            continue
        vec_x -= vec_x_sum/sample_ct
        xy[list_idx] = np.dot(vec_x, vec_y)
        if(not is_one_rr):
            vec_x *= vec_rr
        if(not (is_one_rr and is_const_covar)):
            vec_x -= np.dot(Q, np.dot(np.transpose(Q), vec_x))
        xx[list_idx] = np.dot(vec_x, vec_x)
    
    # Compute final statistics for report.
    beta = xy/xx * cfactor
    se_sq = (snp_var*cfactor)/xx
    chisq = beta * beta / se_sq
    pval = chi2.sf(chisq,1)
    print("Computed association statistics.")
    
    # Write stats to output file.
    print("\nWriting stats into", args.out)
    outfp = open(args.out, "w")
    with open(pvar) as fp:
        var_idx = -1
        for line in fp:
            if line[0] == '#':
                continue
            if line.strip() is not None:
                var_idx += 1
                if var_idx > args.end:
                    break
                if var_idx >= args.begin:
                    list_idx = var_idx - args.begin
                    outfp.write(line.strip() +"\t"+ str(beta[list_idx]) +"\t"+ str(se_sq[list_idx]) +"\t"+ str(chisq[list_idx]) +"\t"+ str(pval[list_idx]) +"\n")
    outfp.close()
    
    print(int(time.time()-start), "s taken for the whole analysis.\n")

