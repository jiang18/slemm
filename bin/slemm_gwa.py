#!/usr/bin/python3

# Jan 11, 2022: enabled the use of reliability
# Jan 13, 2022: corrected genotype centering
# Jan 26, 2022: reduced block size if it is > # of model SNPs on a chr
# Mar 10, 2022: renamed SSGP (ssgp) to SLEMM (slemm)
# May 26, 2022: forced block size to an even number if > # of model SNPs
# Jul 02, 2022: printed the last modified date
# Jul 30, 2022: improved the speed of associations
# Jul 31, 2022: added a header line in output
# Mar 02, 2023: made --maf/--mac behave like PLINK's

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
    parser.add_argument("-c", "--chr", required=True,
                        dest="chr",
                        help="chromosome",
                        metavar="<chr>")
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
    return parser


if __name__ == "__main__":
    start = time.time()
    
    print("*******************************************************************")
    print("* SLEMM-GWA by Jicai Jiang")
    print("* Last Modified: Thu, 02 Mar 2023")
    print("*******************************************************************")
    
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
    snpfile = args.slemm + ".reml.snp.csv"
    if not os.path.isfile(snpfile):
        print("File path {} does not exist.".format(snpfile))
        sys.exit()
    
    # SLEMM gstat file
    # Compute correction factor
    cfactor = None
    block_size = None
    with open(gsfile) as fp:
        next(fp)
        # sum_xx = 0.0
        # sum_xy = 0.0
        sum_r = 0.0
        num_r = 0
        for line in fp:
            elem = line.strip().split()
            x = float(elem[7])
            y = float(elem[8])
            block_size = int(elem[10])
            # sum_xx += x*x
            # sum_xy += x*y
            sum_r += y/x
            num_r += 1
        # cfactor = sum_xy/sum_xx
        cfactor = sum_r/num_r
    print("Correction factor is", cfactor)
    # SLEMM vc file
    # Get variances
    model_snp_ct = None
    snp_var = None
    err_var = None
    with open(vcfile) as fp:
        fp.readline()
        model_snp_ct = int(fp.readline().strip().split(',')[1])
        snp_var = float(fp.readline().strip().split(',')[1])
        err_var = float(fp.readline().strip().split(',')[1])
    print("Block size is", block_size)
    print("SNP var is", snp_var)
    print("ERR var is", err_var)
    tau = err_var/snp_var
    
    # SLEMM snp file
    ksnp_set = set()
    with open(snpfile) as fp:
        next(fp)
        for line in fp:
            elem = line.strip().split(',')
            if elem[1] == args.chr:
                ksnp_set.add(elem[0])
    ksnp_ct = len(ksnp_set)
    print("\nCount of model variants on Chr", args.chr, "is", ksnp_ct)
    if ksnp_ct == 0:
        print("Error: Is --chr", args.chr, "correctly specified?\n")
        sys.exit()
    
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
    '''
    ksnp_list = []
    chr_start_idx = 1000000000
    chr_end_idx = 0
    raw_var_ct = 0
    with open(pvar) as fp:
        header = None
        first_cline = None 
        for line in fp:
            if line[0] == '#':
                header = line
            else:
                first_cline = line
                break
        iid_col = 1
        if (header is not None) and (header[1] == 'C'):
            iid_col = 2
        elem = first_cline.strip().split()
        iid = elem[iid_col]
        if iid in ksnp_set:
            ksnp_list.append(raw_var_ct)
        if elem[0] == args.chr:
            if chr_start_idx > raw_var_ct:
                chr_start_idx = raw_var_ct
            if chr_end_idx < raw_var_ct:
                chr_end_idx = raw_var_ct
        for line in fp:
            elem = line.strip().split()
            iid = elem[iid_col]
            if iid is not None:
                raw_var_ct += 1
            if iid in ksnp_set:
                ksnp_list.append(raw_var_ct)
            if elem[0] == args.chr:
                if chr_start_idx > raw_var_ct:
                    chr_start_idx = raw_var_ct
                if chr_end_idx < raw_var_ct:
                    chr_end_idx = raw_var_ct
    ksnp_ct = len(ksnp_list)
    print("\nIdx range for variants on chr",args.chr,"is",chr_start_idx,"-",chr_end_idx)
    tested_var_ct = chr_end_idx - chr_start_idx + 1
    print("Count of variants to be tested is", tested_var_ct)
    print("Count of model variants to be used is", ksnp_ct)
    
    if block_size > ksnp_ct:
        block_size = int(ksnp_ct/2)*2
        print("\nWarning: Block size in SLEMM was too big and reduced to", block_size)
    
    ksnp_idx = np.array(ksnp_list, np.uint32)
    block_step = int(block_size/2)
    block_ext = int(block_step/2)
    block_ct = int((ksnp_ct-block_step)/block_step)
    
    blocks = []
    for i in range(block_ct):
        blocks.append(ksnp_idx[i*block_step:i*block_step+block_size])
    if (ksnp_ct-block_step) % block_step != 0:
        blocks.append(ksnp_idx[(ksnp_ct-block_size):ksnp_ct])
    block_ct = len(blocks)
    
    block_bound = np.empty([block_ct,2], np.uint32)
    block_bound[0,0] = chr_start_idx
    block_bound[0,1] = blocks[0][block_size-block_ext]
    for i in range(1,block_ct):
        block_bound[i,0] = block_bound[i-1,1]
        block_bound[i,1] = blocks[i][block_size-block_ext]
    block_bound[block_ct-1,1] = chr_end_idx + 1
    print("\nCompleted building blocks:")
    print(block_bound)
    
    # read pgen/bed file
    # compute intermediate terms for assoc tests, xx and xy
    print("\nStarting assoc tests...")
    st_xx = np.empty(tested_var_ct, dtype=np.float32)
    st_xy = np.empty(tested_var_ct, dtype=np.float32)
    with pgenlib.PgenReader(pgen.encode(), raw_sample_ct = raw_sample_ct, sample_subset = sample_subset) as pf:
        buf1 = np.empty([block_size,sample_ct], dtype=np.int8)
        buf2 = np.empty([block_size,sample_ct], dtype=np.float32)
        rvars = np.empty(block_size, dtype=np.float32)
        mat_h = np.empty([block_size,block_size], dtype=np.float32)
        mat_c = np.empty([block_size,block_size], dtype=np.float32)
        int8vec_x = np.empty(sample_ct, dtype=np.int8)
        vec_x = np.empty(sample_ct, dtype=np.float32)
        vec_zx = np.empty(block_size, dtype=np.float32)
        
        half_int8 = np.empty([block_step,sample_ct], dtype=np.int8)
        half = np.empty([block_step,sample_ct], dtype=np.float32)
        half_rvars = np.empty(block_step, dtype=np.float32)
        mat_corner = np.empty([block_step,block_step], dtype=np.float32)
        
        print("There are", block_ct, "blocks.")
        for i in range(block_ct):
            print("  Proccessing block", i)
            if i == 0 or i == block_ct-1:
                pf.read_list(blocks[i], buf1)
                buf2 = buf1.astype('float32')
                buf2_rmeans = buf2.mean(axis=1, keepdims=True)
                rvars = ((2.-buf2_rmeans) * buf2_rmeans/2.)[:,0]
                buf2 -= buf2_rmeans
                if(not is_one_rr):
                    buf2 *= vec_rr
                if(not (is_one_rr and is_const_covar)):
                    buf2 -= np.matmul(np.matmul(buf2,Q), np.transpose(Q))
                mat_h = np.matmul(buf2, np.transpose(buf2))
                
                half = buf2[block_step:,]
                half_rvars = rvars[block_step:]
            else:
                buf2[:block_step,] = half
                rvars[:block_step] = half_rvars
                pf.read_list(blocks[i][block_step:], half_int8)
                half = half_int8.astype('float32')
                half_rmeans = half.mean(axis=1, keepdims=True)
                half_rvars = ((2.-half_rmeans) * half_rmeans/2.)[:,0]
                half -= half_rmeans
                if(not is_one_rr):
                    half *= vec_rr
                if(not (is_one_rr and is_const_covar)):
                    half -= np.matmul(np.matmul(half,Q), np.transpose(Q))
                
                buf2[block_step:,] = half
                rvars[block_step:] = half_rvars
                mat_h[:block_step,:block_step] = mat_corner
                mat_h[:block_step,block_step:] = np.matmul(buf2[:block_step,], np.transpose(half))
                mat_h[block_step:,:block_step] = np.transpose(mat_h[:block_step,block_step:])
                mat_h[block_step:,block_step:] = np.matmul(half, np.transpose(half))
            
            np.copyto(mat_corner, mat_h[block_step:,block_step:])
            mat_h.ravel()[::mat_h.shape[1]+1] += model_snp_ct*tau*rvars
            mat_c, low = linalg.cho_factor(mat_h, overwrite_a=True)
            for j in range(block_bound[i,0], block_bound[i,1]):
                if j % 1000 == 0:
                    print("    Scanning variant", j)
                pf.read(j, int8vec_x)
                vec_x = int8vec_x.astype('float32')
                vec_x_sum = vec_x.sum()
                if vec_x_sum <= args.mac or (2*sample_ct - vec_x_sum) <= args.mac:
                    st_xx[j-chr_start_idx] = np.nan
                    st_xy[j-chr_start_idx] = np.nan
                    continue
                vec_x -= vec_x_sum/sample_ct
                st_xy[j-chr_start_idx] = np.dot(vec_x, vec_y)
                if(not is_one_rr):
                    vec_x *= vec_rr
                if(not (is_one_rr and is_const_covar)):
                    vec_x -= np.dot(Q, np.dot(np.transpose(Q), vec_x))
                vec_zx = np.dot(buf2, vec_x)
                st_xx[j-chr_start_idx] = np.dot(vec_x, vec_x) - np.dot(vec_zx, linalg.cho_solve((mat_c, low), vec_zx))
    # compute statistics for report
    st_xx /= tau
    beta = st_xy/st_xx * cfactor
    serr = np.sqrt((snp_var*cfactor)/st_xx)
    chisq = np.square(beta/serr)
    pval = chi2.sf(chisq,1)
    print("Computed association statistics.")
    
    # write stats to output file
    print("\nWriting stats into", args.out)
    outfp = open(args.out, "w")
    with open(pvar) as fp:
        header = None
        first_cline = None 
        for line in fp:
            if line.strip() is None:
                continue
            if line[0] == '#':
                header = line
            else:
                first_cline = line
                break
        if header is not None:
            header = header.strip() + "\tBETA\tSE\tCHISQ\tP\n"
        elif len(first_cline.split()) == 6:
            header = "CHR\tSNP\tCM\tBP\tA1\tA2\tBETA\tSE\tCHISQ\tP\n"
        else:
            header = "CHR\tSNP\tBP\tA1\tA2\tBETA\tSE\tCHISQ\tP\n"
        outfp.write(header)
        
        var_idx = 0
        if chr_start_idx == 0:
            outfp.write(first_cline.strip() +"\t"+ str(beta[0]) +"\t"+ str(serr[0]) +"\t"+ str(chisq[0]) +"\t"+ str(pval[0]) +"\n")
        for line in fp:
            if line.strip() is not None:
                var_idx += 1
                if var_idx > chr_end_idx:
                    break
                if var_idx >= chr_start_idx:
                    list_idx = var_idx - chr_start_idx
                    outfp.write(line.strip() +"\t"+ str(beta[list_idx]) +"\t"+ str(serr[list_idx]) +"\t"+ str(chisq[list_idx]) +"\t"+ str(pval[list_idx]) +"\n")
    outfp.close()
    print(int(time.time()-start), "s taken for the whole analysis.\n")

