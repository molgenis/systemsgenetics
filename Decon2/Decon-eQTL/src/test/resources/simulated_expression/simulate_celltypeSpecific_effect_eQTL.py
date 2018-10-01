import numpy as np
import argparse
import random


parser = argparse.ArgumentParser(description='Simulate gene expression levels using expression ~ cc1 + cc2 + snp:cc1 + snp:cc2')
parser.add_argument('cellcount_file', help='file containing cell counts')
parser.add_argument('genotype_file', help='file containing genotypes')
parser.add_argument('out_dir', help='output directory to write the simulated data to')
parser.add_argument('number_of_samples', help='Number of samples to simulate', type=int)
parser.add_argument('batch', help='Name of the batch')

args = parser.parse_args()


# Use actual (or Decon-Cell predicted) cell counts to simulate the expression levels. Parse the file
# File should be in format of
#          CC1    CC2
# sample1  75     14
# sample2  84     4
def simulate_cellcounts(number_of_samples, batch):
    cellcount_names = []
    samples = []
    cellcount_per_sample = {}
    with open(args.cellcount_file) as input_file:
        cellcount_names = input_file.readline().strip().split('\t')
        
        for line in input_file:
            line = line.strip().split('\t')
            samples.append(line[0])
            cellcount_per_sample[line[0]] = line[1:]

    sample_info = {}    
    samples_tmp = list(samples)

    random.shuffle(samples_tmp)
    random_selected_samples_list = samples_tmp[:number_of_samples]
    random_selected_samples_set = set(random_selected_samples_list)
    batch_name = 'batch'+str(batch)+'.'+str(number_of_samples)
    
    # opening outfiles
    out_beta_info = open(args.out_dir+'/betas/beta_info_cc'+batch_name+'samples.txt','w')
    out_genotype = open(args.out_dir+'/genotypes/genotypes_'+batch_name+'samples.txt','w')
    out_cellcount = open(args.out_dir+'/cellcounts/cellcounts_'+batch_name+'samples.txt','w')
    out_simulatedExpression =  open(args.out_dir+'/expression/simulated_expression_'+batch_name+'samples.txt','w')
    out_snpToTest = open(args.out_dir+'/snpsToTest.txt','w') 
            
    out_beta_info.write('gene\tsnp')
    for cc in cellcount_names:
        out_beta_info.write('\t'+cc+'_beta\t'+cc+':GT_beta')
        out_cellcount.write('\t'+cc)
    out_beta_info.write('\terror\n')
    out_cellcount.write('\n')
            
               #sample_info[batch_name] = [random_selected_samples_list, random_selected_samples_set,
        #                                         out_beta_info, out_genotype, out_cellcount, out_simulatedExpression]     
    with open(args.genotype_file) as input_file:
        # read in all lines of the file so that it can be closed after
        genotype_lines = input_file.read().split('\n')
        genptype_header = genotype_lines[0].strip().split('\t')
        if not genptype_header == samples:
            for index, sample in enumerate(genptype_header):
                print(samples[index], sample)
            raise RuntimeError("header and samples not same order")
    for sample in random_selected_samples_list:
        cellcounts = cellcount_per_sample[sample]
        out_cellcount.write(sample)
        for cellcount in cellcounts:
            out_cellcount.write('\t'+cellcount)
        out_cellcount.write('\n')

    out_simulatedExpression.write('\t'+'\t'.join(random_selected_samples_list)+'\n')
    out_genotype.write('\t'+'\t'.join(random_selected_samples_list)+'\n')

    mu, sigma = 5, 7
    print('simulate betas')
    error = np.random.normal(0,0, len(genotype_lines[1:])+2)
    betas = np.random.normal(mu, sigma, 10000)
    print('done')        
    out_snpToTest.write('gene\tsnp\n')
    # use random genotypes
    
    genotype_lines = genotype_lines[1:]
    for index, line in enumerate(genotype_lines):
        if index % 100 == 0:
            print('processed',index,'lines')
            if index == 500:
                break
        line = line.strip().split('\t')
        snp = line[0]
        if len(snp.strip()) == 0:
            continue
        out_snpToTest.write('gene_'+str(index)+'\t'+snp+'\n')

        out_simulatedExpression.write('gene_'+str(index))
        out_genotype.write(snp)
        out_beta_info.write('gene_'+str(index)+'\t'+snp)
        
        current_cc_betas = betas[random.sample(range(10000), len(cellcount_names))]
        current_cc_gt_betas = betas[random.sample(range(10000), len(cellcount_names))]
        for cc_index, cellcount_name in enumerate(cellcount_names):
            out_beta_info.write('\t'+str(current_cc_betas[cc_index])+'\t'+str(current_cc_gt_betas[cc_index]))

        #error_index =  random.randint(0,10000
        #out_beta_info.write('\t'+str(betas[error_index\)+'\n')
        
        out_beta_info.write('\t'+str(0)+'\n')
        for sample_index, dosage in enumerate(line[1:]):
            sample = samples[sample_index]
            if sample not in random_selected_samples_set:
                continue
            dosage = float(dosage)
            cellcounts = cellcount_per_sample[sample]
            out_genotype.write('\t'+str(dosage))
            # Expression will be made with expression = cc1 + cc2 + snp:cc1 + snp:cc2 + error
            # so start with 0
            expression = 0
            # then add cc and cc*snp. for cc*snp can add the beta
            for cc_index, cellcount in enumerate(cellcounts):
                cc_contribution = current_cc_betas[cc_index] * float(cellcount)
                expression += cc_contribution
                cc_snp_contribution = current_cc_gt_betas[cc_index] * float(cellcount) * float(dosage)
                expression += cc_snp_contribution
            #out_simulatedExpression.write('\t'+str(expression+error[error_index]))
            out_simulatedExpression.write('\t'+str(expression+0))

        out_simulatedExpression.write('\n')
        out_genotype.write('\n')
            
    print('output written to '+args.out_dir+'/')

    out_beta_info.close()
    out_genotype.close()
    out_cellcount.close()
    out_simulatedExpression.close()
    out_snpToTest.close()
        
simulate_cellcounts(args.number_of_samples,args.batch)
