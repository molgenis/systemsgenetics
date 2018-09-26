import numpy as np
import argparse
import random


parser = argparse.ArgumentParser(description='Simulate gene expression levels using expression ~ cc1 + cc2 + snp:cc1 + snp:cc2')
parser.add_argument('cellcount_file', help='file containing cell counts')
parser.add_argument('genotype_file', help='file containing genotypes')
parser.add_argument('out_dir', help='output directory to write the simulated data to')

args = parser.parse_args()


# Use actual (or Decon-Cell predicted) cell counts to simulate the expression levels. Parse the file
# File should be in format of
#          CC1    CC2
# sample1  75     14
# sample2  84     4
def simulate_cellcounts(number_of_samples):
    cellcount_names = []
    samples = []
    cellcount_per_sample = {}
    with open(args.cellcount_file) as input_file:
        cellcount_names = {}
        cellcount_names = input_file.readline().strip().split('\t')
        for line in input_file:
            line = line.strip().split('\t')
            samples.append(line[0])
            cellcount_per_sample[line[0]] = line[1:]
    
    samples_tmp = list(samples)
    random.shuffle(samples_tmp)
    random_selected_samples_list = samples_tmp[:number_of_samples]
    random_selected_samples_set = set(random_selected_samples_list)
    
    with open(args.genotype_file) as input_file, open(args.out_dir+'/expression/simulated_expression_'+str(number_of_samples)+'samples.txt','w') as (
            out_simulatedExpression), open(args.out_dir+'/snpsToTest.txt','w') as (
            out_snpToTest), open(args.out_dir+'/betas/beta_info_cc'+str(number_of_samples)+'samples.txt','w') as (
            out_beta_info), open(args.out_dir+'/genotypes/genotypes_'+str(number_of_samples)+'samples.txt','w') as (
            out_genotype),open(args.out_dir+'/cellcounts/cellcounts_'+str(number_of_samples)+'samples.txt','w') as out_cellcount:
        out_snpToTest.write('gene\tsnp\n')
        beta_header = 'gene\tsnp'
        out_beta_info.write('gene\tsnp')
            
        for cc in cellcount_names:
            out_beta_info.write('\t'+cc+'_beta\t'+cc+':GT_beta')
            beta_header += '\t'+cc+'_beta\t'+cc+':GT_beta'
            out_cellcount.write('\t'+cc)
        out_cellcount.write('\n')
        out_beta_info.write('\terror\n')
        
        for sample in random_selected_samples_list:
            cellcounts = cellcount_per_sample[sample]
            out_cellcount.write(sample)
            for cellcount in cellcounts:
                out_cellcount.write('\t'+cellcount)
            out_cellcount.write('\n')

        
        # read in all lines of the file so we can check how many lines it has for sampling
        genotype_lines = input_file.read().split('\n')
        header = genotype_lines[0]
        out_simulatedExpression.write('\t'+'\t'.join(random_selected_samples_list)+'\n')
        out_genotype.write('\t'+'\t'.join(random_selected_samples_list)+'\n')
        header = header.strip().split('\t')
        if not header == samples:
            for index, sample in enumerate(header):
                print(samples[index], sample)
            raise RuntimeError("header and samples not same order")
        
        for index, line in enumerate(genotype_lines[1:]):
            if index % 100 == 0:
                print('processed',index,'lines')
            line = line.strip().split('\t')
            snp = line[0]
            if len(snp.strip()) == 0:
                continue
            out_snpToTest.write('gene_'+str(index)+'\t'+snp+'\n')
            out_simulatedExpression.write('gene_'+str(index))
            # mean and standard deviation, draw from normal distribution for betas
            mu, sigma = 3, 5
            betas = np.random.normal(mu, sigma, len(cellcount_names)*2)
            # add some error
            error = np.random.randint(4,5)
            out_genotype.write(snp)
            out_beta_info.write('gene_'+str(index)+'\t'+snp)
            for cc_index, cellount in enumerate(cellcount_names):
                out_beta_info.write('\t'+str(betas[cc_index])+'\t'+str(betas[cc_index+len(cellcount_names)]))
            out_beta_info.write('\t'+str(error)+'\n')
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
                    cc_contribution = betas[cc_index] * float(cellcount)
                    expression += cc_contribution
                    cc_snp_contribution = betas[cc_index+len(cellcounts)] * float(cellcount) * float(dosage)
                    expression += cc_snp_contribution
     
                out_simulatedExpression.write('\t'+str(expression+error))
            out_simulatedExpression.write('\n')
            out_genotype.write('\n')
    print('output written to '+args.out_dir+'/')

with open(args.cellcount_file) as input_file:
    n_samples = -1
    for line in input_file:
        n_samples += 1
        
for i in range(2600, n_samples, 100):
    simulate_cellcounts(i)
