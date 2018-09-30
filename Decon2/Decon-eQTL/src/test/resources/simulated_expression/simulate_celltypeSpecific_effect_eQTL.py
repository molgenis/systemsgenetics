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
def simulate_cellcounts(number_of_samples, n_snps=1000):
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

    sample_info = {}    
    samples_tmp = list(samples)
    for n_sample in number_of_samples:
        # repeat sampling 10 times per n_sample
        for x in range(1,11):                
            random.shuffle(samples_tmp)
            random_selected_samples_list = samples_tmp[:n_sample]
            random_selected_samples_set = set(random_selected_samples_list)
            batch_name = 'batch'+str(x)+'.'+str(n_sample)
            out_beta_info = open(args.out_dir+'/betas/beta_info_cc'+batch_name+'samples.txt','w')
            out_genotype = open(args.out_dir+'/genotypes/genotypes_'+batch_name+'samples.txt','w')
            out_cellcount = open(args.out_dir+'/cellcounts/cellcounts_'+batch_name+'samples.txt','w')
            
            out_simulatedExpression =  open(args.out_dir+'/expression/simulated_expression_'+batch_name+'samples.txt','w')
            
            out_beta_info.write('gene\tsnp')
            for cc in cellcount_names:
                out_beta_info.write('\t'+cc+'_beta\t'+cc+':GT_beta')
                out_cellcount.write('\t'+cc)
            out_beta_info.write('\terror\n')
            out_cellcount.write('\n')
            
            
            sample_info[batch_name] = [random_selected_samples_list, random_selected_samples_set,
                                                 out_beta_info, out_genotype, out_cellcount, out_simulatedExpression]
        
    
    with open(args.genotype_file) as input_file, open(args.out_dir+'/snpsToTest.txt','w') as out_snpToTest:
                    
        for batch_name in sample_info:
            for sample in sample_info[batch_name][0]:
                cellcounts = cellcount_per_sample[sample]
                sample_info[batch_name][4].write(sample)
                for cellcount in cellcounts:
                    sample_info[batch_name][4].write('\t'+cellcount)
                sample_info[batch_name][4].write('\n')

        # read in all lines of the file so we can check how many lines it has for sampling
        genotype_lines = input_file.read().split('\n')
        header = genotype_lines[0]
        for batch_name in sample_info: 
            s_data = sample_info[batch_name]
            s_data[5].write('\t'+'\t'.join(s_data[0])+'\n')
            s_data[3].write('\t'+'\t'.join(s_data[0])+'\n')
        header = header.strip().split('\t')
        if not header == samples:
            for index, sample in enumerate(header):
                print(samples[index], sample)
            raise RuntimeError("header and samples not same order")
        mu, sigma = 3, 5
        # add some error
        print('simulate betas')
        error = np.random.normal(0,0, len(genotype_lines[1:])+2)
        betas = np.random.normal(mu, sigma, len(cellcount_names)*2*(len(genotype_lines[1:])+2))
        beta_index = 0
        print('done')
        
        out_snpToTest.write('gene\tsnp\n')
        # use random genotypes
        random.shuffle(genotype_lines)
        for index, line in enumerate(genotype_lines[1:]):
            if index >= n_snps:
                break
            beta_index += 1
            if index % 10 == 0:
                print('processed',index,'lines')
            line = line.strip().split('\t')
            snp = line[0]
            if len(snp.strip()) == 0:
                continue
            out_snpToTest.write('gene_'+str(index)+'\t'+snp+'\n')

            
            for n_samples in sample_info:
                s_data = sample_info[n_samples]
                # mean and standard deviation, draw from normal distribution for betas
                
                
            
                s_data[5].write('gene_'+str(index))
                s_data[3].write(snp)
                s_data[2].write('gene_'+str(index)+'\t'+snp)
            
                for cc_index, cellcount_name in enumerate(cellcount_names):
                    s_data[2].write('\t'+str(betas[cc_index]*beta_index)+'\t'+str(betas[cc_index+len(cellcount_names)*beta_index]))

                s_data[2].write('\t'+str(error[beta_index])+'\n')

            for sample_index, dosage in enumerate(line[1:]):
                sample = samples[sample_index]
                
                for n_samples in sample_info:
                    s_data = sample_info[n_samples]
                    if sample not in s_data[1]:
                        continue
                    
                    dosage = float(dosage)
                    cellcounts = cellcount_per_sample[sample]
                    s_data[3].write('\t'+str(dosage))
                    # Expression will be made with expression = cc1 + cc2 + snp:cc1 + snp:cc2 + error
                    # so start with 0
                    expression = 0
                    # then add cc and cc*snp. for cc*snp can add the beta
                    for cc_index, cellcount in enumerate(cellcounts):
                        cc_contribution = betas[cc_index*beta_index] * float(cellcount)
                        expression += cc_contribution
                        cc_snp_contribution = betas[cc_index+len(cellcounts)*beta_index] * float(cellcount) * float(dosage)
                        expression += cc_snp_contribution
         
                    s_data[5].write('\t'+str(expression+error[beta_index]))
            for n_samples in sample_info:
                s_data = sample_info[n_samples]
                s_data[5].write('\n')
                s_data[3].write('\n')
            
    print('output written to '+args.out_dir+'/')

    for batch_name in sample_info:
        sample_info[batch_name][2].close()
        sample_info[batch_name][3].close()
        sample_info[batch_name][4].close()
        sample_info[batch_name][5].close()
        
with open(args.cellcount_file) as input_file:
    n_samples = -1
    for line in input_file:
        n_samples += 1
    
number_of_samples = []    
for i in range(100, n_samples, 100):
    number_of_samples.append(i)
simulate_cellcounts(number_of_samples,n_snps=3)
