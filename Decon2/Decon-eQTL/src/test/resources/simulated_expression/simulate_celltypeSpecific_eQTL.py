import numpy as np
import argparse
import random
import sys
import os
import scipy.stats
import collections
import logging 

parser = argparse.ArgumentParser(description='Simulate gene expression levels using expression ~ cc1 + cc2 + snp:cc1 + snp:cc2')
parser.add_argument('cellcount_file', help='file containing cell counts')
parser.add_argument('genotype_file', help='file containing genotypes')
parser.add_argument('out_dir', help='output directory to write the simulated data to')
parser.add_argument('number_of_snps', help='Number of snps to simulate', type=int)
parser.add_argument('batch', help='Name of the batch')

args = parser.parse_args()

if not os.path.exists(args.out_dir+'/betas/'):
    os.makedirs(args.out_dir+'/betas/')
if not os.path.exists(args.out_dir+'/genotypes/'):
    os.makedirs(args.out_dir+'/genotypes/')
if not os.path.exists(args.out_dir+'/cellcounts/'):
    os.makedirs(args.out_dir+'/cellcounts/')
if not os.path.exists(args.out_dir+'/expression/'):
    os.makedirs(args.out_dir+'/expression/')
if not os.path.exists(args.out_dir+'/snpsToTest/'):
    os.makedirs(args.out_dir+'/snpsToTest/')

format = '%(asctime)s - %(levelname)s - %(funcName)s\t- %(message)s'
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format=format, datefmt="%Y-%m-%d %H:%M")

class Expression_simulator():
    '''Simulate expression from multiple cell types using real cell count and genotype data'''
    def __init__(self, genotype_file, out_dir, number_of_snps, batch):
        self.number_of_snps = number_of_snps
        self.batch_name = 'batch'+str(batch)
        self.genotype_file = genotype_file
        self.samples = []
        self.out_dir = out_dir
        self.cellcount_per_sample = {}
        self.cellcount_names = None
        self.number_of_celltypes = None
        self.cc_beta = {}
        self.cc_gt_beta = {}
        self.error = {}
        
        # want to keep the order of SNPs the same, so use ordered dict
        self.genotypes = collections.OrderedDict()
        
        
        # read the genotype and cell count files 
        self.__read_cellcount_data()
        self.__read_genotype_data()
        self.__write_cellcounts()
        self.__write_genotypes()
        self.__write_snps_to_test()
        self.__make_betas()
        self.__make_errors()
        self.__write_betas()
        
    def __read_cellcount_data(self):
        '''Read the cellcount data'''
        logging.info('Read cell count data')
        cellcount_names = []
        with open(args.cellcount_file) as input_file:
            self.cellcount_names = input_file.readline().strip().split('\t')
            for line in input_file:
                line = line.strip().split('\t')
                self.samples.append(line[0])
                self.cellcount_per_sample[line[0]] = [float(x) for x in line[1:]]
                self.number_of_celltypes = len(line[1:])
                
    def __read_genotype_data(self):
        '''Read the genotype data'''
        logging.info('Read genotype data')
        with open(self.genotype_file) as input_file:
            genotype_header = input_file.readline().strip().split('\t')
            # Check that the samples from the genotype file are the same as the samples in the cell counts file
            if not genotype_header == self.samples:
                for index, sample in enumerate(genotype_header):
                    if self.samples[index] != sample:
                        logging.info(index, self.samples[index], sample)
                raise RuntimeError("header and samples not same order")
            
            # Read in all lines after header so that SNPs can be shuffled and sub sampled 
            lines = input_file.read().split('\n')
            random.shuffle(lines)
            
            for line in lines[:self.number_of_snps]:
                # would be faster to now also write the genotypes and do the expression simulation
                # but want to keep these separate for readability. 
                line = line.strip().split('\t')
                # self.genotypes is an ordered dict, so the order of SNPs is preserved and can be used later
                self.genotypes[line[0]] = line[1:]
            
    def __write_betas(self):
        outfile = self.out_dir+'/betas/betas_'+self.batch_name+'.txt'
        with open(outfile,'w') as out:
            out.write('gene\tsnp')
            for cc in self.cellcount_names:
                out.write('\t'+cc+'_beta\t'+cc+':GT_beta')
            out.write('\terror\n')
        
            
            for index, snp in enumerate(self.genotypes):            
                out.write('gene_'+str(index)+'\t'+snp)
                for cc_index in range(0, self.number_of_celltypes) :
                    out.write('\t'+str(self.cc_beta[snp][cc_index])+'\t'+str(self.cc_gt_beta[snp][cc_index]))
                out.write('\t'+str(self.error[snp])+'\n')
   
                
        logging.info('Betas written to '+outfile)
                            
    def __write_cellcounts(self):
        outfile = self.out_dir+'/cellcounts/cellcounts_'+self.batch_name+'.txt'
        with open(outfile,'w') as out:
            for cc in self.cellcount_names:
                out.write('\t'+str(cc))
            out.write('\n')

        logging.info('Cellcounts written to '+outfile)

    def __write_snps_to_test(self):
        outfile = self.out_dir+'snpsToTest/snpsToTest_'+self.batch_name+'.txt'
        with open(outfile,'w') as out:
            out.write('gene\tsnp\n')

            for index, snp in enumerate(self.genotypes):
                # Since gene just gets simulated and we use the same order of SNPs, gene name can just be the index
                out.write('gene_'+str(index)+'\t'+snp+'\n')
        logging.info('Snps to test written to '+outfile)
        
    def __write_genotypes(self):
        outfile = self.out_dir+'/genotypes/genotypes_'+self.batch_name+'.txt'
        with open(outfile,'w') as out:
            out.write('\t'+'\t'.join(self.samples)+'\n')

            # since genotypes is an ordered dict we can just loop over keys
            for snp in self.genotypes:
                out.write(snp)
                # When reading in genotype file we check that the samples are in the same order as
                # the samples in the cell counts file, so can just write out now
                for dosage in self.genotypes[snp]:
                    out.write('\t'+dosage)
                out.write('\n')
        logging.info('Genotypes written to '+outfile)

    def __make_errors(self):
        mu, sigma = 0, 250
        error = np.random.normal(mu, sigma, self.number_of_snps)
        for index, snp in enumerate(self.genotypes):
            self.error[snp] = error[index]

    def __make_betas(self):
        '''Get betas for all the different situations. 
        Want to simulate different combinations of betas because some cell counts can be correlated to each other, which
        can have different effects on beta prediciton based on which combination of cc have the same effect
        Per cell type there are 2 betas to simulate: the cc:gt interaction term beta and the cc term beta. For these
            beta terms there is the restriction that beta*cc:gt + 2*beta*cc >= 0
        
        We want to simulate:
            1. All betas strong positive
            2. All interaction term betas strong negative and cc term beta strong positive
                                                           (note: simualted betas direction towards the allele that is
                                                           coded as 2. cc term betas can not be negative as there is
                                                           no direction, and contribution should be positive)
                                                           Also, cc term need to be strong positive to fulfill 
                                                           (interaction beta + 2*cc beta) >= 0 
            3. For each cell type cc * gt beta positive and cc beta positive, while all others 0
            4. For each cell type cc * gt beta positive and 0, while all others 0
            5+6. Like 3 and 4 but with negative betas
            7. Combinations of cell types positive/negative (e.g. cc1*gt and cc2*gt positive, rest 0)
            8. Combinations of different opposite effects
            9. All of them mixed
            
            Each of these want to do with different range of errors'''
        logging.info('simulate betas')
        number_of_samples = len(self.samples)
        # 1. all betas postive. Want to have different divisions of strong/weak positive effects.
        mu, sigma = 15, 5 # mean and standard deviation
        number_of_betas_to_simulate = self.number_of_celltypes*self.number_of_snps        
        # simulate as many betas in one go as it is faster
        cc_betas = np.random.normal(mu, sigma, number_of_betas_to_simulate)
        cc_gt_betas = np.random.normal(mu, sigma, number_of_betas_to_simulate)
        
        # then to make it more understandable and to test that correct values got simulated
        # separate them for different cell types
        for index, snp in enumerate(self.genotypes):
            # Every n_cell type window of the beta list is for one sample, e.g. with 3 cell types:
            # beta = [  3, 4, 4,    3, 4, 3,    2, 4, 5]
            #         --sample1-- --sample2-- --sample3--
            # So calculate the index using the sample index and the number of cell types
            start_index = index * self.number_of_celltypes
            end_index = start_index + self.number_of_celltypes
            self.cc_beta[snp] = cc_betas[start_index:end_index]
            self.cc_gt_beta[snp] = cc_gt_betas[start_index:end_index]
            # make sure that for B1*cc + B2*cc*gt -> B1 + (2*B2) >= 0
            for beta_index, cc_beta in enumerate(self.cc_beta[snp]):
                cc_gt_beta = self.cc_gt_beta[snp][beta_index]
                
                combined_betas = (2*cc_beta)+cc_gt_beta
                if combined_betas < 0:
                    logging.info('cc_beta: '+str(cc_beta)+', cc_gt_beta: '+str(cc_gt_beta)+', 2*cc_beta + cc_gt_beta = '+str(combined_betas))
                    raise RuntimeError('2*cc_beta ('+str(2*cc_beta)+') + cc_gt_beta ('+str(cc_gt_beta)+') should be >= 0, but was: '+str(combined_betas))
    
    # Use actual (or Decon-Cell predicted) cell counts to simulate the expression levels. Parse the file
    # File should be in format of
    #          CC1    CC2
    # sample1  75     14
    # sample2  84     4
    def simulate_expression(self):
        '''Simulate expression by summing over the simulated expression per cell type + error'''
        outfile = self.out_dir+'/expression/simulated_expression_'+self.batch_name+'.txt'
        with open(outfile,'w') as out:
            out.write('\t'+'\t'.join(self.samples)+'\n')

            for index, snp in enumerate(self.genotypes):
                out.write('gene_'+str(index))
                if index % 100 == 0:
                    logging.info('processed '+str(index)+' SNPs')
                
                # When reading the genotype and cell count data we checked that the header is in same order, so an just use the index
                for sample_index, sample in enumerate(self.samples):
                    dosage = self.genotypes[snp][sample_index]
                    cellcounts = self.cellcount_per_sample[sample]
                    # Expression will be made with expression = cc1 + cc2 + snp:cc1 + snp:cc2 + error
                    # so start with 0
                    expression = 0
                    # then add cc and cc*snp. for cc*snp can add the beta
                    for cc_index, cellcount in enumerate(cellcounts):
                        cc_contribution = self.cc_beta[snp][cc_index] * cellcount
                        expression += cc_contribution
                        cc_snp_contribution = self.cc_gt_beta[snp][cc_index] * float(cellcount) * float(dosage)
                        expression += cc_snp_contribution
     
                    out.write('\t'+str(expression+self.error[snp]))
                out.write('\n')
        logging.info('Simulated expression written to '+outfile)
               
    


if __name__ == '__main__':
    expression_simulator = Expression_simulator(args.genotype_file, args.out_dir, args.number_of_snps, args.batch)
    expression_simulator.simulate_expression()