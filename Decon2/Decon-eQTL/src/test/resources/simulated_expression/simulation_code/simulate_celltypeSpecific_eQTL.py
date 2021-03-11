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
        self.cc_beta_per_type = {'all_positive':{}, 'all_negative':{}}
        self.cc_gt_beta_per_type = {'all_positive':{}, 'all_negative':{}}
        self.error = {}
        self.error_size = [1,5,10,20,40,80,160,250,1000]
        self.beta_types = ['all_positive', 'all_negative','mixed','some_same']
        # want to keep the order of SNPs the same, so use ordered dict
        self.genotypes = collections.OrderedDict()
        
        
        # read the genotype and cell count files 
        self.__read_cellcount_data()
        self.__read_genotype_data()
        self.__make_betas()
        self.__make_errors()
        self.__simulate_expression()
        self.__write_betas()
        self.__write_cellcounts()
        self.__write_genotypes()
        self.__write_snps_to_test()
        
    def __read_cellcount_data(self):
        '''Read the cellcount data'''
        logging.info('Read cell count data')
        cellcount_names = []
        with open(args.cellcount_file) as input_file:
            self.cellcount_names = input_file.readline().strip().split('\t')
            for line in input_file:
                line = line.strip().split('\t')
                
                # skip negative cell counts
                all_cc = [float(x) for x in line[1:]]
                contains_negative = False
                for index, cc in enumerate(all_cc):
                    if cc < 0:
                        contains_negative = True
                        all_cc[index] = abs(all_cc[index])
                if contains_negative:
                    logging.info(line[0]+' contains negative cellcount (changing it to abs(value)): '+', '.join(line[1:]))
                
                self.samples.append(line[0])
                self.cellcount_per_sample[line[0]] = all_cc
                self.number_of_celltypes = len(all_cc)
                
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
            lines = input_file.read().rstrip().split('\n')
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
            out.write('\terror\tbeta_type\terror_sigma\n')
            
            for beta_type in self.beta_types:
                for sigma in self.error_size:
                    for index, snp in enumerate(self.genotypes):
                        gene = 'gene_'+str(index)+'_'+beta_type+'_'+str(sigma)
                        out.write(gene+'\t'+snp)
                        for cc_index in range(0, self.number_of_celltypes):
                            out.write('\t'+str(self.cc_beta_per_type[beta_type][snp][cc_index])+'\t'+str(self.cc_gt_beta_per_type[beta_type][snp][cc_index]))
                        out.write('\t'+str(self.error[sigma][snp])+'\t'+beta_type+'\t'+str(sigma)+'\n')
   
        logging.info('Betas written to '+outfile)
                            
    def __write_cellcounts(self):
        outfile = self.out_dir+'/cellcounts/cellcounts_'+self.batch_name+'.txt'
        with open(outfile,'w') as out:
            for cc in self.cellcount_names:
                out.write('\t'+str(cc))
            out.write('\n')
            for sample in self.samples:
                out.write(sample)
                for i in range(0, self.number_of_celltypes):
                    out.write('\t'+str(self.cellcount_per_sample[sample][i]))
                out.write('\n')

        logging.info('Cellcounts written to '+outfile)

    def __write_snps_to_test(self):
        outfile = self.out_dir+'snpsToTest/snpsToTest_'+self.batch_name+'.txt'
        with open(outfile,'w') as out:
            out.write('gene\tsnp\n')

            for beta_type in self.beta_types:
                for error_type in self.error_size:
                    for index, snp in enumerate(self.genotypes):
                        gene = 'gene_'+str(index)+'_'+beta_type+'_'+str(error_type)
                        # Since gene just gets simulated and we use the same order of SNPs, gene name can just be the index
                        out.write(gene+'\t'+snp+'\n')
        logging.info('Snps to test written to '+outfile)
        
    def __write_genotypes(self):
        outfile = self.out_dir+'/genotypes/genotypes_'+self.batch_name+'.txt'
        permuted_outfile = self.out_dir+'/genotypesPermuted/genotypes_'+self.batch_name+'.txt'
        with open(outfile,'w') as out, open(permuted_outfile,'w') as out_permuted:
            out.write('\t'+'\t'.join(self.samples)+'\n')
            
            # also write permuted genotypes
            shuffled_samples = copy.deepcopy(self.samples)
            random.shuffle(shuffled_samples)
            out_permuted.write('\t'+'\t'.join(self.samples)+'\n')
            
            
            # since genotypes is an ordered dict we can just loop over keys
            for snp in self.genotypes:
                out.write(snp)
                out_permuted.write(snp)
                # When reading in genotype file we check that the samples are in the same order as
                # the samples in the cell counts file, so can just write out now
                for dosage in self.genotypes[snp]:
                    out.write('\t'+dosage)
                    out_permuted.write('\t'+dosage)
                out.write('\n')
                out_permuted.write('\n')
        logging.info('Genotypes written to '+outfile)
        logging.info('Permuted genotypes written to '+permuted_outfile)

    def __make_errors(self):
        mu = 0
        for sigma in self.error_size:
            error = np.random.normal(mu, sigma, self.number_of_snps)
            self.error[sigma] = {}
            for index, snp in enumerate(self.genotypes):
                self.error[sigma][snp] = error[index]

    def __gaussian_betas(self, mu, sigma):
        '''Simulate betas where all betas follow the same distribution (but for cc betas always positive)'''
        number_of_betas_to_simulate = self.number_of_celltypes*self.number_of_snps        
        # simulate as many betas in one go as it is faster
        cc_betas = np.random.normal(abs(mu), sigma, number_of_betas_to_simulate)
        cc_betas = [abs(x) for x in cc_betas]
        cc_gt_betas = np.random.normal(mu, sigma, number_of_betas_to_simulate)
        
        cc_beta = {}
        cc_gt_beta = {}
        
        # then to make it more understandable and to test that correct values got simulated
        # separate them for different cell types
        for index, snp in enumerate(self.genotypes):
            # Every n_cell type window of the beta list is for one sample, e.g. with 3 cell types:
            # beta = [  3, 4, 4,    3, 4, 3,    2, 4, 5]
            #         --sample1-- --sample2-- --sample3--
            # So calculate the index using the sample index and the number of cell types
            start_index = index * self.number_of_celltypes
            end_index = start_index + self.number_of_celltypes
            cc_beta[snp] = cc_betas[start_index:end_index]
            cc_gt_beta[snp] = cc_gt_betas[start_index:end_index]
            
            for beta_index, cc_beta_value in enumerate(cc_beta[snp]):
                # make sure that for B1*cc + B2*cc*gt -> B1 + (2*B2) >= 0
                combined_betas = cc_beta[snp][beta_index] + (2*cc_gt_beta[snp][beta_index])
                x = 0
                while combined_betas < 0:
                    x+=1
                    if x == 20:
                        exit()
                    cc_beta[snp][beta_index] = cc_beta[snp][beta_index]*2
                    combined_betas = cc_beta[snp][beta_index] + (2*cc_gt_beta[snp][beta_index])     
        return(cc_beta, cc_gt_beta)
    
    def __same_betas(self, mu, sigma):
        '''Simulate betas where all betas follow the same distribution (but for cc betas always positive)'''
        number_of_betas_to_simulate = self.number_of_celltypes*self.number_of_snps        
        # simulate as many betas in one go as it is faster. Copy some of the betas so that some cell types get the same gt beta
        cc_betas = np.random.normal(abs(mu), sigma, number_of_betas_to_simulate)
        cc_betas = [abs(x) for x in cc_betas]
        cc_gt_betas = np.random.normal(mu, sigma, number_of_betas_to_simulate)        
        
        cc_beta = {}
        cc_gt_beta = {}
        
        # then to make it more understandable and to test that correct values got simulated
        # separate them for different cell types
        for index, snp in enumerate(self.genotypes):
            # Every n_cell type window of the beta list is for one sample, e.g. with 3 cell types:
            # beta = [  3, 4, 4,    3, 4, 3,    2, 4, 5]
            #         --sample1-- --sample2-- --sample3--
            # So calculate the index using the sample index and the number of cell types
            start_index = index * self.number_of_celltypes
            end_index = start_index + self.number_of_celltypes
            
            # make 2 cell types have the same beta by select random index of value to change and one to change it with
            # this makes it so that e.g. the gt beta of cell type 1 is same as cell type 4 (randomly selected)
            # take -1 because end_index isi closed form 
            to_change = random.randint(start_index, end_index-1)
            change_with = random.randint(start_index, end_index-1)
            while to_change == change_with:
                change_with = random.randint(start_index, end_index-1)
            try:
                cc_gt_betas[to_change] = cc_gt_betas[change_with]
            except IndexError:
                logging.info('to_change: '+str(to_change))
                logging.info('change_with: '+str(change_with))
                logging.info('len(cc_gt_betas): '+str(len(cc_gt_betas)))
                logging.info('start_index: '+str(start_index))
                logging.info('end_index: '+str(end_index))
                raise

            cc_beta[snp] = cc_betas[start_index:end_index]
            cc_gt_beta[snp] = cc_gt_betas[start_index:end_index]
            for beta_index, cc_beta_value in enumerate(cc_beta[snp]):
                # make sure that for B1*cc + B2*cc*gt -> B1 + (2*B2) >= 0
                combined_betas = cc_beta[snp][beta_index] + (2*cc_gt_beta[snp][beta_index])
                while combined_betas < 0:
                    cc_beta[snp][beta_index] = cc_beta[snp][beta_index]*2
                    combined_betas = cc_beta[snp][beta_index] + (2*cc_gt_beta[snp][beta_index])     
        return(cc_beta, cc_gt_beta)
    
    
           
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
                                                           (2*interaction beta + cc beta) >= 0 
            3. mix of effects positive and negative effects
            4. same effects
            Each of these want to do with different range of errors'''
        logging.info('simulate betas')
        number_of_samples = len(self.samples)
        # 1. all betas postive. 
        cc_beta, cc_gt_beta = self.__gaussian_betas(15, 5)
        self.cc_beta_per_type['all_positive'] = cc_beta
        self.cc_gt_beta_per_type['all_positive'] = cc_gt_beta
        # 2. all cc_gt betas negative
        cc_beta, cc_gt_beta = self.__gaussian_betas(-15, 5)
        self.cc_beta_per_type['all_negative'] = cc_beta
        self.cc_gt_beta_per_type['all_negative'] = cc_gt_beta
        # 3. Mix of effects
        cc_beta, cc_gt_beta = self.__gaussian_betas(0, 20)
        self.cc_beta_per_type['mixed'] = cc_beta
        self.cc_gt_beta_per_type['mixed'] = cc_gt_beta        
        # 4. Same effects, make sure the effects are the same for some cell types
        cc_beta, cc_gt_beta = self.__same_betas(0, 20)
        self.cc_beta_per_type['some_same'] = cc_beta
        self.cc_gt_beta_per_type['some_same'] = cc_gt_beta        
    
    # Use actual (or Decon-Cell predicted) cell counts to simulate the expression levels. Parse the file
    # File should be in format of
    #          CC1    CC2
    # sample1  75     14
    # sample2  84     4
    def __simulate_expression(self):
        '''Simulate expression by summing over the simulated expression per cell type + error'''
        outfile = self.out_dir+'/expression/simulated_expression_'+self.batch_name+'.txt'
        with open(outfile,'w') as out:
            out.write('\t'+'\t'.join(self.samples)+'\n')

            for beta_type in self.beta_types:
                logging.info('simulating for beta type '+beta_type)
                for error_type in self.error_size:
                    for index, snp in enumerate(self.genotypes):
                        out.write('gene_'+str(index)+'_'+beta_type+'_'+str(error_type))
                    
                        # When reading the genotype and cell count data we checked that the header is in same order, so an just use the index
                        for sample_index, sample in enumerate(self.samples):
                            dosage = self.genotypes[snp][sample_index]

                            cellcounts = self.cellcount_per_sample[sample]
                            # Expression will be made with expression = cc1 + cc2 + snp:cc1 + snp:cc2 + error
                            # so start with 0
                            expression = 0
                            # then add cc and cc*snp. for cc*snp can add the beta
                            for cc_index, cellcount in enumerate(cellcounts):
                                cc_contribution = self.cc_beta_per_type[beta_type][snp][cc_index] * float(cellcount)
                                expression += cc_contribution
                                cc_snp_contribution = self.cc_gt_beta_per_type[beta_type][snp][cc_index] * float(cellcount) * float(dosage)
                                expression += cc_snp_contribution
                                if cc_contribution+cc_snp_contribution < 0:
                                    logging.info('beta, errror type: '+beta_type+', '+str(error_type))
                                    logging.info('snp: '+snp)
                                    logging.info('cc_index: '+str(cc_index))
                                    logging.info('cellcount: '+str(cellcount))
                                    logging.info('cc_beta: '+str(self.cc_beta_per_type[beta_type][snp][cc_index]))
                                    logging.info('cc_gt_beta: '+str(self.cc_gt_beta_per_type[beta_type][snp][cc_index]))
                                    logging.info('cc_contribution: '+str(cc_contribution))
                                    logging.info('cc_snp_contribution: '+str(cc_snp_contribution))                                               
                                    raise RuntimeError('cc contribution should never be < 0')
                            
                            
                            if expression < 0:
                                raise RuntimeError('Expression should never be lower than 0')
                            if expression+self.error[error_type][snp] < 0:
                                # since expression should never be 0, if with the error it is lower than zero flip the sign of the error
                                # Will give a bias in errors towards positive errors so possibly has to be done differently
                                self.error[error_type][snp] = abs(self.error[error_type][snp])
                            out.write('\t'+str(expression+self.error[error_type][snp]))
                        out.write('\n')
        logging.info('Simulated expression written to '+outfile)
               
    


if __name__ == '__main__':
    expression_simulator = Expression_simulator(args.genotype_file, args.out_dir, args.number_of_snps, args.batch)
