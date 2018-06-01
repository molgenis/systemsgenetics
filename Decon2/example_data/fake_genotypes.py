import random

# Very simple fake genotype file to showcase how to run Decon method
# Also write the snp-gene file for which combination needs to be tested
with open('count.table.simulated.txt') as input_file, open('gene_snp_file.txt','w') as out, open('fake_genotype_dosages.txt','w') as out2:
    x = 0
    out.write('gene\tsnp\n')
    out2.write(input_file.readline())
    for line in input_file:
        # Don't want to test all combinations so only add 10% of them. 
        if random.random() > 0.9:
            expression_values = line.split('\t')[1:]
            # only add if not all expression levels are 0
            # because otherwise Decon-eQTL will give an error due to no variation in the data
            if sum([float(x) for x in expression_values]) == 0:
                continue
            x += 1
            out.write(line.split('\t')[0]+'\tsnp_'+str(x)+'\n')
            out2.write('snp_'+str(x))
            for i in range(0, len(line.split('\t'))-1):
                out2.write('\t'+str(random.randint(0,2)))
            out2.write('\n')
