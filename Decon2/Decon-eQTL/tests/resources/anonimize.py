import random
expression = []
with open('expTable_Corrected_addMean_head.txt') as input_file:
    input_file.readline()
    for line in input_file:
        for element in line.strip().split('\t')[1:]:
            expression.append(element)

with open('dsgTable_testing_head.txt') as dsg_table, \
      open('genotype_dosages.txt','w') as out, \
      open('expression_levels.txt','w') as out_expression:
    header = dsg_table.readline().strip().split('\t')
    for Sample in range(1, len(header),1):
        out.write('\t')
        # give each individual a random Sample number
        Sample = 'Sample_'+str(Sample)
        out.write(Sample)
        out_expression.write('\t')
        out_expression.write(Sample)
    out.write('\n')
    out_expression.write('\n')
    number_of_lines = 0
    for line in dsg_table:
        number_of_lines += 1
        if number_of_lines == 50:
            break
        out.write('genotype_'+str(number_of_lines))
        out_expression.write('genotype_'+str(number_of_lines))
        for i in range(1,3,1):
            for element in line.strip().split('\t')[1:]:
                random_number = random.sample(range(1,3,1),1)[0]
                if random_number == 1:
                    element = float(element)-1
                    if element < 0:
                        element = 0
                elif random_number == 2:
                    element = float(element)+1
                    if element > 2:
                        element = float(element)-1
                elif random_number != 3:
                    raise RuntimeError('should only range 1-3, was'+str(random_number))
                if i == 2:
                    out.write('\t'+str(element))
                    out_expression.write('\t'+random.choice(expression))
        out.write('\n')
        out_expression.write('\n')


with open('counts.txt') as input_file, open('cellcounts.txt','w') as out:
    x = 0
    out.write(input_file.readline())
    for line in input_file:
        x += 1 
        line = '\t'.join(line.split('\t')[1:])
        out.write('Sample_'+str(x)+'\t'+line)
