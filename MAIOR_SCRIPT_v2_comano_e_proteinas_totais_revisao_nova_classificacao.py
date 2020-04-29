
from collections import Counter
import urllib.request
from bs4 import BeautifulSoup as bs
import re
import timeit
import datetime
import sys
sys.setrecursionlimit(1000000)

antibiotic_classes=['ocorrencia_aminoglycoside_antibiotic', 'ocorrencia_various_substrates', 'ocorrencia_poorly_documented', 'ocorrencia_multidrug', 'ocorrencia_tetracycline_antibiotic', 'ocorrencia_beta-lactam_antibiotic', 'ocorrencia_phenicol_antibiotic', 'ocorrencia_fluoroquinolone_antibiotic', 'ocorrencia_macrolide_antibiotic', 'ocorrencia_glycopeptide_antibiotic', 'ocorrencia_peptide_antibiotic', 'ocorrencia_nitroimidazole_antibiotic']
mechanisms=['ocorrencia_antibiotic_target_alteration', 'ocorrencia_antibiotic_inactivation', 'ocorrencia_ABC_Transporter', 'ocorrencia_poorly_documented', 'ocorrencia_antibiotic_efflux', 'ocorrencia_antibiotic_target_replacement', 'ocorrencia_MFS_efflux', 'ocorrencia_antibiotic_target_protection', 'ocorrencia_RND_Efflux', 'ocorrencia_target_overexpression', 'ocorrencia_SMR_efflux']
protein_types=['ocorrencia_RND_mfp', 'ocorrencia_Modifying_Enzyme', 'ocorrencia_Transporter_protein', 'ocorrencia_poorly_documented', 'ocorrencia_Two-component_regulatory_system', 'ocorrencia_Hydrolase', 'ocorrencia_Gene_Modulating_Resistance', 'ocorrencia_Mutation_of_target', 'ocorrencia_Target_protection', 'ocorrencia_Hydroxylase']




patric_file="/home/willian/Documentos/Projetos/Resistencia_sepse/para_artigo/revisão_artigo/nova_classificacao/Blood_dataset.tsv"
resfams_table='/home/willian/Documentos/Projetos/Resistencia_sepse/para_artigo/revisão_artigo/nova_classificacao/Resfams_classificcoes_nova.csv'
out_file='/home/willian/Documentos/Projetos/Resistencia_sepse/para_artigo/revisão_artigo/nova_classificacao/results_tables/resultados_leitura_hmm.tsv_teste_script'
results_hmm_path='/home/willian/Documentos/Projetos/Resistencia_sepse/para_artigo/novo_dataset/resultados_hmm/resultados_hmm/'



################### Extracting infos from Patric table ########################

###############################################################################
infos={}
a=0
year={}
patric_cds={}

for line in open(patric_file, 'r'):
	temp=line.split('\t')
	year[temp[0]]=temp[36]
	patric_cds[temp[0]]=temp[31]

for line in open(patric_file, 'r'):
#for line in open(sys.argv[1]):
	if a ==0:
		a+=1
		continue
	line=line.rstrip()
	temp=line.split('\t')
	lista=[temp[1], temp[3], temp[4], temp[17]]
	infos[temp[0]]=lista

############ Storing infos of antibiotic, proteins and mechanisms ##################

####################################################################################
protein_type={}
antibiotic_class={}
mechanism={}
a=0

for line in open(resfams_table, 'r'):
	if a ==0:
		a+=1
		continue
	else:
		line=line.rstrip()
		temp=line.split('\t')
		antibiotic_class[temp[0]]=temp[2]
		protein_type[temp[0]]=temp[4]
		mechanism[temp[0]]=temp[3]


#################### Storing infos from taxons  ####################################

####################################################################################

#creating an output file
 
out=open(out_file, 'w')
out.write('id\tname\tGenome_status\tAssembly_accession\tyear\tTotal_protein_number\ttaxon\tPhylum\tclass\torder\tfamily\tgenus\toccurrences_aminoglycoside_antibiotic\toccurrences_various_substrates\toccurrences_poorly_documented\toccurrences_multidrug\toccurrences_tetracycline_antibiotic\toccurrences_beta-lactam_antibiotic\toccurrences_phenicol_antibiotic\toccurrences_fluoroquinolone_antibiotic\toccurrences_macrolide_antibiotic\toccurrences_glycopeptide_antibiotic\toccurrences_peptide_antibiotic\toccurrences_nitroimidazole_antibiotic\toccurrences_antibiotic_target_alteration\toccurrences_antibiotic_inactivation\toccurrences_ABC_Transporter\toccurrences_poorly_documented\toccurrences_antibiotic_efflux\toccurrences_antibiotic_target_replacement\toccurrences_MFS_efflux\toccurrences_antibiotic_target_protection\toccurrences_RND_efflux\toccurrences_target_overexpression\toccurrences_SMR_efflux\toccurrences_RND_mfp\toccurrences_Modifying_Enzyme\toccurrences_Transporter_protein\toccurrences_poorly_documented\toccurrences_Two-component_regulatory_system\toccurrences_Hydrolase\toccurrences_Gene_Modulating_Resistance\toccurrences_Mutation_of_target\toccurrences_Target_protection\thidroxylase\n')
tax_phylum={}
tax_classe={}
tax_order={}
tax_family={}
tax_genus={}
gene=[]

w=0
cont2={}
dictionary={}
classification={}

taxons=[]
for line in open(patric_file, 'r'):
    temp_patric=line.split('\t')
    taxons.append(temp_patric[3])
    
taxons=set(taxons)

print("Storing Infos of taxons")

len_tax=len(taxons)

r=0
tax_dict={'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}}
for tx in taxons:
	fp = urllib.request.urlopen("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="+tx)
	r+=1
	if r%100==0:
		print("Storing tax infos "+str(r)+" of "+str(len_tax))
	soup = bs(fp, "html.parser")
	try:
		temp=soup.find_all('form')
	except:
		print(tx)
	try:
		axt=str(temp[0])
	except:
		print(temp)

	m = re.findall(r'title="(.+?)<', axt)
	classification_list=[]
	check_taxonomy={}
	for n in m:
		n=n.replace('"', "").split(">")
		classification_list.append(n)
		if n[0] not in classification:
			classification[n[0]]=[]
		check_taxonomy[n[0]]=n[1]
		classification[n[0]].append(n[1])
		dictionary[tx]=classification_list
		if n[0] == 'phylum':
			tax_dict['phylum'][tx]=n[1]
		if n[0] == 'class':
			tax_dict['class'][tx]=n[1]
		if n[0] == 'order':
			tax_dict['order'][tx]=n[1]
		if n[0] == 'family':
			tax_dict['family'][tx]=n[1]
		if n[0] == 'genus':
			tax_dict['genus'][tx]=n[1]
	if tx not in tax_dict['class']:
		tax_dict['class'][tx]='unclassified '+tax_dict['phylum'][tx]
	if tx not in tax_dict['order']:
		tax_dict['order'][tx]='unclassified '+tax_dict['class'][tx]
	if tx not in tax_dict['family']:
		tax_dict['family'][tx]='unclassified '+tax_dict['class'][tx]
	if tx not in tax_dict['genus']:
		tax_dict['genus'][tx]='unclassified '+tax_dict['family'][tx]



print("reading HMM files")

for elemento in infos:
	

############################ Reading HMM files ###############################

##############################################################################
        
	name=infos[elemento][0]
	tx=infos[elemento][1]
	Genome_status=infos[elemento][2]
	accession=infos[elemento][3]
	genes_file=[]
	antibiotic_class_file=[]
	protein_type_file=[]
	mechanism_file=[]
	evalue={}

	resfams={}
	lineA=elemento+'\t'+name+'\t'+Genome_status+'\t'+accession+'\t'+year[elemento]+'\t'+patric_cds[elemento]+'\t'+tx+'\t'+tax_dict['phylum'][tx]+'\t'+tax_dict['class'][tx]+'\t'+tax_dict['order'][tx]+'\t'+tax_dict['family'][tx]+'\t'+tax_dict['genus'][tx]

#Reading every HMM result file
#Place here your path to result files
	for line in open(results_hmm_path+elemento.replace('"', '')+'._resultado', 'r'):

		
		if line.startswith("#"):
 			continue
		temp=line.split(" ")
		while "" in temp:
 			temp.remove("")
		if float(temp[6]) > (float(temp[5])/10):
 			continue
		if temp[0] not in evalue:
 			evalue[temp[0]]=temp[4]
 			resfams[temp[0]]=temp[3]
		else:
 			if float(temp[4]) < float(evalue[temp[0]]):
                 		evalue[temp[0]]=temp[4]
                 		resfams[temp[0]]=temp[3]
	for coisa in resfams:
		genes_file.append(resfams[coisa])
		gene.append(resfams[coisa])
	for elemento in genes_file:
		if ',' in antibiotic_class[elemento]:
 			antibiotic_class_file.append('multidrug')
		else:
 			antibiotic_class_file.append(antibiotic_class[elemento])
		protein_type_file.append(protein_type[elemento])
		mechanism_file.append(mechanism[elemento])
	c1=Counter(antibiotic_class_file)

	for elemento in antibiotic_classes:
		elemento=elemento.replace('ocorrencia_', '').replace('_', ' ')
		lineA+='\t'+str(c1[elemento])

	c2=Counter(mechanism_file)
	for elemento in mechanisms:
		elemento=elemento.replace('ocorrencia_', '').replace('_', ' ')
		lineA+='\t'+str(c2[elemento])
        
        
	c3=Counter(protein_type_file)

	for elemento in protein_types:
		elemento=elemento.replace('ocorrencia_', '').replace('_', ' ')
		lineA+='\t'+str(c3[elemento])

	lineA+='\t'+str(len(genes_file))+'\t'+str(genes_file)+'\n'
	out.write(lineA)


