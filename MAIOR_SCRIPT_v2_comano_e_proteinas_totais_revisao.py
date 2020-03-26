
from collections import Counter
import urllib.request
from bs4 import BeautifulSoup as bs
import re
import timeit
import datetime

antibiotic_classes=['ocorrencia_aminoglycoside_antibiotic', 'ocorrencia_various_substrates', 'ocorrencia_poorly_documented', 'ocorrencia_multidrug', 'ocorrencia_tetracycline_antibiotic', 'ocorrencia_beta-lactam_antibiotic', 'ocorrencia_phenicol_antibiotic', 'ocorrencia_fluoroquinolone_antibiotic', 'ocorrencia_macrolide_antibiotic', 'ocorrencia_glycopeptide_antibiotic', 'ocorrencia_peptide_antibiotic', 'ocorrencia_nitroimidazole_antibiotic']
mechanisms=['ocorrencia_antibiotic_target_alteration', 'ocorrencia_antibiotic_inactivation', 'ocorrencia_ABC_Transporter', 'ocorrencia_poorly_documented', 'ocorrencia_antibiotic_efflux', 'ocorrencia_antibiotic_target_replacement', 'ocorrencia_MFS_efflux', 'ocorrencia_antibiotic_target_protection', 'ocorrencia_RND_Efflux', 'ocorrencia_target_overexpression']
protein_types=['ocorrencia_RND_mfp', 'ocorrencia_Modifying_Enzyme', 'ocorrencia_Transporter_protein', 'ocorrencia_poorly_documented', 'ocorrencia_Two-component_regulatory_system', 'ocorrencia_Hydrolase', 'ocorrencia_Gene_Modulating_Resistance', 'ocorrencia_Mutation_of_target', 'ocorrencia_Target_protection', 'ocorrencia_Hydroxylase']


################### Extracting infos from Patric table ########################

###############################################################################
infos={}
a=0
year={}
patric_cds={}

patric_file="/home/willian/Documentos/Projetos/Resistencia_sepse/para_artigo/revisão_artigo/script_revisto/testes/tabela_sample.tsv"
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

resfams_table='../../Resfams_com_classes_antibioticos_e_mecanismo - Resfams_com_classes_antibioticos_e_mecanismo.tsv'
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
 
out_file='out_table.tsv'
out=open(out_file, 'w')
out.write('id\tname\tGenome_status\tAssembly_accession\tyear\tTotal_protein_number\ttaxon\tPhylum\tclass\torder\tfamily\tgenus\toccurrences_aminoglycoside_antibiotic\toccurrences_various_substrates\toccurrences_poorly_documented\toccurrences_multidrug\toccurrences_tetracycline_antibiotic\toccurrences_beta-lactam_antibiotic\toccurrences_phenicol_antibiotic\toccurrences_fluoroquinolone_antibiotic\toccurrences_macrolide_antibiotic\toccurrences_glycopeptide_antibiotic\toccurrences_peptide_antibiotic\toccurrences_nitroimidazole_antibiotic\toccurrences_antibiotic_target_alteration\toccurrences_antibiotic_inactivation\toccurrences_ABC_Transporter\toccurrences_poorly_documented\toccurrences_antibiotic_efflux\toccurrences_antibiotic_target_replacement\toccurrences_MFS_efflux\toccurrences_antibiotic_target_protection\toccurrences_RND_efflux\toccurrences_target_overexpression\toccurrences_RND_mfp\toccurrences_Modifying_Enzyme\toccurrences_Transporter_protein\toccurrences_poorly_documented\toccurrences_Two-component_regulatory_system\toccurrences_Hydroxylase\toccurrences_Gene_Modulating_Resistance\toccurrences_Mutation_of_target\toccurrences_Target_protection\thidroxylase\n')
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
start = timeit.default_timer()
for elemento in infos:
	w+=1
	progress=(w/len(infos))*100
	stop = timeit.default_timer()
	estimative=((100-progress)*(float(stop)-float(start)))/(progress*60)
	print('%.2f' %progress,"%",'\tProcessing time: ', '%.1f' %((float(stop)-float(start))/60),'m\tEstimated time remaining: ', datetime.timedelta(minutes=estimative) ,end='\r' )
	temp=infos[elemento][1].split(";")
	tx=temp[0]
	cont2[tx]=cont2.get(tx, 0)+1


	fp = urllib.request.urlopen("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="+tx)
	soup = bs(fp, "html.parser")
	try:
		temp=soup.find_all('form')
	except:
		print(tx)
	try:
		a=str(temp[0])
	except:
		print(temp)

	m = re.findall(r'title="(.+?)<', a)
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
			tax_phylum[elemento]=n[1]
		if n[0] == 'class':
			tax_classe[elemento]=n[1]
		if n[0] == 'order':
			tax_order[elemento]=n[1]
		if n[0] == 'family':
			tax_family[elemento]=n[1]
		if n[0] == 'genus':
			tax_genus[elemento]=n[1]
	if elemento not in tax_classe:
		tax_classe[elemento]='unclassified '+tax_phylum[elemento]
	if elemento not in tax_order:
		tax_order[elemento]='unclassified '+tax_classe[elemento]
	if elemento not in tax_family:
		tax_family[elemento]='unclassified '+tax_order[elemento]
	if elemento not in tax_genus:
         tax_genus[elemento]='unclassified '+tax_family[elemento]

############################ Reading HMM files ###############################

##############################################################################
        
	name=infos[elemento][0]
	taxon_organism=infos[elemento][1]
	Genome_status=infos[elemento][2]
	accession=infos[elemento][3]
	genes_file=[]
	antibiotic_class_file=[]
	protein_type_file=[]
	mechanism_file=[]
	evalue={}

	resfams={}
	lineA=elemento+'\t'+name+'\t'+Genome_status+'\t'+accession+'\t'+year[elemento]+'\t'+patric_cds[elemento]+'\t'+taxon_organism+'\t'+tax_phylum[elemento]+'\t'+tax_classe[elemento]+'\t'+tax_order[elemento]+'\t'+tax_family[elemento]+'\t'+tax_genus[elemento]

#Reading every HMM result file
#Place here your path to result files
	for line in open('/home/willian/Documentos/Projetos/Resistencia_sepse/para_artigo/novo_dataset/resultados_hmm/resultados_hmm/'+elemento.replace('"', '')+'._resultado', 'r'):

		
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

	c1=Counter(mechanism_file)
	for elemento in mechanisms:
		elemento=elemento.replace('ocorrencia_', '').replace('_', ' ')
		lineA+='\t'+str(c1[elemento])
        
        
	c1=Counter(protein_type_file)

	for elemento in protein_types:
		elemento=elemento.replace('ocorrencia_', '').replace('_', ' ')
		lineA+='\t'+str(c1[elemento])

	lineA+='\t'+str(len(genes_file))+'\t'+str(genes_file)+'\n'
	out.write(lineA)

lista_class=Counter(classification['class']).most_common()
lista_phylum=Counter(classification['phylum']).most_common()
lista_genus=Counter(classification['genus']).most_common()
lista_order=Counter(classification['order']).most_common()
lista_family=Counter(classification['family']).most_common()


print("\n#phylum:")
for i,j in lista_phylum:
	print(i,"\t",j)

print("\n#class:")
for i,j in lista_class:
	print(i,"\t",j)

print("\n#order:")
for i,j in lista_order:
	print(i,"\t",j)

print("\n#family:")
for i,j in lista_family:
	print(i,"\t",j)

print("\n#genus:")
for i,j in lista_genus:
	print(i,"\t",j)

