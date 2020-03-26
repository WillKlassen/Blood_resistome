import sys

out = open('filtered_table.txt', 'w')
a=0
for line in open(sys.argv[1], 'r'):
	temp=line.split('\t')
	if a == 0:
		b=0
		for elemento in temp:
			if 'blood' in elemento:
				print(b,'\t',elemento)
			b+=1
		a+=1
	else:
		if len(temp) >= 48:
			temp[33]=temp[33].lower()
			if (('blood' in temp[33]) or ('blood' in temp[34]) or ('blood' in temp[48]) or ('blood' in temp[49])):
				if temp[4] == '':
					continue
				if temp[3] == '':
					continue
				out.write(line)

# 