def parser(filename):
	listeSequences = []
	try:
		fichier = open(filename,"r")
		sequence = ""
		for ligne in fichier:
			if(ligne[0] != ">"):
				sequence+=ligne
			else:
				if(sequence != ""):
					listeSequences.append(sequence)
					sequence = ""
		return listeSequences
	except IOError as e:
		print(e)
		return None

liste = (parser("msaresultd-MUSCLE.fasta"))
for i in range(1):
	print(liste[i])

