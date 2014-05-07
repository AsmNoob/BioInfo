from math import sqrt

def parser(filename):
	listeSequences = []
	try:
		fichier = open(filename,"r")
		sequence = ""
		for ligne in fichier:
			if(ligne[0] != ">"):
				sequence+=ligne.strip("\n")
			else:
				if(sequence != ""):
					listeSequences.append(sequence)
					sequence = ""
		return listeSequences
	except IOError as e:
		print(e)
		return None

def calculAA(liste):
	listeAcidesAmines = []
	for sequence in liste:
		for aa in sequence:
			if(aa not in listeAcidesAmines and aa != "-"):
				listeAcidesAmines.append(aa)
	return listeAcidesAmines

def MatriceAAPosition(liste):
	listeAA = calculAA(liste)
	mat = [[0 for i in range(len(liste[0]))] for j in range(len(listeAA))]
	for sequence in (liste):
		for indicePosition,aa in enumerate(sequence):
			mat[listeAA.index(aa),indicePosition]+=1
	return mat

def calculDesFrequences(matrice,nombreSequences):
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			matice[i,j] /= nombreSequences
	return matrice

def ajoutPseudoCounts(matrice,nombreSequences,listeSwissProt,listeAA):
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			matrice[i,j] = (nombreSequences - 1)*(matrice[i,j]) + (sqrt(nombreSequences))*(listeSwissProt[listeAA[i]])/(sqrt(nombreSequences)+(nombreSequences - 1))
	return matrice		

def matrice_M(matrice,listeSwissProt,listeAA):
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			matrice[i,j] = log(matrice[i,j]/listeSwissProt[listeAA[i]])


listeSwissProt = {"A":8.28, "Q":3.94,"L":9.67,"S":6.50,"R":5.53,"E":6.76,"K":5.85,"T":5.32,"N":4.05,"G":7.09,"M":2.43,"W":1.07,"D":5.45,"H":2.27,"F":3.86,"Y":291,"C":1.36,"I":5.99,"4.68","V":6.87,"B":0,"Z":0,"X":0}
	
MatriceAAPosition(parser("msaresultd-MUSCLE.fasta"))
