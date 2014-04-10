def BLOSUM(listeFichiers,pourcentage):
	listeMatrices = []
	for nomFichier in listeFichiers:
		#######listeSequences#######
		listeSequences = recupererSequences(nomFichier)
		#######listeProteinesDifferentes#######
		#listeProteine = listeProteines(listeSequences)
		#print(listeProteine)
		#######groupesDifferents#######
		listeGroupes = trouverGroupesDifferents(listeSequences,pourcentage)
		#######ajoutListeMatrices#######
		listeMatrices.append(matriceFrequencesPonderees(listeGroupes,['C', 'P', 'S', 'V', 'Y', 'R', 'L', 'M', 'Q', 'W', 'N', 'D', 'F', 'E', 'K', 'A', 'T', 'H', 'G', 'I']))
	#######moyenneDesMatrices#######
	matriceResultante = moyenneMatrices(listeMatrices)
	#######probabilitésD'Occurences#######
	matriceOccurences = calculProbabiliteOccurence(matriceResultante)
	####### #######

def calculProbabiliteOccurence(matrice):
	sommePartieSuperieure = calculPartieSuperieureMatrice(matrice)
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			matrice[i][j]/= sommePartieSuperieure
	return matrice

def calculPartieSuperieureMatrice(matrice):
	somme = 0
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			if(j >= i):
				print("i")
				somme+=matrice[i][j]
	return somme


#Moyenne des matrices passées en paramètres
def moyenneMatrices(listeMatrices):
	matrice = []
	for i in range(len(listeMatrices[0])):
		matrice.append([0 for m in range(len(listeMatrices[0][0]))])
		for j in range(len(listeMatrices[0][0])):
			for n in range(len(listeMatrices)):
				matrice[i][j] += listeMatrices[n][i][j]
			matrice[i][j] /= len(listeMatrices)
	return matrice

#Récupère les séquences dans un fichier donné
def recupererSequences(nomFichier):
	try:
		fichier = open(nomFichier,"r")
		listeSequences = fichier.readlines()
		for i in range(len(listeSequences)):
			listeSequences[i] = listeSequences[i].strip("\n")
		fichier.close()
		return listeSequences
	except (IOError):
		print ("IOError") 
		return None

#Place les sequences dans des groupes différents selon leur similitudes
def trouverGroupesDifferents(listeSequences, pourcentage):
	listeGroupes = []
	for sequence in listeSequences:
		if(listeGroupes == []):
			liste = [sequence]
			listeGroupes.append(liste)
		else:
			groupeTrouve = False
			for groupe in listeGroupes:
				if(not groupeTrouve):
					for i in range(len(groupe)):
						if(comparaisonSequences(sequence,groupe[i], pourcentage) and not groupeTrouve):
							groupe.append(sequence)
							groupeTrouve = True
			if(not groupeTrouve):
				listeM = [sequence]
				listeGroupes.append(listeM)
				groupeTrouve = True
	return listeGroupes

#Compare 2 sequences et si leur taux d'dentité est supérieur au pourcentage donné.
def comparaisonSequences(seq1,seq2,pourcentage):
	elementDeSimilitude = 0
	longueurSequence = len(seq1)
	res = False
	for i in range(longueurSequence):
		if(seq1[i] == seq2[i]):
			elementDeSimilitude+=1
	if((elementDeSimilitude/longueurSequence) >= pourcentage):
		res = True
	return res

#Calcul de la matrice des frequences pondérées
def matriceFrequencesPonderees(listeGroupes,listeProteines):
	matrice = []
	for i in range(len(listeProteines)):
		matrice.append([0 for m in range(len(listeProteines))])
		for j in range(len(listeProteines)):
			#print("proteine1: "+ listeProteines[i]+" proteine2: "+ listeProteines[j] )
			#print(calculFrequence(listeProteines[i],listeProteines[j],listeGroupes))
			matrice[i][j] = calculFrequence(listeProteines[i],listeProteines[j],listeGroupes)
	#print(matrice)
	return matrice


def colonneProteines(listeGroupes,indice):
	colonne = []
	for groupe in listeGroupes:
		for sequence in groupe:
			colonne.append(sequence[indice])
	return colonne

#Calcule la fréquence pondérée pour 2 proteines données 
def calculFrequence(proteine1,proteine2,listeGroupes):
	somme = 0
	for indice in range(len(listeGroupes[0][0])): #taille d'une sequence
		colonne = colonneProteines(listeGroupes,indice)#colonne de proteines sur tous les groupes
		for i in range(len(colonne)):
			if(colonne[i] == proteine1):
				if(proteine1 == proteine2):
					for j in range(i,len(colonne)):
						if(j != i and colonne[j] == proteine2 and calculGroupe(i,listeGroupes) != calculGroupe(j,listeGroupes)):
							#if(proteine1 == "e"):
								#print(calculGroupe(i,listeGroupes))
								#print(calculGroupe(j,listeGroupes))
								#print(calculPoidsGroupe(i,listeGroupes))
								#print(calculPoidsGroupe(j,listeGroupes))
								#print("-----------------------")
							somme+=(calculPoidsGroupe(i,listeGroupes)*(calculPoidsGroupe(j,listeGroupes)))
				else:
					for j in range(len(colonne)):
						if(j != i and colonne[j] == proteine2 and calculGroupe(i,listeGroupes) != calculGroupe(j,listeGroupes)):
							somme+=(calculPoidsGroupe(i,listeGroupes)*(calculPoidsGroupe(j,listeGroupes)))
	return somme


def calculGroupe(indiceColonne,listeGroupes):
	numeroGroupe = 0
	positionActuelle = len(listeGroupes[numeroGroupe]) #taille du premier groupe
	for groupe in listeGroupes:
		if(indiceColonne >= positionActuelle):
			numeroGroupe+=1
			positionActuelle+= len(listeGroupes[numeroGroupe])
	return numeroGroupe




#calcul du poids de la sequence, en utilisant son indice ds la colonne et la taille des différents groupes
def calculPoidsGroupe(indiceColonne,listeGroupes):
	numeroGroupe = 0
	positionActuelle = len(listeGroupes[numeroGroupe]) #taille du premier groupe
	for groupe in listeGroupes:
		if(indiceColonne >= positionActuelle):
			numeroGroupe+=1
			positionActuelle+= len(listeGroupes[numeroGroupe])
	return (1/(len(listeGroupes[numeroGroupe])))



def listeProteines(listeSequences):
    liste = []
    for sequence in listeSequences:
        for proteine in sequence:
            if (proteine not in liste):
                liste.append(proteine)
    return liste


"""sequences = ["atckq","atcrn","asckn","sscrn","sdceq","secen","tecrq"]
listeGroupes = [["atckq","atcrn","asckn","sscrn"],["sdceq","secen"],["tecrq"]]
print(listeGroupes)
print(listeProteines(sequences))
print(matriceFrequencesPonderees(listeGroupes,["a","c","d","e","k","n","q","r","s","t"]))"""

#BLOSUM(["PR00109A.txt","PR00109B.txt","PR00109C.txt","PR00109D.txt","PR00109E.txt"],0.5)