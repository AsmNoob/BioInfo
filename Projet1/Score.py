class Score(list):
    #On initialise le self, comme la matrice contenant la matrice lue.
    def __init__(self, nomFichier):
        with open(nomFichier) as _file:
            _file = _file.read().strip().split('\n') 
            self.extend([[(element) for element in line.split()] for line in _file])
    
    #Renvoie la valeur à une coordonnée précise
    def getValeurIndice(self,id1,id2):
    	return self[id1][id2]

    #Renvoie les coordonnées entre 2 protéines
    def getIndiceProteines(self,protein1,protein2):
    	return (self[0].index(protein1),self[0].index(protein2))

    #Renvoie la valeur se trouvant à l'intersection de ces deux protéines
    def getValeurProteines(self,protein1,protein2):
    	return self[self[0].index(protein1)][self[0].index(protein2)]

    #Place la valeur donnée en paramètre aux indices donnés
    def setValeurIndice(self,id1,id2,valeur):
    	self[id1][id2] = valeur

    #Place la valeur donnée en paramètre à l'intersection de ces deux protéines
    def setValeurProteines(self,proteine1,proteine2,valeur):
    	self[self[0].index(protein1)][self[0].index(protein2)] = valeur


    #Affichage de la matrice
    def affichageMatrice(self):
    	for i in range(len(self)):
    		print()
    		for j in range(len(self[0])):
    			print(self[i][j], end = ' ')




