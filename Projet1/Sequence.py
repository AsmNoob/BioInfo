class Sequence:
	def __init__(self,sequence):
		self.sequence = sequence
 
	def getAcideAmine(self,indice):
		try:
			return self.sequence[indice]
		except(IndexError):
			print("IndexError")
 
	def longueurSequence(self):
		return len(self.sequence)
 
	def affichage(self):
		print(self.data)