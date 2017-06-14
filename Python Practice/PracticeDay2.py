import numpy as np 

# Problem 1
# prob1a = np.ones((4,4), dtype = np.int16)
# prob1a[2][3] = 2
# prob1a[3][1] = 6
# print(prob1a)

# prob1b = np.zeros((5,5))
# n = len(prob1b)
# for i in range(0, n):
# 	for j in range(0, n):
# 		if i == j:
# 			prob1b[i][j] = i + 1
# print(prob1b)

# Problem 2
# prob2 = np.zeros((16,64), dtype = np.int16)
# colsize = len(prob2)
# rowsize = len(prob2[1])

# def Fib(n):
# 	total = 0
# 	if n==1:
# 		return total + 0
# 	elif n==2:
# 		return total + 1
# 	else:
# 		return total + Fib(n-1)+Fib(n-2)
# fibseq = []
# for i in range(1,17):
# 	fibseq.append(Fib(i))
# print(fibseq)

# for i in range(0, colsize):
# 	for j in range(0, rowsize):
# 		prob2[i][j] = fibseq[i]
# print(prob2)

#Problem 3
# unirandnums = np.random.random_integers(2, 16, 20)
# boo = np.zeros((len(unirandnums)), dtype = bool)
# for i in range(0, len(unirandnums)):
# 	if unirandnums[i] >= 5 & unirandnums[i] <= 10 & (unirandnums[i]%2 == 1):
# 		# print(boo[i])
# 		boo[i] = True
# print(boo)

#Problem 4
