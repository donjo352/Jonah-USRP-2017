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
	# if unirandnums[i] >= 5 & unirandnums[i] <= 10 & (unirandnums[i]%2 == 1):
# 		# print(boo[i])
# 		boo[i] = True
# print(boo)

#Problem 4
import matplotlib.pyplot as plt

# def gaussian2D(x, y, mean_xy, cov):
# 	return np.exp((-1.0/2)*np.transpose(x - mean_xy)*(cov**-1.0)*(x - mean_xy))

x = np.arange(0.,10.,0.1)
y = np.arange(0.,10.,0.1)
x,y = np.meshgrid(x,y)
mean_xy = np.array([5.,4.])
cov = np.array([[1.,1.],
             [1.,2.]])**2.
print(x.shape)
print(y.shape)
print(mean_xy.shape)
print(cov.shape)

# f = gaussian2D(x.ravel(), y.ravel(), mean_xy, cov)
# f.reshape(100,100)

# plt.imshow(f.reshape(100,100))
# plt.show()