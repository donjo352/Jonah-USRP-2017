# def addme(x, y):
# 	return x + y
# addme(4,3)

# def addme2(x, y=2):
# 	return x + y 
# addme2(3)

# d = dict(x=20, y=30)
# print(addme2(**d))

# def myfunc(*args, **kwargs):
# 	print(args)
# 	print(kwargs)
# # myfunc(1,2,3,4,5)
# myfunc(1,2,3, a=10, b=20)

# def myplot(x, y, **kwargs):

# 	x *= 10
# 	y *= 100
# 	print(kwargs)
# 	plt.plot(x, y, **kwargs)

# import matplotlib.pyplot as plt 
# import numpy as np

# x = np.arange(10)
# y = np.arange(10)
# print(x, y)
# myplot(x, y, color='red', lw=5)
# plt.show()

s = 'Jonah I. Donnenfield'
# print(type(s))
# print(s.split())
# print(s.upper())
from math import sqrt
class Coord(object):

	def __init__(self, x, y): #__init__ is a constructor
		print("I am running the constructor")
		self.x = x
		self.y = y

	def distance_from_origin(self):
		return sqrt(self.x**2 + self.y**2)

	def __add__(self, other):
		assert type(other) is Coord, "Must add another Coord"

		newx = self.x + other.y
		newy = self.y + other.x
		return Coord(self.x + other.x)

c = Coord(4, 5)
c2 = Coord(10, 20)
import numpy as np
# print(np.arange(12).reshape((3,4)))
# print(c + 2)
# print(c.x, c.y)
# print(type(c))
# print(c.distance_from_origin())

import matplotlib.pyplot as plt
x = np.linspace(0, 1, 100)
y = x**2
print(y[5])
plt.plot(x,y)
plt.show()





