class SomeClass(object):

	def __init__(self, **kwargs):
		self.kwargs = kwargs

	def __getitem__(self, key):
		return self.kwargs[key]

s = SomeClass(a=5., b=0, c="whatever")
print(s.kwargs)
print(s["a"])

class Temperature(object):

	def __init__(self, value_kelvin):
		self.value_kelvin = value_kelvin

	def __getitem__(self, key):
		if key == 'celsius':
			return self.value_kelvin - 273

		elif key == 'fahrenheit':
			return (self.value_kelvin - 273.)*9/5 +32

		elif key == 'kelvin':
			return self.value_kelvin

		else: 
			raise KeyError('Invalid Key.')

t = Temperature(298.)
print(t["celsius"])