#Binary Converter

'''
Conversion of binary to real number
'''
#------------------------------------------------------------------------------
def fromFixedPoint(w: int, b:int, bits:[int]):
   # w: width of the binary representation
   # b: binary point
   integer = [bits[i] for i in range(w-b)]
   integer_value = [-integer[i]*2**(w-b-1-i) if i==0 else 
                   integer[i]*2**(w-b-1-i) for i in range(w-b)]
   fractional = [bits[i] for i in range(w-b, w)]
   fractional_value = [fractional[i]*2**(-1-i) for i in range(b)]
   return sum(integer_value) + sum(fractional_value)
#------------------------------------------------------------------------------
print(fromFixedPoint(10, 3, [0, 1, 0, 1, 1, 0, 0, 1, 1, 0]))
print(fromFixedPoint(10, 5, [1, 0, 0, 1, 0, 1, 0, 1, 1, 1]))
print(fromFixedPoint(8, 2, [1, 0, 1, 0, 1, 0, 1, 1]))
