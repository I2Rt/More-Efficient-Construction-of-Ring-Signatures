import math

# delta = 1.0044
pi = 3.1415926535
e = 2.718281828

q = math.pow(2, 32) - 99

eta = 1 # s
xi = 2  # c

n = 4
k = 7
d = 64

kd = k * d
nd = n * d

beta = math.sqrt((eta*eta * k + xi * xi)*d)


delta = math.pow( math.pow(pi * beta, 1/beta) * beta / (2 * pi * e) , 1/(2 * (beta-1) ) )


print("delta2", delta)

print("beta", beta)