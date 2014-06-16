import math

def test(n, beta):
  minl = 2
  if beta > 0:
    minl = minl + math.ceil(math.log(beta, 2));
  for l0 in range(minl, 1000):
    for l1 in range(minl, 1000):
      q0 = 2**l0
      q1 = 2**l1
      A = (q0/2 - 1)*(n**2)*l0*q0*l1*beta*(n*l1+1)
      B = 1/2 + n*l1*(beta+1/2)
      res = A/q1 + B/q0
      if res < 1/4:
        recrypt_mult = n*l0*(2**(2*l0)) / 1000000
        # print((n**3.5)*l0*(2**(2*l0))*(l1**2) / 1000000000) # total seconds (assuming 10^9 operations per second)
        # print((n**3.5)*l0*(2**(2*l0))*(l1**2)) # total operations
        recrypt_time = math.ceil((n**3.5)*l0*(2**(2*l0))*(l1**2) / 1000000000 / 60)
        if recrypt_time < 60*3: # less than 3 hours -> show in minutes
          recrypt = str(recrypt_time)+" min"
        elif recrypt_time < 60*72: # less than 72 hours -> show in hours
          recrypt = str(math.ceil(recrypt_time/60))+" h"
        else: # else show in days
          recrypt = str(math.ceil(recrypt_time/60/24))+" d"        
        # print((n**3)*l0*(l1**2)*(2**l0)/8/1024) (used KB's)
        memory_usage = math.ceil((n**3)*l0*(l1**2)*(2**l0)/8/1024/1024)
        if memory_usage < 1024: # less than 1 GB -> show in MB
          memory = str(memory_usage) + " MB"
        else: # else show in GB
          memory = str(math.ceil(memory_usage/1024)) + " GB"
        print(n,beta,l0,l1,recrypt_mult,recrypt,memory)
        return

for n in range(2, 12): # [10, 50, 100]:
  beta = math.ceil(math.sqrt(n))
  # beta = 0
  test(n, beta)