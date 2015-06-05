# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 20:30:42 2015

@author: jonesad
"""
#factor an int
def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

#return 2 most balanced factors of n
def balFact(n):
  fac=factors(n)
  lBest=[n,1]
  for elem in fac:
    if abs(elem-n/elem) < lBest[0]-lBest[1]:
      lBest=[max(elem, n/elem),min(elem,n/elem)]
  return lBest