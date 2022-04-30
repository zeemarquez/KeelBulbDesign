def y_naca(x,t):
  return 5*t*(0.2969*(x**0.5) - 0.126*x - 0.3516*(x**2) + 0.2843*(x**3) - 0.1015*(x**4) )

def arange(start,end,N):
  step = (end-start)/(N-1)
  return [start + i*step for i in range(N) ]

def genNACAxyz(xx,c,z=0,N=200):
  t = xx/100

  xtab = []
  y_pos = []
  y_neg = []

  for x in arange(0,1.0,N):

    xtab.append(x*c)
    y_ = y_naca(x,t)*c
    y_pos.append(y_)
    y_neg.append(-y_)

  xtab = xtab + xtab[::-1]
  ytab = y_pos + y_neg[::-1]
  ztab = [z for i in range(len(xtab))]
  return (xtab,ytab,ztab)

def find_closest(y_find,y_array):
  y0 = 0
  i = 0
  for y in y_array:
    if y_find == y:
      return (y,i)
    elif (y_find >= y0)and(y_find < y):
      return (y,i)
    i += 1