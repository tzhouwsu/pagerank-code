# here I want to create the Adjacency matrix, given a set of occupy/unoccupy site list
# the PBC box of sites looks like:
#   1   2   3   4   1
#   5   6   7   8   5
#   9  10  11  12   9
#  13  14  15  16  13
#   1   2   3   4   1
import numpy as np

# here is to set site-occupy-list, 1 is occupied, 0 is unoccupied
sitelist = [0] * 16
sitelist[0] = 1
sitelist[1] = 1
sitelist[2] = 1
sitelist[3] = 1
sitelist[4] = 1
sitelist[5] = 1
sitelist[6] = 1
sitelist[7] = 1
sitelist[8] = 1
sitelist[9] = 1
sitelist[10] = 0
sitelist[11] = 0
sitelist[12] = 1
sitelist[13] = 1
sitelist[14] = 0
sitelist[15] = 0

# the weight for two-sites as adjacent, if it is second nearest, the weight is half
weight_UU = 1.000
weight_UO = 0.660
weight_OO = 0.330
damp = 0.0

# here is to define neighbor-list (it is fixed given the PBC), each site has 6 nearest neighbors
neighblist = []
neighblist.append([2,4,5,8,13,14]) # neighbor list for site 1
neighblist.append([1,3,5,6,14,15]) # site 2
neighblist.append([2,4,6,7,15,16]) # site 3
neighblist.append([1,3,7,8,13,16]) # site 4
neighblist.append([1,2,6,8,9,12]) # site 5
neighblist.append([2,3,5,7,9,10]) # site 6
neighblist.append([3,4,6,8,10,11]) # site 7
neighblist.append([1,4,5,7,11,12]) # site 8
neighblist.append([5,6,10,12,13,16]) # site 9
neighblist.append([6,7,9,11,13,14]) # site 10
neighblist.append([7,8,10,12,14,15]) # site 11
neighblist.append([5,8,9,11,15,16]) # site 12
neighblist.append([1,4,9,10,14,16]) # site 13
neighblist.append([1,2,10,11,13,15]) # site 14
neighblist.append([2,3,11,12,14,16]) # site 15
neighblist.append([3,4,9,12,13,15]) # site 16


# calculate the adjacency matrix
Adj = []
for i in range(0,16) :
   ai = [0.0] * 16
   sj = [0] * 16
   for k in neighblist[i] :
      sj[k-1] = 1   # j is the nearest neighbor of i
      for j in neighblist[k-1] :
         if sj[j-1] == 0 :
            sj[j-1] = 2  # j is the second nearest neighbor of i
   for j in range(0,16) :
      if sj[j] == 1 :
         if sitelist[i] == 1 and sitelist[j] == 1 :
            ai[j] = weight_OO
         elif (sitelist[i] == 1 and sitelist[j] == 0) or (sitelist[i] == 0 and sitelist[j] == 1):
            ai[j] = weight_UO
         elif sitelist[i] == 0 and sitelist[j] == 0 :
            ai[j] = weight_UU
      elif sj[j] == 2 :
         if sitelist[i] == 1 and sitelist[j] == 1 :
            ai[j] = weight_OO * damp
         elif (sitelist[i] == 1 and sitelist[j] == 0) or (sitelist[i] == 0 and sitelist[j] == 1):
            ai[j] = weight_UO * damp
         elif sitelist[i] == 0 and sitelist[j] == 0 :
            ai[j] = weight_UU * damp
      else :
         ai[j] = 0.0
   Adj.append(ai)
#   print sj
#   print ai

# print(np.matrix(Adj))

# print the adjacency matrix
fout = open("adj-matrix","w")
for i in range(0,16) :
   fout.write("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n" % (Adj[i][0],Adj[i][1],Adj[i][2],Adj[i][3],Adj[i][4],Adj[i][5],Adj[i][6],Adj[i][7],Adj[i][8],Adj[i][9],Adj[i][10],Adj[i][11],Adj[i][12],Adj[i][13],Adj[i][14],Adj[i][15])) 

fout.close()


