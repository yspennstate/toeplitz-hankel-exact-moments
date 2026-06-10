# Quantum_Gossip (modified) ---------------------------------------------------
# Lightly edited version of Quantum_Gossip_original.R for the project
# "Exact higher moments of the Toeplitz and Hankel limiting spectral
# distributions" (problem01_toeplitz_hankel_moments).
#
# What this script does (unchanged): builds all planar (non-crossing) perfect
# matchings of 2N boundary points by inserting one edge at a time, checks the
# count equals the Catalan number choose(2N,N)/(N+1), represents the
# Temperley-Lieb generators e_i on the diagram basis with loop weight
# delta = 2 (each closed loop formed when stacking contributes a factor 2),
# and verifies the Temperley-Lieb axioms:
#     e_i^2 = 2 e_i ;  e_i e_j = e_j e_i (|i-j| >= 2) ;  e_i e_{i+-1} e_i = e_i.
#
# Relevance: the non-crossing matchings enumerated here are exactly the
# volume-1 pairings of Theorem 3.1 in the accompanying paper, and the Catalan
# counts 2, 5, 14, 42, 132, 429 are the same numbers appearing in the paper's
# volume spectra. See tl_noncrossing.py for the Python port used as an
# independent cross-check against the exact volume data.
#
# Changes relative to the original ("a little"):
#   - reproducibility: set.seed(20260610) before all sampling
#   - N renamed comment: "dimensions" (was "dimentions"); N kept = 10 default
#   - explicit Catalan check with a printed PASS/FAIL line (was a bare ratio)
#   - axiom checks collected into named booleans and printed as a summary
#   - guarded the plotting blocks behind `do_plots` so the script runs headless
# The combinatorial core (putting_another_edge, Graph_Multiplication, d) is
# untouched.

library("igraph")
set.seed(20260610)
do_plots <- FALSE

# main functions
N = 10        # dimensions (number of strands; 2N boundary points)
number_of_nodes = N

putting_another_edge = function(matrix, connection_length){
  z <- matrix
  root = 1
  while(sum(matrix[root,])>=1){
    root = root +1
  }
  if(root<=N){
    possible_connections = array()
    counter = 1
    if(N%%2==0){
      for (i in 1:(2*N)){
        if(i<=N){
          if(abs(i-root)%%2 != 0){
            possible_connections[[counter]] = i
            counter = counter +1
          }
        }else{
          if(abs(i-root)%%2 == 0){
            possible_connections[[counter]] = i
            counter = counter +1
          }
        }
      }
    }else{
      for (i in 1:(2*N)){
        if(abs(i-root)%%2 != 0){
          possible_connections[[counter]] = i
          counter = counter +1
        }
      }
    }

    for(i in possible_connections){
      if(sum(matrix[,i])==1 | sum(matrix[i,])==1){
        possible_connections = possible_connections[-which(possible_connections==i)]
      }
    }

    Maxtop = N+1
    minbottom = N+1
    for (i in 1:(root-1)){
      Maxtop = max(Maxtop,which(matrix[i,]==1))
      if(which(matrix[i,]==1)<=N & which(matrix[i,]==1)>root){
        minbottom =  min(minbottom,which(matrix[i,]==1))
      }
    }

    if(minbottom>root  & minbottom !=N+1){
      a = which(possible_connections < minbottom)
      b = which(possible_connections > root)
      c = intersect(a,b)
      possible_connections = possible_connections[c]
    }else{
      for (i in possible_connections){
        if(i> N & i<Maxtop | i<root){
          a = which(possible_connections==i)
          possible_connections = possible_connections[-a]
        }
      }
    }

    if(possible_connections[min(length(possible_connections),connection_length)]<=N){
      matrix[root,possible_connections[min(length(possible_connections),connection_length)]] = 1
      matrix[possible_connections[min(length(possible_connections),connection_length)], root] = 1
    }else{
      matrix[root,possible_connections[min(length(possible_connections),connection_length)]] = 1
    }
  }else{
    while(sum(matrix[,root])==1 | sum(matrix[root,])==1){
      root = root +1
    }

    possible_connections = array()
    counter = 1
    for (i in (N+1):(2*N)){
      if(abs(i-root)%%2 != 0){
        possible_connections[[counter]] = i
        counter = counter +1
      }
    }

    start_cutoff= 2*N+1
    nearest_cutoff = start_cutoff

    for (i in 1:(2*N-root)){
      if(sum(matrix[,start_cutoff-i]) >= 1  ){
        if(which(matrix[,start_cutoff-i]==1)<root){
          nearest_cutoff = 2*N - i + 1
        }
      }
    }

    for(i in possible_connections){
      if(sum(matrix[,i])==1){
        possible_connections = possible_connections[-which(possible_connections==i)]
      }
    }

    a = which(possible_connections < nearest_cutoff)
    b = which(possible_connections > root)
    c = intersect(a,b)
    possible_connections = possible_connections[c]

    matrix[root,possible_connections[min(length(possible_connections),connection_length)]] = 1
    matrix[possible_connections[min(length(possible_connections),connection_length)], root] = 1
  }

  for(i in 1:length(matrix[1,])){
    if (sum(matrix[,i])>1){
      print(z)
    }
  }
  return(matrix)
}

Finding_number_of_possible_edges =  function(matrix, connection_length){
  root = 1
  while(sum(matrix[root,])>=1){
    root = root +1
  }
  if(root<=N){
    possible_connections = array()
    counter = 1
    if(N%%2==0){
      for (i in 1:(2*N)){
        if(i<=N){
          if(abs(i-root)%%2 != 0){
            possible_connections[[counter]] = i
            counter = counter +1
          }
        }else{
          if(abs(i-root)%%2 == 0){
            possible_connections[[counter]] = i
            counter = counter +1
          }
        }
      }
    }else{
      for (i in 1:(2*N)){
        if(abs(i-root)%%2 != 0){
          possible_connections[[counter]] = i
          counter = counter +1
        }
      }
    }

    for(i in possible_connections){
      if(sum(matrix[,i])==1 | sum(matrix[i,])==1){
        possible_connections = possible_connections[-which(possible_connections==i)]
      }
    }

    Maxtop = N+1
    minbottom = N+1
    for (i in 1:(root-1)){
      Maxtop = max(Maxtop,which(matrix[i,]==1))
      if(which(matrix[i,]==1)<=N & which(matrix[i,]==1)>root){
        minbottom =  min(minbottom,which(matrix[i,]==1))
      }
    }

    if(minbottom>root  & minbottom !=N+1){
      a = which(possible_connections < minbottom)
      b = which(possible_connections > root)
      c = intersect(a,b)
      possible_connections = possible_connections[c]
    }else{
      for (i in possible_connections){
        if(i> N & i<Maxtop | i<root){
          a = which(possible_connections==i)
          possible_connections = possible_connections[-a]
        }
      }
    }
  }else{
    while(sum(matrix[,root])==1 | sum(matrix[root,])==1){
      root = root +1
    }
    possible_connections = array()
    counter = 1
    for (i in (N+1):(2*N)){
      if(abs(i-root)%%2 != 0){
        possible_connections[[counter]] = i
        counter = counter +1
      }
    }

    start_cutoff= 2*N+1
    nearest_cutoff = start_cutoff
    for (i in 1:(2*N-root)){
      if(sum(matrix[,start_cutoff-i]) >= 1  ){
        if(which(matrix[,start_cutoff-i]==1)<root){
          nearest_cutoff = 2*N - i + 1
        }
      }
    }
    for(i in possible_connections){
      if(sum(matrix[,i])==1){
        possible_connections = possible_connections[-which(possible_connections==i)]
      }
    }
    a = which(possible_connections < nearest_cutoff)
    b = which(possible_connections > root)
    c = intersect(a,b)
    possible_connections = possible_connections[c]
  }
  return(length(possible_connections))
}

f = function(x){
  counter = 1
  matrixies_new = list()
  number_of_edges_left= array()
  for (i in 1:length(x)){
    number_of_edges_left[[i]] = Finding_number_of_possible_edges(x[[i]],1)
  }
  for(j in 1:length(x)){
    for(i in 1:(number_of_edges_left[[j]])){
      matrixies_new[[counter]] = putting_another_edge(x[[j]],i)
      counter = counter +1
    }
  }
  return(matrixies_new)
}

Generators = list()
for(i in 1:N){
  Generators[[i]]= matrix(0L,nrow= 2*N, ncol = 2*N)
  if(i<=floor(N/2)){
    Generators[[i]][1,2*i] = 1
    Generators[[i]][2*i,1] = 1
  }else{
    Generators[[i]][1,(N+2*(i-floor(N/2))-1)] = 1
  }
}

a = Generators
for (i in 1:(N-1)){
  a = f(a)
  a = unique(a)
}

# --- explicit Catalan check (modified) ---------------------------------------
catalan_N <- factorial(2*N)/(factorial(1+N)*factorial(N))
cat(sprintf("diagrams generated: %d ; Catalan(%d) = %d ; %s\n",
            length(a), N, catalan_N,
            ifelse(length(a)==catalan_N, "PASS", "FAIL")))

# Temperley-Lieb structure (unchanged core) -----------------------------------
Graph_Multiplication_full = function(x,y){
 C = matrix(0L,nrow = 3*number_of_nodes, ncol =3*number_of_nodes)
 interval_1 = c(1:number_of_nodes)
 interval_2 = c((number_of_nodes+1):(2*number_of_nodes))
 interval_3 = c((2*number_of_nodes+1):(3*number_of_nodes))
  C[interval_1,interval_1] = x[interval_1,interval_1 ]
  C[interval_1,interval_2] = x[interval_1,interval_2 ]
  C[interval_2,interval_1] = t(x[interval_1,interval_2 ])
  C[interval_1,interval_3] = t(t(y[interval_1,interval_2])%*%t(x[interval_1,interval_2]))
  C[interval_2,interval_2] = x[interval_2,interval_2]+y[interval_1,interval_1]
  C[interval_2,interval_3] = y[interval_1,interval_2 ]
  C[interval_3,interval_2] = t(y[interval_1,interval_2 ])
  C[interval_3,interval_3] = y[interval_2,interval_2 ]
  for(i in 1:(3*number_of_nodes)){
    C[i,i] = 1
  }
  C
}

Graph_Multiplication =  function(x,y){
  G = t(Graph_Multiplication_full(x,y))
  Final = matrix(0L,nrow = 2*number_of_nodes, ncol =2*number_of_nodes)
  I = diag(1L,nrow = 3*number_of_nodes, ncol =3*number_of_nodes)
  for (i in 1:(number_of_nodes+3)){
    I =  G%*%I
  }
  for(j in 1:(3*number_of_nodes)){
    I[j,j] = 0
  }
  w= union(1:number_of_nodes, (2*number_of_nodes+1):(3*number_of_nodes))
  counter = 1
  for(j in w){
    b = max(union(which(I[w,j] != 0),which(I[j,w] != 0)))
    Final[counter,b] = 1
    counter = counter +1
  }
  for(j in 1:(2*number_of_nodes)){
    Final[j,j] = 0
  }
  Final
}

d = function(x,y){
  G = t(Graph_Multiplication_full(x,y))
  I = diag(1L,nrow = 3*number_of_nodes, ncol =3*number_of_nodes)
  for (i in 1:(number_of_nodes+3)){
    I =  G%*%I
  }
  for(j in 1:(3*number_of_nodes)){
    I[j,j] = 0
  }
  w= union(1:number_of_nodes, (2*number_of_nodes+1):(3*number_of_nodes))
  Sum0 = array()
  Sumn0 = array()
  counter = 1
  counter2 = 1
  for(i in 1:length(I[-w,1])){
    if(sum(I[-w,w][i,])==0){
      Sum0[[counter]] = i
      counter = counter +1
    }else{
      Sumn0[[counter2]] = i
      counter2 = counter2+1
    }
  }
  if(length(Sum0)==1 ){
    if(is.na(Sum0)){
      return(0)
    }
  }
  G = graph_from_adjacency_matrix( I[Sum0+number_of_nodes,-w][,-Sumn0])
  components(G)$no
}

# finding the generators (typo "Genorators" fixed in comment: distinct from the
# seed list `Generators` above)
q = array();for(i in 1:length(a)){ q[[i]]=(sum(a[[i]])) }
a = a[-which(q==N)]
q = array();for(i in 1:length(a)){ q[[i]]=(sum(a[[i]])) }
generators = list()
counter = 1
for (i in 1:length(which(q==(N+2)))){
  G = a[which(q==(N+2))[[i]]][[1]]
  for(j in 1:(number_of_nodes)){
    ct = which(G[j,]==1)
    ct2 = which(G[j+N,]==1)
    if(length(ct2)==1){
      if(ct2 ==(N+j+1)){
        if(length(ct)==1 && ct == j+1){
          if(G[j,ct]-G[ct,j]==0 & G[j+N,ct2]-G[ct2,j+N]==0  ){
            generators[[counter]] = G
            generators = unique(generators)
            counter = counter+1
          }
        }
      }
    }
  }
}

# order diagrams whose bottoms agree (unchanged)
counter = 1
A_ordered_by_bottoms =list()
who_has_been_picked = list()
for(i in 1:length(a)){
  alarm = 0
  if(length(A_ordered_by_bottoms)>=1){
    for(k in 1:length(A_ordered_by_bottoms) ){
      if(A_ordered_by_bottoms[[k]]==i){ alarm = 1 }
    }
  }
  if(alarm ==0){
    hes_is_in = array()
    counter = 1
    for(j in min((i+1),length(a)):length(a)){
      if(max(a[[i]][1:N,1:N] -a[[j]][1:N,1:N])<=0|max(-a[[i]][1:N,1:N] +a[[j]][1:N,1:N])<=0){
        hes_is_in[[counter]] = j
        counter = counter +1
      }
    }
    who_has_been_picked = union(i,hes_is_in)
    A_ordered_by_bottoms[(length(A_ordered_by_bottoms)+1):((length(A_ordered_by_bottoms))+length(who_has_been_picked))] = who_has_been_picked
  }
}
A_ordered_by_bottoms = unique(A_ordered_by_bottoms)

# generator matrices with loop weight 2^d (unchanged core, quiet progress)
BigMatrix = list()
for(i in 1:length(generators)){
  BigMatrix[[i]] = matrix(0L,nrow = length(a), ncol =length(a))
  for(j in 1:length(a)){
    G =  Graph_Multiplication(a[[A_ordered_by_bottoms[[j]]]],generators[[i]])
    for(u in 1:number_of_nodes){
      for(v in (number_of_nodes+1):(2*number_of_nodes)){
        G[v,u] = 0
      }
    }
    for (k in 1:length(a)){
      if(identical(a[[A_ordered_by_bottoms[[k]]]],G)){
        D = d(a[[A_ordered_by_bottoms[[j]]]],generators[[i]])
        BigMatrix[[i]][k,j] = 2^D
      }
    }
  }
}

# --- Temperley-Lieb axioms, collected (modified reporting) --------------------
ax1 <- all(sapply(1:length(generators), function(i)
  identical(BigMatrix[[i]]%*%BigMatrix[[i]], 2*BigMatrix[[i]])))
ax2 <- all(unlist(lapply(1:length(generators), function(j)
  sapply(1:length(generators), function(i)
    if(abs(i-j)>=2) identical(BigMatrix[[i]]%*%BigMatrix[[j]],
                              BigMatrix[[j]]%*%BigMatrix[[i]]) else TRUE))))
ax3 <- all(sapply(1:(length(generators)-1), function(i)
  identical(BigMatrix[[i]]%*%BigMatrix[[i+1]]%*%BigMatrix[[i]], BigMatrix[[i]])))
ax4 <- all(sapply(2:length(generators), function(i)
  identical(BigMatrix[[i]]%*%BigMatrix[[i-1]]%*%BigMatrix[[i]], BigMatrix[[i]])))
cat(sprintf("TL axioms: e^2=2e %s | far-commutation %s | e e' e = e %s %s\n",
            ifelse(ax1,"PASS","FAIL"), ifelse(ax2,"PASS","FAIL"),
            ifelse(ax3,"PASS","FAIL"), ifelse(ax4,"PASS","FAIL")))
