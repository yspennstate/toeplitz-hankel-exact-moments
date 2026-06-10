library("igraph")

#main functions 
N = 10        #dimentions 
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
plot_braid = function(matrix){
  kappa = 10
  alpha = .5
  theta = 10
  plot(0,0,xlim =c(80,100*N+20),ylim = c(-1,100), cex = .1, axes=FALSE, frame.plot=TRUE)
  for(i in 1:(2*N)){
    x = array()
    y = array()
    if(i<=N){
      if(which(matrix[i,]==1)<=N){
        connected_to = which(matrix[i,]==1)
        x[1:(100*abs(which(matrix[i,]==1)-i))] = 100*seq(from = min(i,which(matrix[i,]==1)),to = max(which(matrix[i,]==1),i), length =(100*abs(which(matrix[i,]==1)-i)) )
        for (t in 1:length(x)){
          y[[t]] =  kappa*(abs(i-(connected_to))^alpha)*(1-(1/(50*abs(i-(which(matrix[i,]==1)) ))^2)*(x[[t]] - .5*(x[[1]] + 100*max(which(matrix[i,]==1),i)))^2)^.5
        }
        
      }else{
        connected_to = which(matrix[i,]==1)
        if(i==(connected_to-N)){
        x[[1]]=100*i
        y[[1]] = 0
        x[[2]]=100*connected_to-100*N
        y[[2]] = 100
        }else{
          if(i<(connected_to-N)){
          x[1:(100*abs(connected_to-N-i))] = 100*seq(from =i,to = connected_to-N, length =(100*abs(connected_to-N-i)) )
          }else{
            x[1:(100*abs(connected_to-N-i))] = 100*seq(from =i,to = connected_to-N, length =(100*abs(connected_to-N-i)) )
            x = x[length(x):1]
          }
          
          k = atan(-50/theta)/(100*i-100*.5*(i+connected_to-N))
          for (t in 1:length(x)){
            y[[t]] =  theta*tan(k*(x[[t]]-100*.5*(i+connected_to-N))) + 50
          }
        }
      }
      
    }else{
      if(which(matrix[,i]==1)<=N){
       #x[[1]]=100*i-100*N
      # y[[1]] = 100
      # x[[2]]=100*which(matrix[,i]==1)
       #y[[2]] = 0
      }else{
        connected_to = which(matrix[,i]==1)
        x[1:(100*abs(which(matrix[i,]==1)-i))] = 100*seq(from = min(i,which(matrix[i,]==1)),to = max(which(matrix[i,]==1),i), length =(100*abs(which(matrix[i,]==1)-i)) )
        x = x - 100*N
        for (t in 1:length(x)){
          y[[t]] =100 - kappa*(abs(i-(connected_to))^alpha)*(1-(1/(50*abs(i-(which(matrix[i,]==1)) ))^2)*(x[[t]] - .5*(x[[1]] + 100* max(which(matrix[i,]==1),i)-100*N))^2)^.5
        }
      }
    }
    lines(x,y, col = rgb(1,0,0), lwd = 2)
  }
  
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
  print(length(a))
  a = unique(a)
  print(length(a))
}

factorial(2*N)/(factorial(1+N)*factorial(N))
length(a)
                    par(mfrow=c(2,5), mai = c(.1, .1, 0.1, 0.1))
                    
                    for(i in sample(c(1:length(a)), 10, replace = FALSE)){
                      plot_braid(a[[i]])
                    }
                    
                    for(i in 1:length(a)){
                      plot_braid(a[[i]])
                    }

#Full multiplication matrix
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
#condenced multiplication
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
#number of loops
d = function(x,y){
  
  
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
    d = 0
    }
    
  }else{
 G = graph_from_adjacency_matrix( I[Sum0+number_of_nodes,-w][,-Sumn0])
  d = components(G)$no
  }
  d
}


#plotting full multiplication 
plot_multiplication = function(u,v){
  kappa = 10
  alpha = .5
  theta = 10
  line_thickness = 2
  matrix = Graph_Multiplication_full(u,v)
  matrix = matrix - diag(1L,nrow = 3*number_of_nodes, ncol =3*number_of_nodes)
 
  plot(-1,0,xlim =c(99,100*N),ylim = c(-1,200), cex = .1, axes=FALSE, frame.plot=TRUE)
  
  for(i in 1:(3*number_of_nodes-1)){
    
    x = array()
    y = array()
    if(i<=N){
      connected_to = min(which(matrix[i,]>=1))
      if(connected_to<=N){
     
        x[1:(100*abs(connected_to-i))] = 100*seq(from = min(i,connected_to),to = max(connected_to,i), length =(100*abs(connected_to-i)) )
        for (t in 1:length(x)){
          y[[t]] = kappa*(abs(i-(connected_to))^alpha)*(1-(1/(50*abs(i-(connected_to) ))^2)*(x[[t]] - .5*(x[[1]] + 100*max(connected_to,i)))^2)^.5
        }
        
      }else{
             
            if(i==(connected_to-N)){
              x[[1]]=100*i
              y[[1]] = 0
              x[[2]]=100*connected_to-100*N
              y[[2]] = 100
            }else{
              if(i<(connected_to-N)){
                x[1:(100*abs(connected_to-N-i))] = 100*seq(from =i,to = connected_to-N, length =(100*abs(connected_to-N-i)) )
              }else{
                x[1:(100*abs(connected_to-N-i))] = 100*seq(from =i,to = connected_to-N, length =(100*abs(connected_to-N-i)) )
                x = x[length(x):1]
              }
              
              k = atan(-50/theta)/(100*i-100*.5*(i+connected_to-N))
              for (t in 1:length(x)){
                y[[t]] = theta*tan(k*(x[[t]]-100*.5*(i+connected_to-N))) + 50
              }
            }
      }
      
    }else{
      
          if(i<=(2*N)){ #second_layer
            connected_to = which(matrix[i,]>=1)
            if(length(connected_to)>=1 & min(connected_to)<i){
            connected_to = connected_to[-c(which(connected_to<i))]
            }
            

          if(length(connected_to)>=1){
            if(max(matrix[i,connected_to])==1){
              for (ct in connected_to){
                        if(ct<=(2*N)){
                          if(u[i,ct]){
                            up_or_down =0
                          }else{
                            up_or_down = 1
                          }
                              if(up_or_down%%2==0){
                              x[1:(100*abs(ct-i))] = 100*seq(from = min(i,ct),to = max(ct,i), length =(100*abs(ct-i)) )
                              x = x - 100*N
                              for (t in 1:length(x)){
                                y[[t]] =100 - kappa*(abs(i-(ct))^alpha)*(1-(1/(50*abs(i-ct ))^2)*(x[[t]] - .5*(x[[1]] + 100*max(ct,i)-100*N))^2)^.5
                              }
                              }else{
                                x[1:(100*abs(ct-i))] = 100*seq(from = min(i,ct),to = max(ct,i), length =(100*abs(ct-i)) )
                                x = x - 100*N
                                for (t in 1:length(x)){
                                  y[[t]] =100 + kappa*(abs(i-(ct))^alpha)*(1-(1/(50*abs(i-ct ))^2)*(x[[t]] - .5*(x[[1]] + 100*max(ct,i)-100*N))^2)^.5
                                }
                              }
                          lines(x,y,col =  rgb(1-up_or_down,0,up_or_down), lwd = line_thickness)
                        }else{
                          x = array()
                          y = array()
                        
                         
                          if(i==(ct-N)){
                            x[[1]]=100*i-N*100
                            y[[1]] = 100
                            x[[2]]=100*ct-100*2*number_of_nodes
                            y[[2]]  = 200
                          }else{
                            if(i<(ct-N)){
                              x[1:(100*abs(ct-N-i))] = 100*seq(from =i-N,to = ct-2*N, length =(100*abs(ct-N-i)) )
                            }else{
                              x[1:(100*abs(ct-N-i))] = 100*seq(from =i-N,to = ct-2*N, length =(100*abs(ct-N-i)) )
                              x = x[length(x):1]
                            }
                            
                            k = atan(-50/theta)/(100*(i-N)-100*.5*(i+ct-3*N))
                            for (t in 1:length(x)){
                              y[[t]] =  theta*tan(k*(x[[t]]-100*.5*(i+ct-3*N))) + 150
                            }
                          }
                          
                          
                          
                          
                          lines(x,y, col = rgb(0,0,1), lwd = line_thickness)
                        }
              }
            }else{
              ct = connected_to
              for(t in 1:2){
                if(t%%2==0){
                  x[1:(100*abs(ct-i))] = 100*seq(from = min(i,ct),to = max(ct,i), length =(100*abs(ct-i)) )
                  x = x - 100*N
                  for (t in 1:length(x)){
                    y[[t]] =100 -kappa*(abs(i-(ct))^alpha)*(1-(1/(50*abs(i-ct ))^2)*(x[[t]] - .5*(x[[1]] + 100*max(ct,i)-100*N))^2)^.5
                  }
                  lines(x,y, col = rgb(1,0,0), lwd = line_thickness)
                }else{
                  x[1:(100*abs(ct-i))] = 100*seq(from = min(i,ct),to = max(ct,i), length =(100*abs(ct-i)) )
                  x = x - 100*N
                  for (t in 1:length(x)){
                    y[[t]] =100 +kappa*(abs(i-(ct))^alpha)*(1-(1/(50*abs(i-ct ))^2)*(x[[t]] - .5*(x[[1]] + 100*max(ct,i)-100*N))^2)^.5
                  }
                  lines(x,y, col = rgb(0,0,1), lwd = line_thickness)
                  
                }
               
              }
            }
          }
          }else{
            if(sum(matrix[i,]) != 0 ){
              connected_to = which(matrix[i,]>=1)
              if(length(connected_to)>=1){
                if(min(connected_to)<i){
                connected_to = connected_to[-which(connected_to<i)]
                }
              }
              if(length(connected_to)>=1){
            ct = connected_to
            x[1:(100*abs(ct-i))] = 100*seq(from = min(i,ct),to = max(ct,i), length =(100*abs(ct-i)) )
            x = x - 200*N
            for (t in 1:length(x)){
              y[[t]] =200 - kappa*(abs(i-(ct))^alpha)*(1-(1/(50*abs(i-ct))^2)*(x[[t]] - .5*(x[[1]] + 100*max(ct,i)-100*2*N))^2)^.5
            }
            lines(x,y, col = rgb(0,0,1), lwd = line_thickness)
              }
            }
          }
    }
    if(length(connected_to)>=1){
   if(i<=N){
    lines(x,y,col =  rgb(1,0,0), lwd = line_thickness)
   }
    }
  }
}

                      
                      
                      Graph_Multiplication_full(a[[1000]],a[[50000]])
                      Graph_Multiplication(a[[100]],a[[50]])
                      A =sample(c(1:length(a)), 1, replace = FALSE)
                      B = sample(c(1:length(a)), 3, replace = FALSE)
                      par(mfrow=c(1,6), mai = c(.1, .1, 0.1, 0.1))
                      for(i in A){
                       
                        for(j in B){
                          i = sample(c(1:length(a)), 1, replace = FALSE)
                      plot_multiplication(a[[i]],a[[j]])
                      plot_braid(Graph_Multiplication(a[[i]],a[[j]]))
                        }
                      }
                      plot_braid(Graph_Multiplication(a[[length(a)]],a[[length(a)]]))
                      d(a[[1000]],a[[50000]])

#finding the generators (diff from Genorators)

q = array();for(i in 1:length(a)){
  q[[i]]=(sum(a[[i]]))
 
}
a = a[-which(q==N)]
q = array();for(i in 1:length(a)){
  q[[i]]=(sum(a[[i]]))
  
}
generators = list()
counter = 1
for (i in 1:length(which(q==(N+2)))){
  G = a[which(q==(N+2))[[i]]][[1]]
  L = array()
  counter2 = 1
  for(j in 1:(number_of_nodes)){
    ct = which(G[j,]==1)
    ct2 = which(G[j+N,]==1)
    if(length(ct2)==1){
    if(ct2 ==(N+j+1)){
      if(ct == j+1){
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
#plotting generators
                      par(mfrow=c(1,4), mai = c(.1, .1, 0.1, 0.1))
                      for (i in 1:length(generators)){
                      plot_braid(generators[[i]])
                      }
                      
                      par(mfrow=c(1,8), mai = c(.1, .1, 0.1, 0.1))
                      for (i in 1:length(generators)){
                        for (j in c(2,7,9,12)){
                          plot_multiplication(a[[A_ordered_by_bottoms[[j]]]],generators[[1]])
                          plot_braid(Graph_Multiplication(a[[A_ordered_by_bottoms[[j]]]],generators[[1]]))
                        }
                      }
                      
 #generators[[length(generators)+1]] = a[which(q==N)][[1]]




#######First we must order those whose bottoms are the same
counter = 1
A_ordered_by_bottoms =list()
who_has_been_picked = list()

for(i in 1:length(a)){
  alarm = 0
  if(length(A_ordered_by_bottoms)>=1){
  for(k in 1:length(A_ordered_by_bottoms) ){
    if(A_ordered_by_bottoms[[k]]==i){
      alarm = 1
    }
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

#plotting ordered
                  par(mfrow=c(3,5), mai = c(.1, .1, 0.1, 0.1))
                for(i in 1:length(a)){
                  plot_braid(a[[A_ordered_by_bottoms[[i]]]])
                  }



##################
BigMatrix = list()
n <- length(generators)
pb <- txtProgressBar(min = 0, max = n, style=3)
for(i in 1:length(generators)){
  BigMatrix[[i]] = matrix(0L,nrow = length(a), ncol =length(a))
  RowNames = list()
  ColNames = list()
  for (n in 1:length(a)){
    RowNames[[n]] =paste("S",n)
  }
  for(n in 1:length(a)){
    ColNames[[n]] = paste("S",n) 
  }
  row.names(BigMatrix[[i]]) = RowNames
  colnames(BigMatrix[[i]]) = ColNames
  row =array()
  counter = 1
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
      
        row[[counter]]= k
        counter = counter +1
      }
    }
  }
  row = unique(row)
  #BigMatrix[[i]]= BigMatrix[[i]][row[order(row)],]
  setTxtProgressBar(pb, i)
}
BigMatrix[[1]]

BigMatrix[[1]]%*%BigMatrix[[1]]





#These are the axioms 


for(i in  1:length(generators)){
  if(identical(BigMatrix[[i]]%*%BigMatrix[[i]],2*BigMatrix[[i]])){
     print ("yes")}else{print ("no")}
  
}
  
  


for(j in  1:length(generators)){

for(i in  1:length(generators)){if(abs(i-j)>=2){
    if(identical(BigMatrix[[i]]%*%BigMatrix[[j]],BigMatrix[[j]]%*%BigMatrix[[i]])){
      print ("yes")}else{print ("no")}}
                                                          }
}




for(i in  1:(length(generators)-1)){
  if(identical(BigMatrix[[i]]%*%BigMatrix[[i+1]]%*%BigMatrix[[i]],BigMatrix[[i]])){
    print ("yes")}else{print ("no")}
  
}

for(i in  2:length(generators)){
  if(identical(BigMatrix[[i]]%*%BigMatrix[[i-1]]%*%BigMatrix[[i]],BigMatrix[[i]])){
    print ("yes")}else{print ("no")}
  
}




