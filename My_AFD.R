My_AFD<-function(data_table,Chien,Loup,genre,names,centree=T,Reduit=T) {#,genre,names)
                
  
  #variables globales
  
  
  n<-nrow(data_table)
  p<-ncol(data_table)
  n1<-nrow(Chien)
  n2<-nrow(Loup)
  
  nbr_classes<-2
  
  g<-colMeans(data_table)
  g1<-colMeans(Chien)
  g2<-colMeans(Loup)
  
  if(centree)      #Centralisation des variables
  {  
    
    data_table<-data_table-t(replicate(n,g))
    
    Chien<-Chien-t(replicate(n1,g1))
    
    Loup<-Loup-t(replicate(n2,g2))
  }
  
  if(Reduit)          #Réduction  des données
  {
    
    
    sd_vector<-GMCM:::colSds(data_table)
    data_table<-data_table/(t(replicate(n,sd_vector)))
    
    sd_vector_C<-GMCM:::colSds(Chien)
    Chien<-Chien/(t(replicate(n1,sd_vector_C)))
    
    sd_vector_L<-GMCM:::colSds(Loup)
    Loup<-Loup/(t(replicate(n2,sd_vector_L)))
    
  }
  
  Chien<-data_table[1:n1,]
  Loup<-data_table[(n1+1):n,]
  
  g<-colMeans(data_table)
  g1<-colMeans(Chien)
  g2<-colMeans(Loup)
  
  
  
  # La matrice de variance-covariance V 
  V<-matrix(0,nrow=p,ncol=p)
  for ( i in 1:n)
  {
    V<-V+(data_table[i,]-g)%*%t(data_table[i,]-g)
    
  }
  V<-V/n
  
  #La matrice inter-classes B (Between)
  nk<-c(n1,n2)
  B<-(g1-g)%*%t(g1-g)*n1+(g2-g)%*%t(g2-g)*n2
  B<-B/n
  
  # Les matrices intra-classes 
  W1<-cov(Chien)*(n1-1)/n1  # on utilise cette astuce parce que  
  W2<-cov(Loup)*(n2-1)/n2   #la covariance sous R est non biasée.
  
  W<-n1*W1+n2*W2
  W<-W/n
  
  
  #compariason des variances
  print(V>B)
  print("les points dans l'espaces sont plus dispersées que les centres des classes")
  
  print(V>W)
  print("les points dans l'espaces sont plus dispersées que les points de chaque classe")
  
  
  #Les variables initiales partagent le groupe.
  
  print("les variables initiales partagent le groupe")
  print(V/(B+W))
  
  
  #---------Question 2: étude de V-1 et B -----------------------------------------
  
  # La matrice B est une matrice symétrique  positive ==> Diagonalisable 
  #et de dimension p*p
  # La matrice V est une matrice symétrique  positive ==> Diagonalisable
  #et de dimension p*p
  
  V_inverse<-solve(V) # la matrice V_(-1) est une matrice symétrique
  #(car V est symetrique ) et de dimension p*p
  
  
  d<-V_inverse%*%B    # La matrice n'est pas forcément diagonalisable :
  # On peut pas la diagonaliser directement
  
  #----------------Question 3: Représentation de chaque groupe d'individu--------
  
  #Diagonalisation de la matrice (V-1B)
  #Passage par l'intermidiare C
  
  C<-matrix(0,nrow=p,ncol=nbr_classes)
  C[,1] <-(g1-g)*sqrt(n1/n)
  C[,2]<-(g2-g)*sqrt(n2/n)
  
  print("la décompsition de B")
  
  print(B/(C%*%t(C)))
  
  #on remarque bien que B=C%*%C'
  
  eigen_matrix<-eigen(t(C)%*%V_inverse%*%C)
  
  
  lambda<-eigen_matrix$values
  barplot(lambda)
  vectorss<-eigen_matrix$vectors
  vectors<-V_inverse%*%C%*%vectorss
  
  
  #Les axes factoriels
  F1<-data_table%*%vectors[,1]
  F2<-data_table%*%vectors[,2]
  F1C<-F1[1:n1]
  F2C<-F2[1:n1]
  
  
  #-----------------Question 6: rapport de corrélation----------------------------
  #On prend  un seul axe factoriel car on a que deux classes (Chiens et Loups) donc:
  
  
  
  Z<-data_table%*%vectors[,1]
  SCT<-sum((Z-mean(Z))*(Z-mean(Z)))
  SCE<-n1*(mean(Z[1:n1])-mean(Z))*(mean(Z[1:n1])-mean(Z))
  SCE<-SCE+ n2*(mean(Z[(n1+1):n])-mean(Z))*(mean(Z[(n1+1):n])-mean(Z))
  mu1 <-SCE/SCT
  
  
  Z<-data_table%*%vectors[,2]
  SCT<-sum((Z-mean(Z))*(Z-mean(Z)))
  SCE<-n1*(mean(Z[1:n1])-mean(Z))*(mean(Z[1:n1])-mean(Z))
  SCE<-SCE+ n2*(mean(Z[(n1+1):n])-mean(Z))*(mean(Z[(n1+1):n])-mean(Z))
  mu2 <-SCE/SCT
  
  plot(F1C,F2C,pch=15,col="red",
       xlab=paste("Axe 1(",round(mu1*100,2),"%)"),
       ylab=paste("Axe 2(",round(mu2*100,2),"%)"),
       
       xlim=c(min(F1)-sd(F1)/4,max(F1)+ sd(F1)/4),
       ylim = c(min(F2)-sd(F2)/4,max(F2)+sd(F2)/4))
  F1L<-F1[(n1+1):n]
  F2L<-F2[(n1+1):n]
  points(F1L,F2L,pch=18,col="blue")
  abline(v=(min(F1L)+max(F1C))/2,lty=2)
  title("Discirmination des classes Chien et Loup")
  text(F1, y = F2, labels = genre, offset = 0.25,cex = 1)
  
  
  if(centree & Reduit)
  {
    library("plotrix")
    r<-matrix(0,nrow=p,ncol=2)
   
    r[,1]<-sqrt(lambda[1])*vectors[,1]
    r[,2]<-sqrt(abs(lambda[2]))*vectors[,2]
    
    plot.new()
    plot.window(xlim=c(-1.2,1.2),ylim=c(-1.2,1.2))
   # plot.window(#xlab=paste("Axe 1(",round(mu1*100,2),"%)"),
                #ylab=paste("Axe 2(",round(mu2*100,2),"%)"),
    #            xlim=c(-1.2,1.2),
     #           ylim=c(-1.2,1.2)#c(-max(abs(r[,2]))-abs(sd(r[,2])),max(abs(r[,2]))+abs(sd(r[,2])))
                
                
       #         )
                
    
    abline(h=0,lty=2)
    abline(v=0,lty=2)
    arrows(integer(6),integer(6),r[,1],1e+23*r[,2],lwd=2)
    draw.circle(0,0,1)
    title("Graph of quantitative variables")
    #text(r[,1], y = 1e+23*r[,2], labels = names, offset = 0.3,cex = 1)
    
    
  }
  
  return(r)
  
}
