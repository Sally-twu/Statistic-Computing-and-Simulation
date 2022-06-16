Sys.time()
set.seed(12345)
#easy
# s <- matrix(0,ncol=9,nrow=9)
# s[1,c(1,4,5,8,9)]=c(4,7,6,1,2)
# s[2,c(1:3,5,6,8,9)]=c(7,8,2,1,5,4,9)
# s[3,c(2:8)]=c(6,1,4,2,8,3,7)
# s[4,c(1,2,6:9)]=c(5,1,4,7,6,3)
# s[5,c(1,8,9)]=c(3,2,8)
# s[6,c(3:5)]=c(6,9,3)
# s[7,c(1,5)]=c(2,4)
# s[8,c(1:9)]=c(1,4,7,5,9,3,2,8,6)
# s[9,c(1)]=c(6)
#medium
# s <- matrix(0,ncol=9,nrow=9)
# s[1,c(6,7)]=c(2,1)
# s[2,c(1,3,5,7,8)]=c(2,3,1,6,9)
# s[3,c(4,5,8)]=c(4,3,8)
# s[4,c(3,4,6,8)]=c(1,5,8,4)
# s[5,c(3,4,6:8)]=c(9,7,1,5,6)
# s[6,c(1:4,6:8)]=c(5,7,2,6,3,9,1)
# s[7,c(3,6,7,9)]=c(8,4,3,6)
# s[8,c(1:4,7:9)]=c(1,2,4,3,8,5,9)
# s[9,c(2,6,7)]=c(6,9,4)
# #hard
s <- matrix(0,ncol=9,nrow=9)
s[1,c(3,8)]=c(7,8)
s[2,c(3,5,7,9)]=c(8,2,9,5)
s[3,c(2,5,7)]=c(9,5,1)
s[4,c(4:8)]=c(4,3,5,7,2)
s[5,c(1,3,5,7,8)]=c(7,3,6,4,9)
s[6,c(2,7:9)]=c(2,5,3,6)
s[7,c(3,4,6,7)]=c(5,6,4,2)
s[8,c(1,8)]=c(2,4)
s[9,c(4,6,7)]=c(2,3,8)
#計算重複值
target <- function(s){
  tar <- sum(apply(s,1,duplicated)+apply(s,2,duplicated))#前:對列回傳是否為重複值
  for(r in 1:9){#3*3九宮格是否有重複值
    bloa <- (1:3)+3*(r-1)%%3
    blob <- (1:3)+3*trunc((r-1)/3)
    tar <- tar+sum(duplicated(as.vector(s[bloa,blob]))) 
  }
  return(tar)
}
#隨機塞值
pool=array(TRUE,dim=c(9,9,9))
for (i in 1:9)
  for (j in 1:9){
    if (s[i,j]>0) pool[i,j,-s[i,j]]=FALSE
  }
for(t in 1:100){# random order for visit of all sites
  for(i in sample(1:81)){
    if(s[i]==0){
      a=((i-1)%%9)+1
      b=trunc((i-1)/9)+1
      boxa=3*trunc((a-1)/3)+1
      boxa=boxa:(boxa+2)
      boxb=3*trunc((b-1)/3)+1
      boxb=boxb:(boxb+2)
      for (u in (1:9)[pool[a,b,]]){#eliminates impossible values
        pool[a,b,u]=(sum(u==s[a,])+sum(u==s[,b])
                     +sum(u==s[boxa,boxb]))==0
      }
      if (sum(pool[a,b,])==1){ # only one possible case, solution found!
        s[i]=(1:9)[pool[a,b,]]
        #print(s)
      }
      if (sum(pool[a,b,])==0){ # solution does not exist, exit!
        print("wrong sudoku")#sum(pool[a,b,])==0可以放0種數字
        break()
      }
    }
  }
}

Niter <- 10^5
Nsteps <- 10^1#影響其收斂效果
lmax <- 10^5
temps <- exp(sqrt(seq(1,log(lmax)^2,le=Niter+1)))#le:切幾分
lcur <- temps[1]#2.718282e+00
cur <- s
for(r in (1:81)[s==0]){#隨機填數字
  cur[r]=sample((1:9)[pool[r+ 81*(0:8)]],1)
}
tarcur <- target(cur)#重複的數字有幾個

plot(c(0,0),col="white",xlim=c(1,Niter),ylim=c(0,tarcur),xlab="iterations",ylab="penalty")

for(t in 1:Niter){
  if (tarcur==0){
    print(t)
    print(cur)
    break()
  }
  nchange <- 0
  for(d in 1:Nsteps){
    prop <- cur
    i=sample((1:81)[as.vector(s)==0],sample(1:sum(s==0),1,pro=1/(1:sum(s==0))))#連抽幾個都是隨機抽樣
    for (r in 1:length(i))
      prop[i[r]]=sample((1:9)[pool[i[r]+81*(0:8)]],1)
    if (log(runif(1))/lcur<tarcur-target(prop)){#如果新的比較好就接受新解
      #如果不除lcur代表他可以接受比較壞的結果
      nchange=nchange+(tarcur>target(prop))
      cur=prop#新的解
      points(t,tarcur,col="forestgreen",cex=.3,pch=19)
      tarcur=target(cur)#新的tarcur
    }
    if (tarcur==0){
      break()
    }
    lcur=sample(c(1,10^(-4)),1,pro=c(1.5*(log(t+1))+1,1))*temps[t+1]#15176.3
  }
}
Sys.time()
