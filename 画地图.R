
rm(list = ls(all=TRUE)) 
library(sp)
library(rgeos)
library(maptools)
library(raster)
library(rgdal)
library(rjson)
library(ggplot2)

### 读入shapefile数据
town<- readOGR(dsn = "C:/Users/lenovo/Desktop/mapdata", layer="TOWN_MOI_1070205", use_iconv = TRUE, encoding = "UTF-8") ### 解决了中文乱码问题
plot(town)

###另两种打开文件的方式
##1. crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")    
## states=readShapePoly("mapdata201803011040/TOWN_MOI_1070205.shp",proj4string=crswgs84,verbose=TRUE)  
##2. states=readShapePoly(file.choose())  
a<-town@data[,c(3,4)]   ######## 观察与核对地区名
#write.csv(a,"E:/台湾空气污染/毕业论文/地图地名核对.csv")
#town_name<-write.csv(a,"town_name.csv")
#town[4,]   ## 判断是否是离岛的地区
out_island<-c(4,5,26,27,28,29,30,31,132,148,149,150,151,152,153,362,363,364,365) ## 去掉离岛地区
town<-town[-out_island,]   ##### 去掉19个离岛地区
str(town@data)

### 计算邻接矩阵
adj.mat.sub <- gTouches(town, byid=TRUE, returnDense=TRUE)
adj.mat.sub[adj.mat.sub==TRUE]<-1
adj.mat.sub[adj.mat.sub==FALSE]<-0
adj.mat.sub
#write.table(adj.mat.sub,"邻接矩阵.txt",col.names = FALSE, row.names = FALSE)


### 判断是否还有离岛地区
sum1<-matrix(0,nrow=349,ncol = 1)
for (i in 1:349){
  sum1[i]<-sum(adj.mat.sub[i,])
  
}
which(sum1==0)
plot(town)


### 画disease mapping 发病率的统计地图
### 1.将估计的p与
p<-read.csv("C:/Users/lenovo/Desktop/diseasemappingresult/2013/result/Rprob.csv",header=TRUE)
p<-p[,2]
p<-round(p,3)
p
ID<-town@data[,"TOWNCODE"]
ID<-as.character(ID)
ID<-as.numeric(ID)
p1<-cbind(ID,p) ### 只要有两列id和数值
### 2.将数值与空间数据连接起来
p1<-as.data.frame(p1)
names(p1)<-c("TOWNCODE","p")
sh2<-merge(town@data,p1,by=intersect(names(town@data),names(p1)),sort=FALSE)


# Set the palette
p <- colorRampPalette(c("white", "red"))(7)
palette(p)

# Scale the total population to the palette
pop <-sh2$p
pop
min(pop)
max(pop)
diff(range(pop))

cols<-matrix(NA,length(pop),1)
#cols[which(pop==0),]<-1
cols[which(pop>=0),]<-1
cols[which(pop>=0.2),]<-2
cols[which(pop>=0.25),]<-3
cols[which(pop>=0.30),]<-4
cols[which(pop>=0.35),]<-5

table(cols)
plot(town, col=cols)
cols
sh2

## 加图例
legend("right", legend = c("<0","0-0.2","0.25-0.30","0.30-0.35",">=0.35"),fill = 1:5)


