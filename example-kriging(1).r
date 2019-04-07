
#### Krig准备工作开始, 导入程序包
library(sp)
library(gstat)

#### simdat需要转换为sp格式, 需要指定坐标
data=read.csv("C:/Users/lenovo/Desktop/diseasemappingresult/2013/data/monitoring point.csv",header=FALSE)
data1=read.csv("C:/Users/lenovo/Desktop/diseasemappingresult/2013/data/disease point.csv",header=FALSE)

coordinates(data) <- ~V1+V2
coordinates(data1) <-~V1+V2
spplot(data,"V3") ### 查看每个点的值   


#### 创建克里格所用的网格, 20m方格的中心位置
xgrid <- seq(ceiling(min(data$V1)), floor(max(data$V1)), by = 0.005) 
ygrid <- seq(ceiling(min(data$V2)), floor(max(data$V2)), by = 0.005) 
basexy <- expand.grid(xgrid, ygrid)
colnames(basexy) <- c("x", "y")
plot(y ~ x, basexy)

coordinates(basexy) <- ~x+y
gridded(basexy) <- TRUE

### 因Kriging需要基于 Variogram以及相应的模型.Variogram横轴是距离, 纵轴是方差. 
#做模型拟合时, 应先查看 Variogram. Variogram根据点与点之间的距离, 分成若干段, 并计算每一段的方差, 
#再通过拟合Variogram模型, 获得Krig所需要的参数. 拟合Variogram模型时, 需要提供模型的初始值,再进行迭代优化, 否则很难获得模型的准确参数. 

######### 若所进行的是普通克里格插值, 则需要插值的变量只与自身取样点的距离相关,
#用公式表示为 alt~1, 这里的~1是variogram的规定 (参见 ?variogram). 
vgm1 <- variogram(V3~1, data)
plot(vgm1) 

### variogram的模型常包括
### Circular
### Spherical
### Exponential
### Gaussian
### Linear
### gstat程序包提供的模型可以通过 vgm() 查看

m <- fit.variogram(vgm1,vgm(80,"Sph",1.5,0.2))
plot(vgm1, model=m)

### 在拟合gm模型时, 需要提供以下几个参数的初始值
##psill: partial sill
##model: 模型的类型
##range: 范围
###nugget: 块金值


###各参数的几何意义如图一所示. 初始值要通过在variogram的图中读取, 尽量接近真实值, 这样后续的迭代才会成功. 

##例如以下例子中: .01 为 psill
##"Sph" 为球状模型
##300 为范围
##nugget为 块金值. 

### 设定vgm模型的初始值
m <- vgm(80,"Sph",1.5, 0.2)

### 进行克里格
krige_res <- krige(V3~1, data, basexy, model = m)

### 查看克里格插值的结果
spplot(krige_res, zcol = "var1.pred", main = "Predictions of altitude based on the randomly sampled data", col.regions = terrain.colors(100))
aa<- krige(V3~1, locations= data, newdata= data1, m <- vgm(80,"Sph",1.5, 0.2) )
aa
write.csv(aa,"C:/Users/lenovo/Desktop/diseasemappingresult/2013/data/pollute_after.csv")
