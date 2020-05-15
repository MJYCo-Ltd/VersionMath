TEMPLATE = subdirs

CONFIG += ordered

### 为整个STK的支撑项
SUBDIRS += Sofa \           # IAU发布的基础算法
           Math \
           Satellite        # 卫星轨道等算法
### end
