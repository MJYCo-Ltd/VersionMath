TEMPLATE = subdirs

CONFIG += ordered

### 为整个STK的支撑项
SUBDIRS +=  Math \
           GisMath \
           Satellite\        # 卫星轨道等算法
           SatelliteToolKit\ #卫星工具类
           TestSTK
### end
