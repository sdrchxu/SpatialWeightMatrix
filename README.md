# SpatialWeightMatrix
## 项目描述
一个能够便捷地使用要素Z值字段计算空间权重矩阵并进行空间自相关分析的工具   
## 特性  
- 支持考虑高程的空间自相关计算
- 可选基于距离的二进制加权、反距离权、高斯权三种空间权重矩阵生成方式
- 支持多shp文件批处理，并创建Global Moran'I分析汇总表
- 权重文件生成过程支持并行计算
- 兼容Arcgis, Arcgis Pro, GeoScene Pro, Geoda
## 使用说明
需要python3解释器，建议python 3.7及以上。    
**请使用release文件夹中的python脚本**
### 使用"cal_moran_pysal"
1. 安装依赖
~~~
pip install geopandas
pip install pandas
pip install numpy
pip install pysal
pip install scipy
pip install splot
pip install matplotlib
~~~
2. 调用模块
~~~python
import sys
sys.path.append(r"{存放cal_moran_pysal.py文件的路径}")
import cal_moran_pysal
~~~
3. 计算Global Moran'I
若要执行对单个要素的空间自相关分析：
~~~python
global_moran(shp_path,field,output_file,distance_function,threshold,std,elevation)
~~~
若要对某个文件夹中的全部空间要素进行分析：
~~~python
global_moran_folder(folder_path, field, output_folder, distance_function, threshold, std, elevation)
~~~
### 使用cal_moran_arcpy(_multiprocess)
**必须为Windows操作系统，需要Arcpy模块**
该脚本可导入Arcgis或GeoScene系列软件的自定义工具箱中使用，获得与Arcgis原生工具箱相似的体验。    
请优先使用带"multiprocess"后缀的脚本，该脚本支持了多核并行计算，目前已知此脚本在部分GeoScene Pro的图形化界面中可能报错：“pickle.PicklingError: Can't pickle <functionprocess weights at0x0000022901A3A168>: attribute lookup   process weights on  mainfailed”。若出现此问题，请在命令行中使用脚本或在图形化界面中使用不带"multiprocess"后缀的脚本。
1. 


