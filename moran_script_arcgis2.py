from pysal.lib import weights
import geopandas as gpd
import pandas as pd
import WeightsUtilities as WU
import numpy as np
from scipy.spatial.distance import cdist
from arcpy.stats import GenerateSpatialWeightsMatrix
from arcpy.stats import SpatialAutocorrelation
import arcpy
import math
import os


def __inverse_weights(gdf,L0,elevation):
    """
    计算反距离权矩阵(私有方法，不应在外部调用)

    参数：
    ------------
    gdf:         shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含权重信息的列表
    """
    # 计算距离矩阵
    if elevation==True:
        coordinates = gdf[['X', 'Y', 'Z']].values
    elif elevation==False:
        coordinates = gdf[['X', 'Y']].values
    distances = cdist(coordinates, coordinates)
    # 应用反距离权重
    spatial_weights = np.where((distances!=0)&(distances <= L0), 1/distances, 0)

    # 生成一个包含权重信息的列表
    weight_info = []
    for i in range(len(spatial_weights)):
        for j in range(len(spatial_weights)):
            weight = spatial_weights[i, j]
            if weight != 0:
                weight_info.append(f"{gdf.iloc[i]['ID']} {gdf.iloc[j]['ID']} {weight}")

    return weight_info




def __gauss_weights(gdf,L0,elevation):
    """
    计算高斯权矩阵(私有方法，不应在外部调用)

    参数：
    ------------
    gdf:    shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含权重信息的列表
    """
    sigma=1.0 #标准方差，控制了函数的曲线在尖峰周围的陡峭程度
    gaussian_const=math.pow(math.pi*2,-0.5) #高斯常数
    # 转为numpy数组
    if elevation==True:
        coordinates = gdf[['X', 'Y', 'Z']].values
    elif elevation==False:
        coordinates = gdf[['X', 'Y']].values
    # 计算距离矩阵
    distances = cdist(coordinates, coordinates)

    # 计算高斯矩阵
    weights1 = np.where(distances<=L0,(np.sqrt(distances/L0)),0)
    weights2 = np.where((distances<=L0)&(distances!=0),(gaussian_const*np.exp(-np.power(weights1,2) / 2)),0)
    # weights = np.where(distances<=L0,(gaussian_const*np.exp(-distances**2 / 2*sigma**2)),0)

    # 生成一个包含权重信息的列表
    weight_info = []
    for i in range(len(weights2)):
        for j in range(len(weights2)):
            weight = weights2[i, j]
            if weight != 0:
                weight_info.append(f"{gdf.iloc[i]['ID']} {gdf.iloc[j]['ID']} {weight}")

    return weight_info



def __threshold_weights(gdf,L0,elevation):
    """
    计算阈值权矩阵(私有方法，不应在外部调用)

    参数：
    ------------
    gdf:    shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含权重信息的列表
    """
    # 提取所需的字段（ID, X, Y, Z）
    if elevation==True:
        coordinates = gdf[['X', 'Y', 'Z']].values
    elif elevation==False:
        coordinates = gdf[['X', 'Y']].values

    #计算要素间的距离
    # coordinates = df[['X', 'Y']].values
    distances = cdist(coordinates, coordinates)

    # 对得到的距离进行判断并赋值
    spatial_weights = np.where((distances<=L0)&(distances!=0), 1, 0)
    # 生成一个包含权重信息的列表
    weight_info = []
    for i in range(len(spatial_weights)):
        for j in range(len(spatial_weights)):
            weight = spatial_weights[i, j]
            if weight != 0:
                weight_info.append(f"{gdf.iloc[i]['ID']} {gdf.iloc[j]['ID']} {weight}")

    return weight_info



def global_moran(shp_path,field,distance_function,generate_report,threshold=float('inf'),std=True,elevation=True):
    """
    主函数，计算Moran'I请调用此函数
    计算全局Moran'I指数，调用该函数后会计算全局Moran'I指数，
    并在指定路径生成txt和png分析结果

    参数：
    -----------
    shp_path:          shapefile文件位置
    field:             输入字段
    distance_function: 空间关系概念化函数，可选threshold/gaussian/inverse
    output_file:       输出路径，输出文件名无需写扩展名
    threshold:         距离阈值，留空则为不设置阈值
    std:               标准化方法，True为行标准化，False为不进行标准化
    elevation:         是否在距离计算中考虑高程影响

    返回值：
    -------------
    None
    """
    arcpy.AddMessage("Reading Shapefile...")
    #获取dataframe
    gdf = gpd.read_file(shp_path)
    #设置当前工作目录
    arcpy.env.workspace=os.getcwd()
    featureset=arcpy.FeatureSet(shp_path)

    arcpy.AddMessage("Calculating Weights Matrix...")
    if distance_function=='threshold':
        spatial_weights_list = __threshold_weights(gdf,threshold,elevation)
    elif distance_function=='gaussian':
        spatial_weights_list = __gauss_weights(gdf,threshold,elevation)  #高斯权函数需要传入一个Shapefile对象
    elif distance_function=='inverse':
        spatial_weights_list = __inverse_weights(gdf,threshold,elevation)
    else:
        raise ValueError("distance_function must be 'threshold', 'gaussian' or 'inverse'")

    # 将权重信息写入到内存空间的txt文件
    txt_file="in_memory/"+os.path.basename(shp_path).replace(".shp",".txt")
    with open(txt_file, 'w',encoding="ascii") as f:
        f.write("ID\n")
        for info in spatial_weights_list:
            f.write(f"{info}\n")

    arcpy.AddMessage("Calculating Moran'I...")
    #计算Moran'I
    SpatialAutocorrelation(featureset,Input_Field=field,Generate_Report=generate_report,Conceptualization_of_Spatial_Relationships="GET_SPATIAL_WEIGHTS_FROM_FILE",
                       Weights_Matrix_File=txt_file,Standardization=std)

    #释放内存
    arcpy.Delete_management(txt_file)


