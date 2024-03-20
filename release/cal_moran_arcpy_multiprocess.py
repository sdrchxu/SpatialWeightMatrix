import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from arcpy.stats import SpatialAutocorrelation
import arcpy
import math
import os
from multiprocessing import Pool, cpu_count

## 该脚本用于调用arcpy输出Moran'I指数报告，支持多线程运算 ##
## 目前在命令行中测试通过，GeoScene Pro 2.1导入测试不通过 ##

def process_weights(args):
    """
    将包含权重的numpy数组转为列表
    """
    i, spatial_weights_length, spatial_weights, gdf = args
    weight_info = []
    for j in range(spatial_weights_length):
        weight = spatial_weights[i, j]
        if weight != 0:
            weight_info.append(f"{int(gdf.iloc[i]['ID'])} {int(gdf.iloc[j]['ID'])} {weight}")
    return weight_info



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
    spatial_weights_length=len(spatial_weights)
    #按CPU核心数分配进程
    with Pool(cpu_count()) as p:
        args = [(i, spatial_weights_length, spatial_weights, gdf) for i in range(spatial_weights_length)]
        weight_info = p.map(process_weights, args)
    #将嵌套列表展开
    weight_info = [item for sublist in weight_info for item in sublist]

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
    spatial_weights_length=len(weights2)
    #按CPU核心数分配进程
    with Pool(cpu_count()) as p:
        args = [(i, spatial_weights_length, weights2, gdf) for i in range(spatial_weights_length)]
        weight_info = p.map(process_weights, args)
    #将嵌套列表展开
    weight_info = [item for sublist in weight_info for item in sublist]

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
    spatial_weights_length=len(spatial_weights)
    #按CPU核心数分配进程
    with Pool(cpu_count()) as p:
        args = [(i, spatial_weights_length, spatial_weights, gdf) for i in range(spatial_weights_length)]
        weight_info = p.map(process_weights, args)
    #将嵌套列表展开
    weight_info = [item for sublist in weight_info for item in sublist]

    return weight_info


def __read_features_to_dataframe(shp_path, z_field, id_field):
    """
    从shapefile中获取X,Y,Z,ID(私有方法，不应在外部调用)   

    参数：
    -------------
    shp_path: shp文件路径
    z_field:  z字段名
    id_field: id字段名   

    返回值：
    -------------
    DataFrame，其中columns=['X', 'Y', 'Z', 'ID']

    """
    data = []
    with arcpy.da.SearchCursor(shp_path, ["SHAPE@XY", z_field, id_field]) as cursor:
        for row in cursor:
            x, y = row[0]
            z = row[1]
            id_value = row[2]
            data.append((x, y, z, id_value))
    df = pd.DataFrame(data, columns=['X', 'Y', 'Z', 'ID'])
    return df


def __read_features_to_dataframe_noele(shp_path,id_field):
    """
    从shapefile中获取X,Y,Z,ID(私有方法，不应在外部调用)   

    参数：
    -------------
    shp_path: shp文件路径
    z_field:  z字段名
    id_field: id字段名   

    返回值：
    -------------
    DataFrame，其中columns=['X', 'Y','ID']

    """
    data = []
    with arcpy.da.SearchCursor(shp_path, ["SHAPE@XY", id_field]) as cursor:
        for row in cursor:
            x, y = row[0]
            id_value = row[1]
            data.append((x, y,id_value))
    df = pd.DataFrame(data, columns=['X', 'Y', 'ID'])
    return df




def global_moran(shp_path,analyze_field,z_field,id_field,distance_function,generate_report,
                threshold=float('inf'),std='None',elevation=False):
    """
    主函数，计算Moran'I请调用此函数
    计算全局Moran'I指数，调用该函数后会计算全局Moran'I指数，
    并返回PeoProcessing工具执行结果

    参数：
    -----------
    shp_path:          shapefile文件位置
    analyze_field:     用于评估空间自相关的字段
    z_field:           高程字段
    id_field:          唯一ID字段
    distance_function: 空间关系概念化函数，可选threshold/gaussian/inverse
    generate_report:   是否输出报告
    threshold:         距离阈值，留空则为不设置阈值
    std:               标准化方法，"Row"为行标准化，"None"为不进行标准化
    elevation:         是否在距离计算中考虑高程影响

    返回值：
    -------------
    全局莫兰指数执行结果
    """
    # arcpy.AddMessage("Reading Shapefile...")
    arcpy.SetProgressorLabel("Reading Shapefile...")
    arcpy.SetProgressorPosition(20)
    # print("Reading Shapefile...")
    #获取dataframe
    if elevation==True:
        gdf = __read_features_to_dataframe(shp_path,z_field,id_field)
    else:
        gdf=__read_features_to_dataframe_noele(shp_path,id_field)

    #设置当前工作目录
    arcpy.env.workspace=os.getcwd()
    featureset=arcpy.FeatureSet(shp_path)

    arcpy.SetProgressorLabel("Calculating Weights Matrix...")
    arcpy.SetProgressorPosition(40)
    # print("Calculating Weights Matrix...")
    if distance_function=='threshold':
        spatial_weights_list = __threshold_weights(gdf,threshold,elevation)
    elif distance_function=='gaussian':
        spatial_weights_list = __gauss_weights(gdf,threshold,elevation)  #高斯权函数需要传入一个Shapefile对象
    elif distance_function=='inverse':
        spatial_weights_list = __inverse_weights(gdf,threshold,elevation)
    else:
        raise ValueError("distance_function must be 'threshold', 'gaussian' or 'inverse'")

    # print("Generate Weight File...")
    arcpy.SetProgressorLabel("Generate Weight File...")
    arcpy.SetProgressorPosition(60)
    # 将权重信息写入到系统临时文件夹的txt文件
    tmp_dir=os.getenv('TEMP')
    txt_file=tmp_dir+"\\"+os.path.basename(shp_path).replace(".shp",".txt")
    with open(txt_file, 'w',encoding="ascii") as f:
        f.write(id_field+"\n")
        for info in spatial_weights_list:
            f.write(f"{info}\n")

    arcpy.SetProgressorLabel("Calculating Moran'I...")
    arcpy.SetProgressorPosition(80)
    # arcpy.AddMessage("Calculating Moran'I...")
    # print("Calculating Moran'I...")

    #计算Moran'I
    result=SpatialAutocorrelation(featureset,Input_Field=analyze_field,Generate_Report=generate_report,Conceptualization_of_Spatial_Relationships="GET_SPATIAL_WEIGHTS_FROM_FILE",
                    Weights_Matrix_File=txt_file,Standardization=std)
    
    for msg in range(0,result.messageCount):
        arcpy.AddReturnMessage(msg)
    #删除临时文件
    arcpy.Delete_management(txt_file)
    arcpy.SetProgressorPosition(100)
    return result


if __name__ == "__main__":
    #控制进度条
    arcpy.SetProgressor("step","Loading Script...")

    #获取输入参数
    shp_path=arcpy.GetParameterAsText(0)
    analyze_field=arcpy.GetParameterAsText(1)
    z_field=arcpy.GetParameterAsText(2)
    id_field=arcpy.GetParameterAsText(3)
    distance_function=arcpy.GetParameterAsText(4)
    generate_report=arcpy.GetParameter(5)
    threshold=arcpy.GetParameter(6)
    std=arcpy.GetParameterAsText(7)
    elevation=arcpy.GetParameter(8)


    result=global_moran(shp_path,analyze_field,z_field,id_field,distance_function,generate_report,threshold,std,elevation)
    #设置输出消息
    # arcpy.AddMessage(result.getMessages())
    


