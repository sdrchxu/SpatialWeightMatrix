import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import math
import os
import arcpy
from multiprocess import Pool, cpu_count

## 该脚本用于多线程生成适于Arcgis或GeoDa空间自相关分析的空间权重矩阵文件 ##

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


def cal_weight_txt(shp_path,out_txt_path,z_field,id_field,distance_function,
                threshold=float('inf'),elevation=False,software='arcgis'):
    """
    主函数，计算Arcgis权重矩阵文件请调用此函数
    计算全局Moran'I指数，调用该函数后会计算全局Moran'I指数，
    并在指定路径生成txt和png分析结果

    参数：
    -----------
    shp_path:          shapefile文件位置
    out_txt_path:      输出的txt权重文件路径
    z_field:           高程字段
    id_field:          唯一ID字段
    distance_function: 空间关系概念化函数，可选threshold/gaussian/inverse
    threshold:         距离阈值，留空则为不设置阈值
    elevation:         是否在距离计算中考虑高程影响

    返回值：
    -------------
    无
    """
    # arcpy.AddMessage("Reading Shapefile...")
    print("Reading Shapefile...")
    #获取dataframe
    if elevation==True:
        gdf = __read_features_to_dataframe(shp_path,z_field,id_field)
    else:
        gdf=__read_features_to_dataframe_noele(shp_path,id_field)

    #设置当前工作目录
    arcpy.env.workspace=os.getcwd()
    featureset=arcpy.FeatureSet(shp_path)

    # arcpy.AddMessage("Calculating Weights Matrix...")
    print("Calculating Weights Matrix...")
    if distance_function=='threshold':
        spatial_weights_list = __threshold_weights(gdf,threshold,elevation)
    elif distance_function=='gaussian':
        spatial_weights_list = __gauss_weights(gdf,threshold,elevation)  #高斯权函数需要传入一个Shapefile对象
    elif distance_function=='inverse':
        spatial_weights_list = __inverse_weights(gdf,threshold,elevation)
    else:
        raise ValueError("distance_function must be 'threshold', 'gaussian' or 'inverse'")

    print("Generate Weight File...")
    # 将权重信息写入到txt文件
    if software=='arcgis':
        with open(out_txt_path, 'w',encoding="ascii") as f:
            f.write("ID\n")
            for info in spatial_weights_list:
                f.write(f"{info}\n")
    elif software=='geoda':

    print("Done!")

def cal_weight_txt_folder(shp_folder,output_folder,z_field,id_field,distance_function,
                threshold=float('inf'),elevation=False):
    """
    计算给定文件夹中所有shapefile文件的Moran'I指数

    参数：
    -----------
    shp_folder:        shapefile文件夹路径
    output_folder:     结果输出文件夹路径
    z_field:           高程字段
    id_field:          唯一ID字段
    distance_function: 空间关系概念化函数，可选threshold/gaussian/inverse
    threshold:         距离阈值，留空则为不设置阈值
    elevation:         是否在距离计算中考虑高程影响

    返回值：
    -------------
    无
    """
    for file_name in os.listdir(shp_folder):
        if file_name.endswith(".shp"):
            shp_path = os.path.join(shp_folder, file_name)
            output_file = os.path.join(output_folder, file_name.replace(".shp",".txt"))
            print(f"Processing {file_name}...")
            cal_weight_txt(shp_path,output_file,z_field,id_field,distance_function,
                threshold,elevation)




if __name__ == "__main__":
    
    cal_weight_txt("F:\大创数据\中间产出的数据\云南省和黄淮海平原已处理好的火点\云南省逐月火点\云南省已处理好的火点_1月.shp",
                   "D:\Lenovo\Desktop\云南大学\大创\程序代码\空间权重矩阵测试\swm测试\云南省已处理好的火点_1月.txt",'Z','ID','gaussian',55000,
                   True)

