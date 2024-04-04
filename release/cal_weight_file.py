import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import math
import os
import arcpy
from multiprocessing import Pool, cpu_count
import gc

## 该脚本用于批量或非批量生成适于Arcgis或GeoDa空间自相关分析的空间权重矩阵文件,支持多线程加速 ##

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


def __inverse_weights(gdf,L0,elevation,thread_num):
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
    with Pool(thread_num) as p:
        args = [(i, spatial_weights_length, spatial_weights, gdf) for i in range(spatial_weights_length)]
        weight_info = p.map(process_weights, args)
    #将嵌套列表展开
    weight_info = [item for sublist in weight_info for item in sublist]

    return weight_info




def __gauss_weights(gdf,L0,elevation,thread_num):
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
    with Pool(thread_num) as p:
        args = [(i, spatial_weights_length, weights2, gdf) for i in range(spatial_weights_length)]
        weight_info = p.map(process_weights, args)
    #将嵌套列表展开
    weight_info = [item for sublist in weight_info for item in sublist]

    return weight_info



def __threshold_weights(gdf,L0,elevation,thread_num):
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
    with Pool(thread_num) as p:
        args = [(i, spatial_weights_length, spatial_weights, gdf) for i in range(spatial_weights_length)]
        weight_info = p.map(process_weights, args)
    #将嵌套列表展开
    weight_info = [item for sublist in weight_info for item in sublist]

    return weight_info


def __read_features_to_dataframe(shp_path, z_field, id_field,z_scale_factor=1):
    """
    从shapefile中获取X,Y,Z,ID(私有方法，不应在外部调用)   

    参数：
    -------------
    shp_path: shp文件路径
    z_field:  z字段名
    id_field: id字段名   
    z_scale_factor: 用于将z值乘以一个系数，以适应空间权重矩阵的计算

    返回值：
    -------------
    DataFrame，其中columns=['X', 'Y', 'Z', 'ID']

    """
    data = []
    with arcpy.da.SearchCursor(shp_path, ["SHAPE@XY", z_field, id_field]) as cursor:
        for row in cursor:
            x, y = row[0]
            z = row[1]*z_scale_factor
            id_value = row[2]
            data.append((x, y, z, id_value))
    df = pd.DataFrame(data, columns=['X', 'Y', 'Z', 'ID'])
    return df


def __read_features_to_dataframe_noele(shp_path,id_field):
    """
    从shapefile中获取X,Y,ID(私有方法，不应在外部调用)
    不读取Z值   

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


def cal_weight_txt(shp_path,out_path,z_field,id_field,distance_function,
                threshold=float('inf'),elevation=False,software='arcgis',thread_num=cpu_count(),
                z_scale_factor=1):
    """
    主函数，计算Arcgis权重矩阵文件请调用此函数
    计算全局Moran'I指数，调用该函数后会计算全局Moran'I指数，
    并在指定路径生成txt和png分析结果

    参数：
    -----------
    shp_path:          shapefile文件位置
    out_path:          输出的权重文件路径，不需要加扩展名
    z_field:           高程字段
    id_field:          唯一ID字段
    distance_function: 空间关系概念化函数，可选threshold/gaussian/inverse
    threshold:         距离阈值，留空则为不设置阈值
    elevation:         是否在距离计算中考虑高程影响
    software:          适用于空间分析软件的格式，可选arcgis/geoda，arcgis:.txt/geoda:.kwt
    thread_num:        使用的线程数，不应设置为大于CPU线程的值，默认使用全部CPU线程
    z_scale_factor:    Z值缩放系数

    返回值：
    -------------
    无
    """
    # arcpy.AddMessage("Reading Shapefile...")
    print("Reading Shapefile...")
    #获取dataframe
    if elevation==True:
        gdf = __read_features_to_dataframe(shp_path,z_field,id_field,z_scale_factor)
    else:
        gdf=__read_features_to_dataframe_noele(shp_path,id_field)

    shp_name=os.path.basename(shp_path).split('.')[0]

    #设置当前工作目录
    arcpy.env.workspace=os.getcwd()

    # arcpy.AddMessage("Calculating Weights Matrix...")
    print("Calculating Weights Matrix...")
    if distance_function=='threshold':
        spatial_weights_list = __threshold_weights(gdf,threshold,elevation,thread_num)
    elif distance_function=='gaussian':
        spatial_weights_list = __gauss_weights(gdf,threshold,elevation,thread_num)  
    elif distance_function=='inverse':
        spatial_weights_list = __inverse_weights(gdf,threshold,elevation,thread_num)
    else:
        raise ValueError("distance_function must be 'threshold', 'gaussian' or 'inverse'")

    print("Generate Weight File...")

    # 将权重信息写入到txt文件
    if software=='arcgis':
        out_path=out_path+'.txt'
        with open(out_path, 'w',encoding="ascii") as f:
            f.write(id_field+"\n")
            for info in spatial_weights_list:
                f.write(f"{info}\n")

    elif software=='geoda':
        out_path=out_path+'.gwt'
        with open(out_path, 'w',encoding="gb2312") as f:
            f.write("0 "+str(len(gdf))+" "+shp_name+" "+id_field+"\n")
            for info in spatial_weights_list:
                f.write(f"{info}\n")
    del spatial_weights_list
    gc.collect()

    print("Done!")



def cal_weight_txt_folder(shp_folder,output_folder,z_field,id_field,distance_function,
                threshold=float('inf'),elevation=False,software="arcgis",thread_num=cpu_count(),
                z_scale_factor=1):
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
    software:          适用于空间分析软件的格式，可选arcgis/geoda
    z_scale_factor:    Z值缩放系数

    返回值：
    -------------
    无
    """
    for file_name in os.listdir(shp_folder):
        if file_name.endswith(".shp"):
            shp_path = os.path.join(shp_folder, file_name)
            output_file = os.path.join(output_folder, file_name.split('.')[0])
            print(f"Processing {file_name}...")
            cal_weight_txt(shp_path,output_file,z_field,id_field,distance_function,
                threshold,elevation,software,thread_num,z_scale_factor)




if __name__ == "__main__":
    #在此处设置参数

    shp_path="F:\大创数据\中间产出的数据\上海市和重庆市处理好的住宅区点\处理好的重庆市主城区住宅区_80抽稀.shp"
    out_path="D:\Lenovo\Desktop\云南大学\大创\程序代码\空间权重矩阵测试\空间权重矩阵中间文件\不考虑高程\\geoda\重庆市住宅区_反距离权_考虑高程_16590"

    # shp_folder='F:\大创数据\中间产出的数据\对云南省和黄淮海平原创建的样方\云南省'
    # output_folder='D:\Lenovo\Desktop\云南大学\大创\程序代码\空间权重矩阵测试\空间权重矩阵中间文件\考虑高程\\geoda'

    z_field='Z'
    id_field='ID'
    distance_function='inverse'
    # threshold=10640  #上海
    threshold=16590  #重庆
    elevation=False
    software='geoda'
    thread_num=16
    z_factor=2

    cal_weight_txt(shp_path,out_path,z_field,id_field,distance_function,
                   threshold,elevation,software,thread_num)




