import geopandas as gpd
import pandas as pd
import numpy as np
from pysal.lib import weights
from scipy.sparse import csr_matrix
import scipy
from scipy.spatial.distance import cdist
from pysal.explore import esda
from splot.esda import plot_moran
import matplotlib.pyplot as plt
import os
import math

def __inverse_weights(gdf,L0,elevation):
    """
    计算反距离权矩阵(私有方法，不应在外部调用)

    参数：
    ------------
    gdf:         shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含空间权重矩阵的DataFrame
    """
    # 提取所需的字段（ID, X, Y, Z）
    # df = gdf[['ID', 'X', 'Y', 'Z']].copy()


    # 计算距离矩阵
    if elevation==True:
        coordinates = gdf[['X', 'Y', 'Z']].values
    elif elevation==False:
        coordinates = gdf[['X', 'Y']].values
    distances = cdist(coordinates, coordinates)
    # 应用反距离权重
    spatial_weights = np.where((distances!=0)&(distances <= L0), 1/distances, 0)
    # 将空间权重矩阵保存到DataFrame中
    spatial_weights_df = pd.DataFrame(spatial_weights)
    return spatial_weights_df




def __gauss_weights(gdf,L0,elevation):
    """
    计算高斯权矩阵(私有方法，不应在外部调用)

    参数：
    ------------
    gdf:    shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含空间权重矩阵的DataFrame
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

    # 创建空间权重矩阵的DataFrame
    spatial_weights_df = pd.DataFrame(weights2, index=gdf['ID'], columns=gdf['ID'])

    return spatial_weights_df



def __threshold_weights(gdf,L0,elevation):
    """
    计算阈值权矩阵(私有方法，不应在外部调用)

    参数：
    ------------
    gdf:    shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含空间权重矩阵的DataFrame
    """
    # 提取所需的字段（ID, X, Y, Z）
    df = gdf[['ID', 'X', 'Y', 'Z']].copy()
    if elevation==True:
        coordinates = gdf[['X', 'Y', 'Z']].values
    elif elevation==False:
        coordinates = gdf[['X', 'Y']].values

    #计算要素间的距离
    # coordinates = df[['X', 'Y']].values
    distances = cdist(coordinates, coordinates)

    # 对得到的距离进行判断并赋值
    spatial_weights = np.where((distances<=L0)&(distances!=0), 1, 0)

    # 将空间权重矩阵保存到DataFrame中
    spatial_weights_df = pd.DataFrame(spatial_weights)
    return spatial_weights_df


def global_moran(shp_path,field,output_file,distance_function,threshold=float('inf'),std=True,elevation=True):
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
    print("Reading Shapefile...")
    gdf = gpd.read_file(shp_path)
    print("Calculating Weights Matrix...")
    if distance_function=='threshold':
        spatial_weights_df = __threshold_weights(gdf,threshold,elevation)
    elif distance_function=='gaussian':
        spatial_weights_df = __gauss_weights(gdf,threshold,elevation)  #高斯权函数需要传入一个Shapefile对象
    elif distance_function=='inverse':
        spatial_weights_df = __inverse_weights(gdf,threshold,elevation)
    else:
        raise ValueError("distance_function must be 'threshold', 'gaussian' or 'inverse'")
    sparse_matrix=scipy.sparse.csr_matrix(spatial_weights_df.values)

    print("Calculating Moran'I...")
    # 创建空间权重对象,WSP 类是pysal中的一个子类，表示"weights spatial",专用于表示处理空间权重矩阵的子类
    wsp = weights.WSP(sparse_matrix)
    #WSP无法直接参与Moran'I指数计算，所以需要将它转为w类
    w=weights.W.from_WSP(wsp)
    if std==True:
        w.transform = 'r' #对空间权重矩阵进行行标准化

    # 如果权重矩阵是对称的，可以使用Sym的子类
    # w = weights.WSP(symmetrize='True', values=spatial_weights_df.values)

    moran = esda.Moran(gdf[field], w)

    print("Saving Result Reports...")
    output_txt = output_file + '.txt'
    output_png = output_file + '.png'

    with open(output_txt, 'w') as f:
        f.write('Moran\'s I: {}\n'.format(moran.I))
        f.write('P-value: {}\n'.format(moran.p_sim))
        f.write('Z-score: {}\n'.format(moran.z_sim))

    plot_moran(moran, zstandard=True, figsize=(20, 8))
    plt.savefig(output_png)
    plt.close()
    print("Done!")
    return moran.I,moran.p_sim,moran.z_sim



def global_moran_folder(folder_path, field, output_folder, distance_function, threshold=float('inf'), std=True, elevation=True):
    """
    计算给定文件夹中所有shapefile文件的Moran'I指数

    参数：
    -----------
    folder_path:       shapefile文件夹位置
    field:             输入字段
    output_folder:     结果输出文件夹路径
    distance_function: 空间关系概念化函数，可选threshold/gaussian/inverse
    threshold:         距离阈值，留空则为不设置阈值
    std:               标准化方法，True为行标准化，False为不进行标准化
    elevation:         是否在距离计算中考虑高程影响

    返回值：
    -------------
    None
    """
    results=[]
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".shp"):
            shp_path = os.path.join(folder_path, file_name)
            output_file = os.path.join(output_folder, os.path.splitext(file_name)[0])
            print(f"Processing {file_name}...")
            I,p_sim,z_sim=global_moran(shp_path, field, output_file, distance_function, threshold, std, elevation)
            results.append({
                "File": file_name,
                "Moran\'s I": I,
                "P-value": p_sim,
                "Z-score": z_sim
            })
    output_csv = os.path.join(output_folder, "spatial_autocorrelation_results.csv")
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False,encoding='GB2312')


def output_gauss_weights(shp_path,L0,elevation,output_file_path):
    """
    输出高斯权函数矩阵

    参数：
    ------------
    shp_path:    shapefile文件路径
    L0:          阈值
    elevation:   是否在距离计算中考虑高程影响
    output_file_path: 输出csv文件的路径

    返回值：
    -------------
    空
    """
    gdf = gpd.read_file(shp_path)
    sigma=1.0 #标准方差，控制了函数的曲线在尖峰周围的陡峭程度
    gaussian_const=math.pow(math.pi*2,-0.5)#高斯常数
    # 转为numpy数组
    if elevation==True:
        coordinates = gdf[['X', 'Y', 'Z']].values
    elif elevation==False:
        coordinates = gdf[['X', 'Y']].values
    # 计算距离矩阵
    distances = cdist(coordinates, coordinates,'euclidean')

    # 计算高斯矩阵
    # weights = np.where(distances<=L0,np.where(distances==0,0,(gaussian_const*np.exp(-distances**2 / 2*sigma**2))),0)
    weights1 = np.where(distances<=L0,(np.sqrt(distances/L0)),0)
    weights2 = np.where(distances<=L0,(gaussian_const*np.exp(-np.power(weights1,2) / 2)),0)

    # 创建空间权重矩阵的DataFrame并输出csv文件
    spatial_weights_df = pd.DataFrame(weights2, index=gdf['ID'], columns=gdf['ID'])
    spatial_weights_df.to_csv(output_file_path)

    # 创建空间权重矩阵的DataFrame并输出txt文件
    # with open(output_file_path, 'w') as f:
    #     for i in range(len(weights2)):
    #         for j in range(len(weights2[i])):
    #             if weights2[i][j] != 0:
    #                 f.write(f"{gdf['ID'][i]} {gdf['ID'][j]} {weights2[i][j]}\n")
    print("done")


def output_inverse_weights(shp_path,L0,elevation,output_path):
    """
    输出反距离权矩阵

    参数：
    ------------
    shp_path:    shapefile文件路径
    L0:          阈值
    output_path: csv输出路径

    返回值：
    -------------
    空
    """
    # 提取所需的字段（ID, X, Y, Z）
    # df = gdf[['ID', 'X', 'Y', 'Z']].copy()

    gdf=gpd.read_file(shp_path)
    # 计算距离矩阵
    if elevation==True:
        coordinates = gdf[['X', 'Y', 'Z']].values
    elif elevation==False:
        coordinates = gdf[['X', 'Y']].values
    distances = cdist(coordinates, coordinates)
    # 应用反距离权重
    spatial_weights = np.where(distances==0, 0, np.where(distances <= L0,1/distances,0))
    # 将空间权重矩阵保存到DataFrame中
    spatial_weights_df = pd.DataFrame(spatial_weights)
    spatial_weights_df.to_csv(output_path)
