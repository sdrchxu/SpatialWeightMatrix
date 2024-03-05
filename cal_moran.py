import geopandas as gpd
import pandas as pd
import numpy as np
from pysal.lib import weights
from scipy.sparse import csr_matrix
import scipy
import shapefile
from scipy.spatial.distance import cdist

def inverse_weights(gdf,L0):
    """
    计算反距离权矩阵

    参数：
    ------------
    gdf:         shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含空间权重矩阵的DataFrame
    """
    # 提取所需的字段（ID, X, Y, Z）
    df = gdf[['ID', 'X', 'Y', 'Z']].copy()

    # 设置距离阈值L0

    # 创建一个空的空间权重矩阵
    spatial_weights = np.zeros((len(df), len(df)))

    # 计算空间权重矩阵
    for i in range(len(df)):
        for j in range(len(df)):
            # 计算两要素间的实际距离
            distance = np.sqrt((df['X'][i] - df['X'][j])**2 + (df['Y'][i] - df['Y'][j])**2 + (df['Z'][i] - df['Z'][j])**2)
            
            # 判断距离是否小于阈值
            if distance <= L0:
                spatial_weights[i, j] = 1

    # 将空间权重矩阵保存到DataFrame中
    spatial_weights_df = pd.DataFrame(spatial_weights, columns=df['ID'], index=df['ID'])
    return spatial_weights_df


def gauss_weights(gdf,L0):
    """
    计算高斯权矩阵

    参数：
    ------------
    gdf:    shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含空间权重矩阵的DataFrame
    """

    # 获取字段X、Y、Z和ID的索引
    fields = [field[0] for field in gdf.fields[1:]]
    x_index = fields.index('X')
    y_index = fields.index('Y')
    z_index = fields.index('Z')
    id_index = fields.index('ID')

    # 获取点坐标和ID
    points = []
    ids = []
    for sr in gdf.shapeRecords():
        points.append((sr.record[x_index], sr.record[y_index], sr.record[z_index]))
        ids.append(sr.record[id_index])

    # 转换为numpy数组
    data = np.array(points)

    # 设置距离阈值
    threshold = L0  # 单位：米

    # 计算空间权重矩阵
    num_points = len(ids)
    weight_matrix = np.zeros((num_points, num_points))

    for i in range(num_points):
        center_point = data[i]
        distances = cdist([center_point], data)[0]
        weights = np.exp(-distances**2 / (2 * threshold**2))
        weight_matrix[i] = weights

    # 创建DataFrame并设置索引和列名
    df_weights = pd.DataFrame(weight_matrix, index=ids, columns=ids)
    return df_weights

def threshold_weights(gdf,L0):
    """
    计算阈值权矩阵

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

    # 创建一个空的空间权重矩阵
    spatial_weights = np.zeros((len(df), len(df)))

    # 计算空间权重矩阵
    for i in range(len(df)):
        for j in range(len(df)):
            # 计算两要素间的实际距离
            distance = np.sqrt((df['X'][i] - df['X'][j])**2 + (df['Y'][i] - df['Y'][j])**2 + (df['Z'][i] - df['Z'][j])**2)
            
            # 判断距离是否小于阈值
            if distance <= L0:
                spatial_weights[i, j] = 1

    # 将空间权重矩阵保存到DataFrame中
    spatial_weights_df = pd.DataFrame(spatial_weights)
    return spatial_weights_df


def global_moran(shp_path,field,distance_function,threshold=float('inf'),std=True):
    """
    计算全局Moran'I指数

    参数：
    -----------
    shp_path:          shapefile文件位置
    field:             输入字段
    distance_function: 空间关系概念化函数，可选threshold/gauss/inverse
    threshold:         距离阈值，留空则为不设置阈值
    std:               标准化方法，True为行标准化，False为不进行标准化
    """
    gdf = gpd.read_file(shp_path)
    if distance_function=='threshold':
        spatial_weights_df = threshold_weights(gdf,threshold)
    elif distance_function=='gauss':
        spatial_weights_df = gauss_weights(gdf,threshold)
    elif distance_function=='inverse':
        spatial_weights_df = inverse_weights(gdf,threshold)
    sparse_matrix=scipy.sparse.csr_matrix(spatial_weights_df.values)

    # 创建空间权重对象,WSP 类是pysal中的一个子类，表示"weights spatial",专用于表示处理空间权重矩阵的子类
    wsp = weights.WSP(sparse_matrix)
    #WSP无法直接参与Moran'I指数计算，所以需要将它转为w类
    w=weights.W.from_WSP(wsp)
    w.transform = 'r' #对空间权重矩阵进行行标准化
    print(type(w))

    # 如果权重矩阵是对称的，可以使用Sym的子类
    # w = weights.WSP(symmetrize='True', values=spatial_weights_df.values)

    # 打印空间权重矩阵
    print(w.sparse)