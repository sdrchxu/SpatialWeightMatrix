import geopandas as gpd
from pysal.model import spreg
from mgwr.gwr import GWR
from mgwr.sel_bw import Sel_BW
from sklearn.neighbors import KernelDensity
from sklearn.linear_model import LinearRegression
import numpy as np
import pandas as pd

def calculate_distance_with_elevation(point1, point2):
    # 计算delta X,Y,Z
    x_diff = point1['X'] - point2['X']
    y_diff = point1['Y'] - point2['Y']
    z_diff = point1['Z'] - point2['Z']

    distance = np.sqrt(x_diff**2 + y_diff**2 + z_diff**2)
    return distance

def calculate_aicc(kde, data):
    # 计算AICc值
    n = len(data)
    log_likelihood = np.sum(kde.score_samples(data.reshape(-1, 1)))
    aicc = -2 * log_likelihood + 2 * kde.n_features_in_ + 2 * (kde.n_features_in_ + 1) / (n - kde.n_features_in_ - 1)
    return aicc

def spatial_weighted_regression(shapefile_path, location_columns, dependent_variable_column):
    # 读取shapefile
    gdf = gpd.read_file(shapefile_path)

    # 计算每个点对其他点的带高程距离
    distances = []
    for i, point1 in gdf.iterrows():
        for j, point2 in gdf.iterrows():
            if i != j:
                distance = calculate_distance_with_elevation(point1, point2)
                distances.append(distance)

    # 转换为NumPy数组
    distances = np.array(distances).reshape(-1, 1)

    # 使用KernelDensity估计概率密度
    kde = KernelDensity(bandwidth=0.1, kernel='gaussian')
    kde.fit(distances)

    # 计算AICc值
    aicc = calculate_aicc(kde, distances)

    print(f"Bandwidth selected based on AICc: {kde.bandwidth}")
    print(f"AICc value: {aicc}")

    # 进行地理加权回归分析
    X_location = gdf[location_columns].values
    y = gdf[dependent_variable_column].values

    kde_bandwidth = kde.bandwidth
    kde_weights = np.exp(kde.score_samples(distances))

    weights_matrix = np.diag(kde_weights)

    X_weighted = np.dot(weights_matrix, X_location)
    y_weighted = kde_weights * y

    # 使用线性回归进行地理加权回归
    model = LinearRegression()
    model.fit(X_weighted, y_weighted)

    coefficients = model.coef_

    print("Spatial Weighted Regression Coefficients:")
    for i, column in enumerate(location_columns):
        print(f"{column}: {coefficients[i]}")
    print(f"Dependent Variable: {dependent_variable_column}")

    # 将结果保存到txt文件
    result_df = pd.DataFrame({'Coefficient': coefficients, 'AICc': aicc},
                             index=location_columns)
    result_df.to_csv('spatial_weighted_regression_result.txt', sep='\t', index_label='Variable')

# 调用函数进行地理加权回归分析
shapefile_path = 'path/to/your/shapefile.shp'
location_columns = ['X', 'Y', 'Z']
dependent_variable_column = 'YourDependentVariable'

spatial_weighted_regression(shapefile_path, location_columns, dependent_variable_column)















# 读取 shapefile 文件
shp_path = "D:\Lenovo\Desktop\云南大学\大创\复杂地形条件对空间关系构建的影响研究\data\云南省火点数据回归\FinshinNet_Point_Shp.shp"
gdf = gpd.read_file(shp_path)

# 提取自变量、因变量和高程数据
y = gdf['dependent_variable']
X = gdf[['independent_variable1', 'independent_variable2']]
elevations = gdf['elevation']  # 高程字段
gdf['geometry'] = gdf['geometry'].centroid.apply(lambda point: (point.x, point.y, point.z) if point.has_z else (point.x, point.y))

# 创建 GWR 模型
bw = Sel_BW(gdf[['geometry']], y, gdf.geometry, spherical=True)
bw = bw.search(criterion='AICc')

gwr_model = GWR(y, X, bw=bw, fixed=False, kernel='gaussian')
gwr_results = gwr_model.fit()

# 输出结果到文本文件
output_path = 'gwr_results.txt'
with open(output_path, 'w') as output_file:
    output_file.write(gwr_results.summary())




