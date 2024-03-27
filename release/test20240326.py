
import cal_moran_pysal


folder_path="F:\大创数据\中间产出的数据\云南省和黄淮海平原已处理好的火点\云南省逐月火点"
field='FRP'
output_folder="D:\Lenovo\Desktop\云南大学\大创\程序代码\空间权重矩阵测试\阈值权函数结果\考虑高程2_云南省逐月_阈值权_FRP_55KM"
distance_function="threshold"
threshold=55000
std=True
elevation=True
z_scale_factor=6

cal_moran_pysal.global_moran_folder(folder_path,field,output_folder,distance_function,threshold,std,elevation,z_scale_factor)