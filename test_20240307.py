import cal_moran_2

folder_path="F:\大创数据\中间产出的数据\云南省和黄淮海平原已处理好的火点\云南省逐月火点"
field='FRP'
distance_function='gauss'
output_folder="D:\Lenovo\Desktop\云南大学\大创\程序代码\空间权重矩阵测试\高斯权函数结果\考虑高程2"
cal_moran_2.global_moran_folder(folder_path,field,output_folder,distance_function,threshold=55000,elevation=True)