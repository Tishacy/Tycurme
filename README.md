<p style="text-align:center; font-family: Georgia; font-style: italic; font-size:45px; font-weight: 700; margin-bottom: 0px">Tycurme</p>
<p style="text-align:center; color:#999; font-family: Georgia; font-style: italic; margin-top: 0px;">Theis 配线法自动拟合程序</p>
<ul style="text-align:left; color:#999; font-family: Georgia; margin:0 0 0 350px; text-decoration:none; padding:0;">
    <li>author: TIANShuo</li>    
    <li>date: 2018/6/12</li>
</ul>

## 操作环境
- python 3.0及其以上版本
- 依赖库：
	```bash
	pandas>=0.24.2
	matplotlib>=3.0.1
	numpy>=1.13.3
	```
	安装依赖：
	```bash
	$ pip3 install -r requirements.txt
	```

## 操作方法

1. 解压文件，并进入解压得到的文件夹

2. 首先在`s-t.xlsx`文件中输入实测数据，并将该文件保存
  - 第一列为时间，以 $\text{min}$ 为单位
  - 第二列为流量，以 $\text{m}^3/\text{h}​$ 为单位
  - 第三列及后面所有列数据为观测孔数据
  	- 第一行数据为观测孔距离抽水井距离 $r$，以m为单位
  	- 其余行数据为降深值 $s$，以m为单位
  > 注意：观测孔数目不限，可以往后任意添加数据

3. 在文件夹中运行`tycurme.py`文件：
  ```bash
  $ python3 tycurme.py
  ```
  界面分为三个区域：图像区、控制区、结果区。

  <img src="./imgs/GUI.png" style="width: 500px;">

  1. 图像区：红色散点为实测数据，黑色曲线为Theis配线标准曲线，采用双对数坐标
  2. 控制区(Control Zone)：
  	- `x_bias`为横坐标偏移量，`y_bias`为纵坐标偏移量
  	- 手动拟合模式，有三种方法可以改变偏移量：
  		- 滑动条控制
  		- 输入框输入
  		- 加减按钮（输入框两侧）
  	- 自动拟合模式：共有三种自动拟合算法
  		- Traverse算法：遍历横、纵偏移量，并取RMSE（均方根误差）最小，得到最佳拟合情况）
  		- Mass Center算法：通过计算实验数据的质心坐标，与等范围段的标准曲线的质心坐标作对比，并取RMSE值最小。
  		- Slope算法：通过计算实验数据的线性拟合方程，与等范围段的标准曲线的线性拟合方程作对比，并取RMSE值最小。
  3. 结果区(Results Zone)：
  	分别显示$W$、$1/u$、$s$、$t/r^2$、$T$、$S$值，并显示所拟合的RMSE值

  > 注意：图像区、控制区、结果区相互关联，即控制区的偏移量改变时，图像与结果会实时变化。
