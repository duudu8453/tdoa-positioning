# 时差定位（TDOA）原理演示

这是一个用于演示时差定位（TDOA，Time Difference of Arrival）原理的交互式可视化项目。该项目提供了一个直观的界面，帮助用户理解时差定位算法的工作原理、误差影响以及几何构型对定位精度的影响。

## 功能特点

- 交互式接收站布置：用户可以自由添加、移动、删除接收站
- 实时定位计算：根据接收站位置和时间差数据实时计算目标位置
- 多种双曲线表示：支持参数方程、等高线、极坐标和距离差直接计算四种方式
- 误差模拟：包含时间测量误差、位置误差、大气延迟和多径效应
- 几何构型分析：计算并显示GDOP（几何精度因子）、覆盖面积等指标
- 动态信号传播演示：可视化信号从发射源到各接收站的传播过程

## 文件结构

- `time_difference_positioning.html`：主要的演示界面，包含完整的交互功能
- `test_tdoa_algorithm.js`：TDOA算法的核心实现，包括最小二乘求解等
- `algorithm_validation.html`：算法验证页面，用于测试定位精度

## 使用方法

1. 打开`time_difference_positioning.html`文件
2. 通过点击添加或拖动调整接收站位置
3. 观察双曲线的绘制和交点
4. 调整各种误差参数观察其影响
5. 点击"开始动画"查看信号传播过程

## 时差定位原理

时差定位（TDOA）是一种利用信号到达不同接收站的时间差来确定信号源位置的技术。其基本原理是：

1. 多个接收站接收同一信号
2. 计算信号到达不同接收站的时间差
3. 每对接收站的时间差确定一条双曲线
4. 多条双曲线的交点即为信号源位置

本项目使用迭代最小二乘法求解TDOA方程组，以获得最优的位置估计。