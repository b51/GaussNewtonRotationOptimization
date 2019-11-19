## GaussNewtonRotationOptimization

### What's this repository for
Optimizatin Rotation with known t_cw, several pair of world pose and camera projection pixel corrdination, get optimized Rotation Matrix R_cw
已知平移矩阵 t_cw 和几对世界坐标到相平面的投影坐标，使用高斯牛顿方法求解旋转矩阵 R_cw
由于平移矩阵已知，所以该方法其实为 PNP 的简单模式

### Reference
#### Rotation Jacobian calculation
https://www.cnblogs.com/gaoxiang12/p/5689927.html
or 视觉 SLAM 十四讲 P.75, P.146

#### GaussNewton Method

高斯牛顿法用来优化最小二乘法

基本原理如下：

$$\mathbf{error = y - f(x)}$$
其中y为真值, f(x) 为计算值

$$\mathbf{F(x)} = \mathbf{error}^T * \mathbf{error}$$

$$\mathbf{error(x+\delta{x}) = error(x) + J\delta{x}}$$

$$\mathbf{F(x+\delta{x}) = (error(x) + J\delta{x})^T * (error(x) + J\delta{x})
= error(x)^Terror(x) + 2error(x)^TJ\delta{x} + \delta{x}^TJ^TJ\delta{x}}
$$
根据泰勒定理, 有
$$\mathbf{F'(x) = 2error(x)^TJ + 2\delta{x}^TJ^TJ}$$

当F(x)的导数为 0 时，F(x) 有极小值, 有
$$\mathbf{J^TJ\delta{x} = -J^Terror}$$

解该方程得到$$\mathbf{\delta{x}}$$ 并使用该值对 x 进行更新 $$\mathbf{x = x + \delta{x}}$$ 直到收敛

### Usage
```bash
$ mkdir build
$ cd build
$ make -j
$ ./GNRotationOptimization
```
