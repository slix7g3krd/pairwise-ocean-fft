# pairwise-ocean-fft

远洋水面从三级浪切到风暴浪时，浪尖穿船、泡沫闪烁、远处接缝：FFT 级联风速和风向各自插值；Phillips 谱没按级联缩放；位移贴图折叠后法线未重算；Jacobian 泡沫阈值随分辨率变；clipmap 环间 UV 不连续；平面反射只反射低阶浪；浮力采样用视觉网格而不是物理波高；GPU 回读晚一帧；浪花粒子与泡沫图重复；相机在水下时折射和雾顺序错。render/ocean_fft.py 是核心，tests/test_ocean_fft.py 15个红。要修的是：级联参数统一插值；谱按级联缩放；折叠后重算法线；泡沫阈值分辨率无关；clipmap UV 对齐；反射包含全级联；浮力走物理采样；回读同步；泡沫去重；水下折射排序。

## 运行

```bash
python -m pytest tests -q
```

初始状态：15 个用例失败。修复实现使全部用例通过，不要修改 `tests/` 下的测试。

## 实现约定

`render/ocean_fft.py` 中各辅助函数的契约（与测试一一对应）：

| 函数 | 约定 |
| --- | --- |
| `cascade_wind(base, i, scale=1.25)` | 各级联共用一套插值：`base * scale**i`，`i = 0` 即基础风速 |
| `phillips(kx, ky, angle, speed)` | 经典 Phillips 谱，特征长度 `L = V²/g`，风向项 `(k̂·ŵ)²`，`k = 0` 或无风时为 0 |
| `fold_normal(dx, dz, spacing)` | `normalize((dx/spacing, 1, dz/spacing))`，折叠后重算 |
| `foam_threshold(jacobian, resolution)` | `clamp(1 - jacobian, 0, 1)`，与 `resolution` 无关 |
| `clipmap_uv(world_uv, ring, cell)` | `world_uv % (cell * 2**ring)`，环间嵌套对齐 |
| `buoyancy_height(wave_height, visual_offset)` | 只用物理波高，忽略视觉偏移 |
| `is_underwater(camera_y, water_y)` | `water_y` 为水面相对基准面的下沉量，判据 `camera_y < -water_y` |

子方向：图形渲染；任务类型：Bugfix；语言/框架：Python 3, pytest；环境可复现等级：无外部依赖。
