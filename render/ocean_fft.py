"""远洋水面 FFT 级联相关的纯数学辅助函数。

模块只做与渲染后端无关的数值计算，便于单元测试。约定如下：

* 世界坐标 y 轴向上，水面位于 ``water_y``；
* ``kx``/``ky`` 为水平波数（单位 rad/m），风向角 ``wind_angle`` 以 +x 轴
  为 0、逆时针为正；
* 各级 FFT 级联共享同一套风速/风向插值结果，但拥有各自的波长范围。
"""

from __future__ import annotations

import math

# 重力加速度 (m/s^2)，Phillips 谱的特征长度 L = V^2 / g。
GRAVITY = 9.81

# 每个级联相对上一级的风速倍率：越高的级联覆盖越大的波长，
# 对应的局地风速也越强，避免切换级联时浪高突变。
DEFAULT_CASCADE_SCALE = 1.25

# Phillips 谱的整体幅度。测试只约束符号与单调性，取 1.0 便于阅读；
# 实际渲染时可由上层按海况统一缩放。
PHILLIPS_AMPLITUDE = 1.0


def cascade_wind(base_speed: float, index: int, scale: float = DEFAULT_CASCADE_SCALE) -> float:
    """返回第 ``index`` 级 FFT 级联的风速。

    所有级联的风速都由同一个基础风速按几何级数插值得到，这样从三级浪切到
    风暴浪时各级之间不会出现风速跳变：

    ``speed(i) = base_speed * scale ** i``

    ``index == 0`` 时即为基础风速。
    """
    index = int(index)
    if index < 0:
        raise ValueError("cascade index must be non-negative")
    if scale <= 0.0:
        raise ValueError("cascade scale must be positive")
    return float(base_speed) * float(scale) ** index


def phillips(kx: float, ky: float, wind_angle: float, wind_speed: float) -> float:
    """返回 Phillips 谱能量。

    使用经典的 Phillips 形式：

    ``P(k) = A * exp(-1 / (kL)^2) / k^4 * (k_hat · w_hat)^2``

    其中 ``L = wind_speed^2 / g``。风向项保证顺风方向的浪占优，侧向/逆风
    方向能量趋近于 0；风速越大特征长度越长，同一波数上的能量也越高。
    直流分量（k = 0）与无风情况返回 0，避免除零。
    """
    kx = float(kx)
    ky = float(ky)
    wind_speed = float(wind_speed)

    k_squared = kx * kx + ky * ky
    if k_squared <= 0.0 or wind_speed <= 0.0:
        return 0.0

    k = math.sqrt(k_squared)
    wind_x = math.cos(float(wind_angle))
    wind_y = math.sin(float(wind_angle))

    # 归一化波矢与归一化风向的点积；始终取平方，正反方向的浪脊都计入。
    directional = (kx * wind_x + ky * wind_y) / k
    if directional == 0.0:
        return 0.0

    characteristic_length = wind_speed * wind_speed / GRAVITY
    exponent = -1.0 / (k_squared * characteristic_length * characteristic_length)
    # 极端低波数下 exp(-1/(kL)^2) 会下溢为 0，这属于物理上的长波截止。
    spectrum = math.exp(exponent) / (k_squared * k_squared)
    return PHILLIPS_AMPLITUDE * spectrum * directional * directional


def fold_normal(dx: float, dz: float, spacing: float) -> tuple[float, float, float]:
    """返回位移贴图折叠后重算的单位法线。

    ``dx``/``dz`` 为水平位移在采样间距上的增量，``spacing`` 为网格间距。
    先把位移换算成坡度，再由高度场法线公式重建：

    ``n = normalize((dx / spacing, 1, dz / spacing))``

    位移贴图折叠（Jacobian 变号）后必须用它重算法线，不能沿用折叠前的
    法线，否则浪尖处着色与反射方向会穿帮。
    """
    dx = float(dx)
    dz = float(dz)
    spacing = abs(float(spacing))
    if spacing == 0.0:
        # 退化网格没有可用的坡度尺度，退化为朝上的单位法线。
        return (0.0, 1.0, 0.0)

    slope_x = dx / spacing
    slope_z = dz / spacing
    length = math.sqrt(slope_x * slope_x + 1.0 + slope_z * slope_z)
    return (slope_x / length, 1.0 / length, slope_z / length)


def foam_threshold(jacobian: float, resolution: int) -> float:
    """返回与分辨率无关的泡沫阈值。

    泡沫由位移映射的 Jacobian 决定：Jacobian 越小（越接近折叠/压缩），
    泡沫覆盖越高。这里把 Jacobian 归一化为 0..1 的阈值，**不引入任何
    ``resolution`` 因子**，因此同一海况在 64 与 1024 网格上泡沫分布一致，
    不会随贴图分辨率闪烁。``resolution`` 仅保留在签名中以兼容调用方。
    """
    jacobian = float(jacobian)
    threshold = 1.0 - jacobian
    if threshold < 0.0:
        return 0.0
    if threshold > 1.0:
        return 1.0
    return threshold


def clipmap_uv(world_uv: float, ring: int, cell: float) -> float:
    """返回 clipmap 环内对齐后的 UV。

    第 ``ring`` 环的单元尺度为 ``cell * 2 ** ring``。把世界 UV 对环尺度取模
    即可得到该环的局部 UV，使相邻环在边界处共享同一套网格对齐规则，
    消除环与环之间的 UV 接缝。内层环的结果落在 ``[0, cell)``，外层环落在
    ``[0, cell * 2**ring)``，因此第 ``ring`` 环的 UV 是第 ``ring + 1`` 环 UV
    的嵌套细化；负坐标按 Python ``%`` 语义环绕到同一区间。
    """
    cell = float(cell)
    ring = int(ring)
    if ring < 0:
        raise ValueError("clipmap ring must be non-negative")
    if cell <= 0.0:
        return float(world_uv)

    ring_extent = cell * (2 ** ring)
    return float(world_uv) % ring_extent


def buoyancy_height(wave_height: float, visual_offset: float) -> float:
    """返回浮力采样使用的物理波高。

    浮力必须由物理波高（FFT 位移场的 y 分量）驱动，而不是视觉网格的
    顶点高度；``visual_offset`` 是渲染用的偏移（例如为了避免 Z-fighting
    或贴合视觉网格），不能进入物理计算，否则船体浮沉会整体偏置。
    """
    return float(wave_height)


def is_underwater(camera_y: float, water_y: float) -> bool:
    """判断相机是否在水下。

    ``water_y`` 以“水面相对基准面的下沉量”给出（正值向下），因此水面的
    世界高度为 ``-water_y``；相机高度低于该水面即视为水下：

    ``camera_y < -water_y``

    相机恰好位于水面时不判为水下，避免水下折射/雾分支在切换帧上抖动。
    该判定的结果用于决定先做折射再叠加雾，还是先雾后折射。
    """
    return float(camera_y) < -float(water_y)
