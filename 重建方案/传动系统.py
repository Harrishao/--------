# -*- coding: utf-8 -*-
"""传动系统设计: 电机(2.2kW, 940r/min) → V带 → 轴I → 两级直齿轮 → 主轴(25r/min)
   布置: 三轴竖排在主轴正下方, 齿轮箱位于凸轮组对侧(Z<0), 带轮在轴I外端, 电机水平置右。"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Circle, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
fp = "/sessions/modest-dazzling-cori/.local/lib/python3.10/site-packages/mplfonts/fonts/NotoSansCJKsc-Regular.otf"
fm.fontManager.addfont(fp); plt.rcParams["font.family"] = fm.FontProperties(fname=fp).get_name()
plt.rcParams["axes.unicode_minus"] = False

# ---- 输入 ----
P, n_m, n_target = 2.2, 940.0, 25.0
i_total = n_m / n_target
# ---- 分配: 带 i1 ~2.36, 齿轮 z1/z2=17/77, z3/z4=17/60 ----
dd1, dd2 = 100.0, 236.0
i1 = dd2 / dd1
z1, z2, z3, z4 = 17, 77, 17, 60
m1, m2 = 2.0, 3.0
i2, i3 = z2/z1, z4/z3
i_act = i1 * i2 * i3
n_I, n_II, n_main = n_m/i1, n_m/(i1*i2), n_m/i_act
a1g, a2g = m1*(z1+z2)/2, m2*(z3+z4)/2
# ---- V带 (A型) ----
KA = 1.2; Pca = KA*P
v = np.pi*dd1*n_m/60000.0
a0 = 400.0
L0 = 2*a0 + np.pi*(dd1+dd2)/2 + (dd2-dd1)**2/(4*a0)
Ld = 1400.0
a_belt = a0 + (Ld-L0)/2
alpha1 = 180 - (dd2-dd1)/a_belt*57.3
P0, dP0, Ka_, KL = 0.95, 0.11, 0.95, 0.96
zb = Pca/((P0+dP0)*Ka_*KL)
zbn = int(np.ceil(zb))
# ---- 转矩 ----
eta_b, eta_g, eta_brg = 0.96, 0.97, 0.99
T_I  = 9550*P*eta_b/n_I
T_II = 9550*P*eta_b*eta_g*eta_brg/n_II
T_O  = 9550*P*eta_b*eta_g**2*eta_brg**2/n_main
print(f"i_total需求={i_total:.2f}  实际={i_act:.3f}  主轴转速={n_main:.2f} r/min (偏差{(n_main-25)/25*100:+.1f}%)")
print(f"带: A型 dd1={dd1:.0f} dd2={dd2:.0f} i1={i1:.3f} v={v:.2f}m/s L0={L0:.0f}→Ld={Ld:.0f} a={a_belt:.0f} 包角={alpha1:.1f}° 根数={zb:.2f}→{zbn}")
print(f"齿轮1: m={m1} z={z1}/{z2} i={i2:.3f} a={a1g:.1f} d={m1*z1:.0f}/{m1*z2:.0f}")
print(f"齿轮2: m={m2} z={z3}/{z4} i={i3:.3f} a={a2g:.1f} d={m2*z3:.0f}/{m2*z4:.0f}")
print(f"转速: 轴I={n_I:.1f} 轴II={n_II:.1f} 主轴={n_main:.2f} r/min")
print(f"转矩: 轴I={T_I:.1f} 轴II={T_II:.1f} 主轴={T_O:.1f} N·m")
# ---- 布置 ----
S_main = np.array([0.0, 0.0])
S_II   = np.array([0.0, -a2g])
S_I    = np.array([0.0, -a2g - a1g])
S_mot  = S_I + np.array([a_belt, 0.0])
print(f"轴II=(0,{-a2g:.1f})  轴I=(0,{-a2g-a1g:.1f})  电机轴=({a_belt:.0f},{-a2g-a1g:.1f})")

# ---- 图: 前视(XY) + 轴向站位(ZY) ----
fig, (ax, az) = plt.subplots(1, 2, figsize=(13.5, 8), gridspec_kw=dict(width_ratios=[1.5, 1]))
# 凸轮组包络
for rr, lab in [(82, "三凸轮最大包络 r≈82(粉筛v4)")]:
    th = np.linspace(0, 2*np.pi, 100)
    ax.fill(rr*np.cos(th), rr*np.sin(th), color="#dde3e8", alpha=.7, zorder=1)
ax.annotate("凸轮组包络 r≈82\n(Z=0~142.5)", (-78, 60), fontsize=8.5, color="#56606a")
ax.plot([0,0],[150,216],color="#888",lw=3,zorder=2); ax.annotate("下冲头", (6,160), fontsize=8, color="#666")
# 齿轮节圆
for c, r, lab, fc in [(S_main, m2*z4/2, f"z4={z4}, d={m2*z4:.0f}", "#cdd9c8"),
                      (S_II,   m2*z3/2, f"z3={z3}, d={m2*z3:.0f}", "#e7ddc8"),
                      (S_II,   m1*z2/2, f"z2={z2}, d={m1*z2:.0f}", "#cdd9c8"),
                      (S_I,    m1*z1/2, f"z1={z1}, d={m1*z1:.0f}", "#e7ddc8")]:
    ax.add_patch(Circle(c, r, fc=fc, ec="#5a6a50", lw=1.3, alpha=.85, zorder=3))
ax.annotate(f"z4 d={m2*z4:.0f} (主轴)", (95, -20), fontsize=9)
ax.annotate(f"z3 d={m2*z3:.0f}\nz2 d={m1*z2:.0f}", (82, -a2g-12), fontsize=9)
ax.annotate(f"z1 d={m1*z1:.0f}", (24, -a2g-a1g-6), fontsize=9)
# 带轮+带
ax.add_patch(Circle(S_I, dd2/2, fc="#d8cfe0", ec="#5a4a6a", lw=1.3, alpha=.85, zorder=2))
ax.add_patch(Circle(S_mot, dd1/2, fc="#d8cfe0", ec="#5a4a6a", lw=1.3, zorder=2))
for s in (+1, -1):
    ax.plot([S_I[0], S_mot[0]], [S_I[1]+s*dd2/2, S_mot[1]+s*dd1/2], color="#5a4a6a", lw=1.4, zorder=2)
ax.annotate(f"大带轮 dd2={dd2:.0f} (轴I)", (S_I[0]-205, S_I[1]-46), fontsize=9)
ax.annotate(f"小带轮 dd1={dd1:.0f}\n电机 Y112M-6\n2.2kW 940r/min", (S_mot[0]-46, S_mot[1]-110), fontsize=9)
ax.annotate(f"A型V带×{zbn}  a={a_belt:.0f}", ((S_I[0]+S_mot[0])/2-60, S_I[1]+78), fontsize=9, color="#5a4a6a")
for c, lab in [(S_main,"O1 主轴 (0,0)"), (S_II,f"轴II (0,{-a2g:.1f})"), (S_I,f"轴I (0,{-a2g-a1g:.1f})"), (S_mot,f"电机轴 ({a_belt:.0f},{-a2g-a1g:.1f})")]:
    ax.add_patch(Circle(c, 6, fc="white", ec="k", lw=1.1, zorder=6))
    ax.plot(*c, marker="+", color="k", ms=8, zorder=7)
    ax.annotate(lab, (c[0]+10, c[1]+8), fontsize=9, zorder=8)
ax.set_title("传动系统布置 前视图（XY）", fontsize=11.5)
ax.set_aspect("equal"); ax.grid(alpha=.25, lw=.5)
ax.set_xlim(-260, 620); ax.set_ylim(-420, 240)
ax.set_xlabel("X (mm)"); ax.set_ylabel("Y (mm)")
# ---- 轴向站位 ZY ----
def shaft(y, z0, z1_, lab):
    az.plot([z0, z1_], [y, y], color="#444", lw=2.5)
    az.annotate(lab, (z1_+6, y-4), fontsize=8.5)
def wheel(z, y, w, r, fc):
    az.add_patch(Rectangle((z-w/2, y-r), w, 2*r, fc=fc, ec="#555", lw=1, alpha=.85))
shaft(0, -120, 160, "主轴")
shaft(-a2g, -150, -20, "轴II")
shaft(-a2g-a1g, -200, -20, "轴I")
for z, w, r, lab in [(0,25,64,"下冲头凸轮"), (67.5,30,69,"上冲头凸轮"), (132.5,20,82,"粉筛凸轮v4")]:
    wheel(z, 0, w, r*0.35, "#c9d4dd"); az.annotate(lab, (z-16, 36), fontsize=7.5, rotation=90)
wheel(-60, 0, 45, m2*z4/2*0.35, "#cdd9c8"); az.annotate("z4", (-64, 36), fontsize=8)
wheel(-60, -a2g, 50, m2*z3/2*0.35, "#e7ddc8"); az.annotate("z3", (-64, -a2g+12), fontsize=8)
wheel(-110, -a2g, 45, m1*z2/2*0.35, "#cdd9c8"); az.annotate("z2", (-114, -a2g+30), fontsize=8)
wheel(-110, -a2g-a1g, 50, m1*z1/2*0.35, "#e7ddc8"); az.annotate("z1", (-114, -a2g-a1g+12), fontsize=8)
wheel(-165, -a2g-a1g, 35, dd2/2*0.35, "#d8cfe0"); az.annotate("大带轮", (-172, -a2g-a1g+44), fontsize=8)
wheel(-195, -a2g-a1g, 30, 140*0.35, "#e0cfcf"); az.annotate("飞轮(暂定φ280×40,\n惯量待盈亏功定)", (-228, -a2g-a1g-86), fontsize=8)
az.axvline(0, color="#aaa", lw=.6, ls="-.")
az.set_title("轴向站位（Z–Y, 半径×0.35 示意）\n齿轮箱在凸轮组对侧 Z<0", fontsize=10.5)
az.set_aspect("equal"); az.grid(alpha=.25, lw=.5)
az.set_xlim(-260, 200); az.set_ylim(-420, 240)
az.set_xlabel("Z (mm)"); az.set_ylabel("Y (mm)")
fig.suptitle(f"电机940 → 主轴{n_main:.2f} r/min   i = {i1:.2f}(带) × {i2:.2f} × {i3:.2f} = {i_act:.2f}（需求37.6）", fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(os.path.join(HERE, "参考图_传动系统.png"), dpi=170)
print("图 OK")
