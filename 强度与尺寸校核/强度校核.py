# -*- coding: utf-8 -*-
"""15吨压片机 关键部件 初次强度与尺寸校核（解析法/材料力学公式）
   假设: 静力分析、忽略摩擦与自重、压片力 Fp=150kN、材料 45钢调质(σs≈355,σb≈600)。
"""
import numpy as np
pi=np.pi
Fp=150e3          # 压片力 N
sigma_s=355.; sigma_b=600.      # 45钢调质 MPa
sig_allow=235.    # 许用弯曲/拉压 [σ]=σs/1.5 MPa
tau_pin=90.       # 销许用剪切 [τ] MPa(保守)
sig_bear=1.5*sig_allow  # 许用挤压 MPa
tauT_shaft=40.    # 轴扭转许用 [τT] MPa(含弯曲影响,设计用)
sigH_gear=550.; sigF_gear=250.  # 齿轮 45钢调质许用 MPa
sigH_cam=1200.    # 凸轮/滚子 表面硬化 许用接触 MPa
print("="*64)

# ========== 1. 上冲头 保压力链 静力分析(保压位坐标 取自结果.json) ==========
E=np.array([0,335.]); C=np.array([8.304,409.539]); D=np.array([205.848,480.789])
B=np.array([248.073,58.967]); O3=np.array([140,457.039]); O2=np.array([100,35.]); P=np.array([68.039,-43.762])
def unit(v): return v/np.linalg.norm(v)
# 小连杆 C-E (二力杆): 竖直分量平衡 150kN
uEC=unit(C-E); F2=Fp/uEC[1]          # 小连杆轴力
FC=-F2*uEC                            # 小连杆作用于上摆杆C点的力(沿C->E)
# 上摆杆 O3: 对O3取矩 = 0, 解长连杆轴力 F1
rC=C-O3; uDB=unit(B-D); rD=D-O3
MC=rC[0]*FC[1]-rC[1]*FC[0]
m_per=rD[0]*uDB[1]-rD[1]*uDB[0]
F1=-MC/m_per                          # 长连杆轴力
FB=F1*uDB                             # 长连杆作用于下摆杆B点的力
# 下摆杆 O2: 对O2取矩, 滚子力(设垂直于滚子臂,上界)
rB=B-O2; MB=rB[0]*FB[1]-rB[1]*FB[0]
Lr=85.; Froll=abs(MB)/Lr              # 凸轮接触力(上界)
print("【1 上冲头保压力链】")
print(" 小连杆轴力 F2=%.1f kN | 长连杆轴力 F1=%.1f kN | 凸轮接触力 Froll≈%.1f kN"%(F2/1e3,F1/1e3,Froll/1e3))

# ========== 2. 销 校核(双剪+挤压) ==========
def pin_check(F,d,Lbear,name):
    A2=2*pi/4*d**2; tau=F/A2
    sb=F/(d*Lbear)
    dmin=np.sqrt(F/(2*tau_pin*pi/4))
    print("  %s: F=%.0fkN d=φ%g 双剪τ=%.0fMPa(许%g) 挤压σ=%.0fMPa(许%g) -> %s, 需≥φ%.0f"%(
        name,F/1e3,d,tau,tau_pin,sb,sig_bear, "合格" if tau<=tau_pin and sb<=sig_bear else "不合格",dmin))
print("【2 高载销(现为d8, 承力杆销)】")
pin_check(F2,8,15,"C销(小连杆-冲头)")
pin_check(F1,8,18,"D销(上摆杆-长连杆)")
pin_check(F1,8,18,"B销(长连杆-下摆杆)")

# ========== 3. 杆件 校核(拉压, 长连杆加压杆稳定) ==========
def bar_check(F,sec,name,L=None):
    A=sec; sig=F/A; dmin=F/sig_allow
    s="拉压σ=%.0fMPa(许%g)->%s 需A≥%.0fmm²"%(sig,sig_allow,"合格" if sig<=sig_allow else "不合格",dmin)
    print("  %s: F=%.0fkN 现截面A=%.0fmm² %s"%(name,F/1e3,A,s))
print("【3 承力杆件(截面为建模辅助初值)】")
bar_check(F2,30*20,"小连杆(设30×20)")
bar_check(F1,30*20,"长连杆(设30×20)")
# 长连杆压杆稳定(两端铰, 欧拉)
E_mod=2.06e5; b,h=20,30; I=b*h**3/12; L1=423.93
Pcr=pi**2*E_mod*I/L1**2
print("  长连杆压杆稳定: I=%.0fmm⁴ 欧拉临界力Pcr=%.1fkN vs 轴力%.1fkN -> %s"%(I,Pcr/1e3,F1/1e3,"稳定" if Pcr>2*F1 else "失稳/裕度不足"))

# ========== 4. 冲头 压缩 ==========
d_punch=34.; A=pi/4*d_punch**2; sig=Fp/A
print("【4 冲头φ34 压缩】 σ=%.0fMPa 安全系数=%.2f -> %s"%(sig,sigma_s/sig,"合格" if sig<=sig_allow else "偏高"))

# ========== 5. 主轴 扭转(T=746.5N·m, 现φ20) ==========
T=746.5e3  # N·mm
d_shaft=20.; Wp=pi/16*d_shaft**3; tau=T/Wp
dmin=(16*T/(pi*tauT_shaft))**(1/3)
print("【5 主轴 扭转】 T=746.5N·m d=φ20 τ=%.0fMPa(许%g) -> %s, 需≥φ%.0f"%(tau,tauT_shaft,"合格" if tau<=tauT_shaft else "不合格",dmin))

# ========== 6. 齿轮 低速级(主轴大齿轮z4, 受最大转矩) ==========
m=3.; z3=17; z4=60; b=50.; d3=m*z3; d4=m*z4; u=z4/z3
Ft=2*T/d4                      # 圆周力 N
sigF=Ft/(b*m*0.38)             # Lewis弯曲(Y≈0.38) MPa
sigH=189.8*np.sqrt(Ft*(u+1)/(b*d3*u))  # 接触(简化, ZE·ZH≈189.8)
print("【6 齿轮低速级 m3 z17/60 b50】 Ft=%.0fN 弯曲σF=%.0fMPa(许%g)->%s 接触σH=%.0fMPa(许%g)->%s"%(
    Ft,sigF,sigF_gear,"合格" if sigF<=sigF_gear else "不合格",sigH,sigH_gear,"合格" if sigH<=sigH_gear else "不合格"))

# ========== 7. 凸轮 接触(赫兹线接触, 取上冲头凸轮保压力) ==========
bw=15.; Fpr=Froll/bw           # N/mm 单位宽载荷
rho_cam=23.; rho_roll=12.; nu=0.3
sigHc=np.sqrt(Fpr*E_mod*(1/rho_cam+1/rho_roll)/(2*pi*(1-nu**2)))
print("【7 上冲头凸轮 接触(若保压力经凸轮)】 σH=%.0fMPa(许%g) -> %s"%(sigHc,sigH_cam,"合格" if sigHc<=sigH_cam else "严重超限"))

# ========== 8. 飞轮 转动惯量(盈亏功法) ==========
sc=0.010   # 假设压实有效行程 10mm, 力线性升到150kN
Wr=0.5*Fp*sc           # 压片功 J
phi_press=np.radians(80)   # 加压相位角程
Md=Wr/(2*pi); Mr=Wr/phi_press
# 盈亏功(梯形累积): 近休累积 +Md*(2π-phi_press)... 解析: ΔW=Md*(2π-phi_press)
dW=Md*(2*pi-phi_press)
delta=0.05
for shaft,nn in [("主轴",25),("轴II",87.9),("轴I",398.3),("电机轴",940)]:
    w=nn*2*pi/60; J=dW/(delta*w**2)
    # 实心铸铁盘 J=0.5 m R²; 取R=0.15(φ300) 反算需厚度
    print("  飞轮装于%s(n=%g): 需 J=%.3f kg·m²"%(shaft,nn,J), end="")
    if shaft=="轴I":
        R=0.15; m_disk=2*J/R**2; rho=7200; th=m_disk/(rho*pi*R**2)
        print("  -> φ300实心盘需厚%.0fmm/质量%.0fkg"%(th*1000,m_disk))
    else: print()
print("  盈亏功 ΔW=%.0f J (压片功 Wr=%.0f J, 假设压实行程10mm)"%(dW,Wr))
print("="*64)
