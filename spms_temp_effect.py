# -*- coding: utf-8 -*-
"""
SPMS FOC - 온도 상승에 따른 λpm(자석 자속) 변화와 토크 영향
대상: 10 pole (극쌍수 P=5) / 12 slot

전제:
  - R  (고정자 저항)   : 온도에 가장 크게 증가  (+0.393 %/C, 구리)
  - Ls (인덕턴스)      : 거의 불변              (공극이 지배)
  - λpm (자석 자속)    : 감소                   (-0.11 %/C, NdFeB)

  Te = 1.5 * P * λpm(T) * iq
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ---- 한글 폰트 ------------------------------------------------------
try:
    _fp = "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"
    font_manager.fontManager.addfont(_fp)
    plt.rcParams["font.family"] = font_manager.FontProperties(fname=_fp).get_name()
except Exception:
    pass
plt.rcParams["axes.unicode_minus"] = False

# ---- 모터 파라미터 (본인 모터 값으로 교체) --------------------------
POLES = 10
P = POLES // 2                 # 극쌍수 = 5
LAMBDA_PM_20 = 0.010           # 20도에서의 자석 쇄교자속 [Wb]
LS = 0.2e-3                    # Ld = Lq [H]  (온도 불변으로 취급)
R_20 = 0.050                   # 20도에서의 상저항 [ohm]
IS_MAX = 20.0                  # 상전류 피크 한계 [A]
VDC = 48.0
VMAX = VDC / np.sqrt(3)

# ---- 온도계수 -------------------------------------------------------
ALPHA_CU = 0.00393             # 구리 저항 온도계수 [1/C]
MAGNETS = {                    # 자석 종류별 Br 온도계수 [%/C]
    "NdFeB (네오디뮴)": -0.11,
    "SmCo (사마륨코발트)": -0.035,
    "Ferrite (페라이트)": -0.19,
}
ALPHA_PM = MAGNETS["NdFeB (네오디뮴)"] / 100.0   # 기본 자석


def lam(T, alpha=ALPHA_PM):
    """온도 T[C]에서의 자석 쇄교자속 [Wb]"""
    return LAMBDA_PM_20 * (1.0 + alpha * (T - 20.0))


def res(T):
    """온도 T[C]에서의 상저항 [ohm]"""
    return R_20 * (1.0 + ALPHA_CU * (T - 20.0))


def kt(T):
    """온도 T[C]에서의 토크상수 [Nm/A_peak]"""
    return 1.5 * P * lam(T)


def max_torque(rpm, T, n=1201):
    """전류/전압 한계 하에서 낼 수 있는 최대 토크 (R 포함)"""
    we = 2 * np.pi * rpm / 60.0 * P
    lm, r = lam(T), res(T)
    if we < 1e-6:
        return 0.0, IS_MAX, 1.5 * P * lm * IS_MAX
    idg = np.linspace(-IS_MAX, 0.0, n)
    # (R*id - we*Ls*iq)^2 + (R*iq + we*Ls*id + we*λpm)^2 = Vmax^2  -> iq 2차방정식
    A = r * idg
    B = we * LS
    C = we * LS * idg + we * lm
    a = B**2 + r**2
    b = 2.0 * (r * C - A * B)
    c = A**2 + C**2 - VMAX**2
    disc = b**2 - 4 * a * c
    iq_v = np.where(disc >= 0, (-b + np.sqrt(np.maximum(disc, 0.0))) / (2 * a), -1.0)
    iq_i = np.sqrt(np.maximum(IS_MAX**2 - idg**2, 0.0))          # 전류 한계
    iq = np.minimum(iq_v, iq_i)
    if iq.max() <= 0:
        return np.nan, np.nan, np.nan
    k = int(np.argmax(iq))
    return idg[k], iq[k], 1.5 * P * lm * iq[k]


fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
fig.suptitle(
    f"온도 상승에 따른 λpm 변화와 토크  |  {POLES}극({P}극쌍) 12슬롯 SPMS  "
    f"(λpm@20C={LAMBDA_PM_20} Wb, R@20C={R_20*1e3:.0f} mΩ, Ls={LS*1e3:.1f} mH 고정)",
    fontsize=13, fontweight="bold")

TEMPS = [20, 60, 100, 140, 180]
cols = plt.cm.coolwarm(np.linspace(0.05, 0.95, len(TEMPS)))

# ====================================================================
# ① Te vs iq : 온도가 오를수록 기울기(Kt)가 눕는다
# ====================================================================
ax = axes[0, 0]
iq = np.linspace(0, IS_MAX, 100)
for T, c in zip(TEMPS, cols):
    ax.plot(iq, kt(T) * iq, lw=2.4, color=c,
            label=f"{T}°C : λpm={lam(T)*1e3:.2f} mWb, Kt={kt(T):.4f}")
t20, t180 = kt(20) * IS_MAX, kt(TEMPS[-1]) * IS_MAX
ax.annotate("", xy=(IS_MAX, t180), xytext=(IS_MAX, t20),
            arrowprops=dict(arrowstyle="<->", color="k", lw=1.6))
ax.text(IS_MAX * 0.80, (t20 + t180) / 2,
        f"동일 iq={IS_MAX:.0f}A 인데\n{t20:.3f} → {t180:.3f} Nm\n({(t180/t20-1)*100:.1f} %)",
        fontsize=10, ha="right", va="center",
        bbox=dict(fc="white", ec="k", alpha=0.85))
ax.set_title("① 같은 전류를 넣어도 온도가 오르면 토크가 준다", fontsize=12)
ax.set_xlabel("q축 전류 iq [A]")
ax.set_ylabel("토크 Te [Nm]")
ax.grid(alpha=0.3)
ax.legend(fontsize=8.5, loc="upper left")

# ====================================================================
# ② 온도 vs λpm (자석 종류별) + R 대비
# ====================================================================
ax = axes[0, 1]
Tv = np.linspace(20, 180, 200)
styles = ["-", "--", "-."]
for (name, coef), ls in zip(MAGNETS.items(), styles):
    rel = (1.0 + coef / 100.0 * (Tv - 20.0)) * 100.0
    ax.plot(Tv, rel, ls, lw=2.4, label=f"λpm : {name}  ({coef:+.3f} %/°C)")
ax.plot(Tv, (res(Tv) / R_20) * 100.0, lw=2.8, color="darkorange",
        label=f"R : 구리 ({ALPHA_CU*100:+.3f} %/°C)")
ax.axhline(100, color="k", lw=0.9)
ax.plot(Tv, np.full_like(Tv, 100.0), ":", lw=2.4, color="seagreen",
        label="Ls : 인덕턴스 (거의 불변)")
ax.annotate(f"R  {res(180)/R_20*100-100:+.0f}%", xy=(180, res(180)/R_20*100),
            xytext=(-95, -5), textcoords="offset points",
            color="darkorange", fontsize=11, fontweight="bold")
ax.annotate(f"λpm(NdFeB) {(1+ALPHA_PM*160)*100-100:+.0f}%",
            xy=(180, (1 + ALPHA_PM * 160) * 100), xytext=(-145, 9),
            textcoords="offset points", color="tab:blue",
            fontsize=11, fontweight="bold")
ax.set_title("② 온도에 따른 파라미터 변화율 (20°C = 100%)", fontsize=12)
ax.set_xlabel("온도 [°C]")
ax.set_ylabel("20°C 대비 [%]")
ax.grid(alpha=0.3)
ax.legend(fontsize=8.5, loc="lower left")

# ====================================================================
# ③ 토크-속도 곡선 (온도별)
# ====================================================================
ax = axes[1, 0]
rpm = np.linspace(0, 11000, 320)
for T, c in zip(TEMPS, cols):
    te = np.array([max_torque(r, T)[2] for r in rpm])
    ax.plot(rpm, te, lw=2.4, color=c, label=f"{T}°C")
    ok = ~np.isnan(te)
    if ok.any():
        ax.plot(rpm[ok][-1], te[ok][-1], "o", color=c, ms=5)
ax.set_title("③ 토크-속도 곡선 : 온도가 오르면 아래로 눌린다", fontsize=12)
ax.set_xlabel("회전속도 [rpm]")
ax.set_ylabel("최대 토크 [Nm]")
ax.grid(alpha=0.3)
ax.legend(fontsize=9, title="권선/자석 온도", title_fontsize=9)
ax.text(300, kt(20) * IS_MAX * 0.18,
        "저속(정토크) : λpm 감소만큼 토크 감소\n"
        "고속(약자속) : λpm이 줄어 오히려 속도는 더 나감",
        fontsize=9, bbox=dict(fc="lightyellow", ec="gray", alpha=0.9))

# ====================================================================
# ④ 같은 토크를 유지하려면? -> 전류 증가 -> 동손 폭증
# ====================================================================
ax = axes[1, 1]
TE_TARGET = 1.0                                   # 유지하고 싶은 토크 [Nm]
iq_need = TE_TARGET / (1.5 * P * lam(Tv))         # 필요한 iq
p_cu = 1.5 * iq_need**2 * res(Tv)                 # 3상 동손 [W]

ax.plot(Tv, iq_need, lw=2.6, color="crimson", label=f"필요 iq (Te={TE_TARGET} Nm 유지)")
ax.set_xlabel("온도 [°C]")
ax.set_ylabel("필요 전류 iq [A]", color="crimson")
ax.tick_params(axis="y", labelcolor="crimson")
ax.grid(alpha=0.3)

ax2 = ax.twinx()
ax2.plot(Tv, p_cu, lw=2.6, ls="--", color="darkorange", label="동손 1.5·iq²·R [W]")
ax2.set_ylabel("동손 [W]", color="darkorange")
ax2.tick_params(axis="y", labelcolor="darkorange")

i0, p0 = iq_need[0], p_cu[0]
i1, p1 = iq_need[-1], p_cu[-1]
ax.text(0.04, 0.72,
        f"20°C → 180°C\n"
        f"  iq   : {i0:.1f} → {i1:.1f} A  ({i1/i0:.2f}배)\n"
        f"  동손 : {p0:.1f} → {p1:.1f} W  ({p1/p0:.2f}배)\n\n"
        f"발열↑ → 온도↑ → 다시 발열↑\n(열 폭주)",
        transform=ax.transAxes, fontsize=10,
        bbox=dict(fc="mistyrose", ec="crimson", alpha=0.9))
ax.set_title("④ 토크를 유지하려 들면 발열이 제곱으로 늘어난다", fontsize=12)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=9, loc="lower right")

plt.tight_layout(rect=(0, 0, 1, 0.955))
plt.savefig("spms_temp_torque.png", dpi=140)
print("saved: spms_temp_torque.png\n")

# ---- 콘솔 요약표 ----------------------------------------------------
print("온도[C]\tλpm[mWb]\tKt[Nm/A]\tR[mOhm]\tTe@20A[Nm]\t토크비[%]\t필요iq@1Nm[A]\t동손@1Nm[W]")
for T in (20, 40, 60, 80, 100, 120, 140, 160, 180):
    lm, r, k = lam(T), res(T), kt(T)
    te = k * IS_MAX
    iqn = TE_TARGET / k
    pcu = 1.5 * iqn**2 * r
    print(f"{T}\t{lm*1e3:.3f}\t{k:.5f}\t{r*1e3:.2f}\t{te:.4f}\t"
          f"{te/(kt(20)*IS_MAX)*100:.1f}\t{iqn:.2f}\t{pcu:.2f}")
