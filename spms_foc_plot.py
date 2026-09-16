# -*- coding: utf-8 -*-
"""
SPMS(표면부착형 영구자석) 모터 FOC 제어 - 토크/dq 전류 관계 그래프
대상: 10 pole (극쌍수 P=5) / 12 slot 모터

핵심 식:  Te = 1.5 * P * lambda_pm * iq      (Ld = Lq 이므로 릴럭턴스 토크 = 0)
          Kt = 1.5 * P * lambda_pm           (P=5 -> Kt = 7.5 * lambda_pm)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ---- 한글 폰트 설정 -------------------------------------------------
for path in ("/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",):
    try:
        font_manager.fontManager.addfont(path)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=path).get_name()
        break
    except Exception:
        pass
plt.rcParams["axes.unicode_minus"] = False

# ---- 모터 파라미터 (본인 모터 값으로 바꾸세요) ----------------------
POLES = 10                 # 극수
P = POLES // 2             # 극쌍수 = 5
LAMBDA_PM = 0.010          # 자석 쇄교자속 [Wb]
LS = 0.2e-3                # Ld = Lq = Ls [H]
IS_MAX = 20.0              # 상전류 피크 한계 [A]
VDC = 48.0                 # DC 링크 전압 [V]
VMAX = VDC / np.sqrt(3)    # 상전압 피크 한계 (SVPWM) [V]

KT = 1.5 * P * LAMBDA_PM   # 토크 상수 [Nm/A_peak]


def solve_operating_point(we, n_grid=2001):
    """주어진 전기각속도에서 전류/전압 한계를 만족하는 최대 토크 동작점."""
    if we < 1e-6:
        return 0.0, IS_MAX, KT * IS_MAX
    v_lim = VMAX / we                       # 전압 제한원 반지름 (자속 단위)
    id_grid = np.linspace(-IS_MAX, 0.0, n_grid)
    iq_cur = np.sqrt(np.maximum(IS_MAX**2 - id_grid**2, 0.0))          # 전류 한계
    term = v_lim**2 - (LS * id_grid + LAMBDA_PM) ** 2
    iq_vol = np.where(term > 0, np.sqrt(np.maximum(term, 0.0)) / LS, -1.0)  # 전압 한계
    iq = np.minimum(iq_cur, iq_vol)
    iq = np.where(iq_vol < 0, -1.0, iq)
    if iq.max() <= 0:
        return np.nan, np.nan, np.nan
    k = int(np.argmax(iq))
    return id_grid[k], iq[k], KT * iq[k]


fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(
    f"SPMS FOC 토크-dq전류 관계  |  {POLES}극({P}극쌍) 12슬롯, "
    f"λpm={LAMBDA_PM} Wb, Ls={LS*1e3:.2f} mH, Is_max={IS_MAX:.0f} A, Vdc={VDC:.0f} V",
    fontsize=13, fontweight="bold",
)

# ====================================================================
# (1) Te = 1.5*P*λpm*iq  :  iq에 정비례하는 직선
# ====================================================================
ax = axes[0, 0]
iq = np.linspace(0, IS_MAX, 200)
for lam in (0.005, 0.010, 0.020, 0.030):
    kt = 1.5 * P * lam
    ax.plot(iq, kt * iq, lw=2, label=f"λpm={lam:.3f} Wb  (Kt={kt:.3f} Nm/A)")
ax.plot(IS_MAX, KT * IS_MAX, "o", color="crimson", ms=9, zorder=5)
ax.annotate(f"iq={IS_MAX:.0f}A → {KT*IS_MAX:.2f} Nm",
            xy=(IS_MAX, KT * IS_MAX), xytext=(-150, 10),
            textcoords="offset points", color="crimson", fontsize=10,
            arrowprops=dict(arrowstyle="->", color="crimson"))
ax.set_title("① 토크는 iq에만 비례 :  Te = 1.5·P·λpm·iq", fontsize=12)
ax.set_xlabel("q축 전류 iq [A]")
ax.set_ylabel("토크 Te [Nm]")
ax.grid(alpha=0.3)
ax.legend(fontsize=9)

# ====================================================================
# (2) dq 평면 : 전류 제한원 + 전압 제한원 + 동작 궤적
# ====================================================================
ax = axes[0, 1]
th = np.linspace(0, 2 * np.pi, 400)
ax.plot(IS_MAX * np.cos(th), IS_MAX * np.sin(th), "k--", lw=1.8,
        label=f"전류 제한원 (Is={IS_MAX:.0f} A)")

# 등토크 선 (iq = const) -> 토크는 id와 무관함을 보여줌
for te in (0.5, 1.0, 1.5):
    iq_line = te / KT
    if iq_line <= IS_MAX:
        ax.axhline(iq_line, color="gray", ls=":", lw=1)
        ax.text(-IS_MAX * 0.98, iq_line + 0.3, f"Te={te} Nm (등토크선)",
                fontsize=8, color="gray")

# 전압 제한원 (속도별) : 중심 id = -λpm/Ls
id_c = -LAMBDA_PM / LS
rpm_list = [5000, 6000, 7000, 8500]
colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(rpm_list)))
for rpm, c in zip(rpm_list, colors):
    we = 2 * np.pi * rpm / 60 * P
    r = VMAX / we / LS
    ax.plot(id_c + r * np.cos(th), r * np.sin(th), color=c, lw=1.6,
            label=f"전압 제한원 @ {rpm:,} rpm")
ax.text(-IS_MAX * 1.65, 2.0, f"전압 제한원 중심\nid = -λpm/Ls = {id_c:.0f} A\n(화면 밖)", fontsize=8, color="dimgray")

# 동작 궤적
rpm_traj = np.linspace(0, 14000, 400)
ids, iqs = [], []
for rpm in rpm_traj:
    we = 2 * np.pi * rpm / 60 * P
    d, q, _ = solve_operating_point(we)
    ids.append(d)
    iqs.append(q)
ax.plot(ids, iqs, color="crimson", lw=3, label="동작 궤적 (id=0 → 약자속)")

ax.set_title("② dq 전류 평면 : id는 토크와 무관, 전압 한계용", fontsize=12)
ax.set_xlabel("d축 전류 id [A]")
ax.set_ylabel("q축 전류 iq [A]")
ax.set_xlim(-IS_MAX * 1.7, IS_MAX * 0.4)
ax.set_ylim(-3, IS_MAX * 1.2)
ax.grid(alpha=0.3)
ax.axhline(0, color="k", lw=0.8)
ax.axvline(0, color="k", lw=0.8)
ax.legend(fontsize=8, loc="upper left")

# ====================================================================
# (3) 속도-토크 / 속도-출력
# ====================================================================
ax = axes[1, 0]
rpm = np.linspace(0, 16000, 600)
te, pw = [], []
for r in rpm:
    we = 2 * np.pi * r / 60 * P
    _, _, t = solve_operating_point(we)
    te.append(t)
    pw.append(np.nan if np.isnan(t) else t * 2 * np.pi * r / 60)
te = np.array(te)
pw = np.array(pw)

ax.plot(rpm, te, color="crimson", lw=2.5, label="토크 [Nm]")
ax.set_xlabel("회전속도 [rpm]")
ax.set_ylabel("토크 [Nm]", color="crimson")
ax.tick_params(axis="y", labelcolor="crimson")
ax.grid(alpha=0.3)

# 기저속도 표시
valid = ~np.isnan(te)
base_idx = np.where(valid & (te < te[0] * 0.995))[0]
if len(base_idx):
    rpm_base = rpm[base_idx[0]]
    ax.axvline(rpm_base, color="navy", ls="--", lw=1.5)
    ax.text(rpm_base, te[0] * 0.55, f"  기저속도 ≈ {rpm_base:,.0f} rpm\n  (여기부터 약자속)",
            fontsize=9, color="navy")

ax2 = ax.twinx()
ax2.plot(rpm, pw, color="seagreen", lw=2, ls="--", label="출력 [W]")
ax2.set_ylabel("출력 [W]", color="seagreen")
ax2.tick_params(axis="y", labelcolor="seagreen")
ax.set_title("③ 속도-토크 / 속도-출력 곡선", fontsize=12)

# ====================================================================
# (4) 속도에 따른 id, iq 지령
# ====================================================================
ax = axes[1, 1]
ids_p, iqs_p = [], []
for r in rpm:
    we = 2 * np.pi * r / 60 * P
    d, q, _ = solve_operating_point(we)
    ids_p.append(d)
    iqs_p.append(q)
ax.plot(rpm, iqs_p, color="crimson", lw=2.5, label="iq* (토크 생성)")
ax.plot(rpm, ids_p, color="navy", lw=2.5, label="id* (약자속, 토크 기여 0)")
ax.axhline(0, color="k", lw=0.8)
if len(base_idx):
    ax.axvline(rpm_base, color="gray", ls="--", lw=1.5)
    ax.text(rpm_base * 0.35, IS_MAX * 0.45, "id = 0 구간\n(MTPA)", fontsize=10, ha="center")
    ax.text(rpm_base * 1.35, IS_MAX * 0.45, "약자속 구간\n(id < 0)", fontsize=10, ha="center")
ax.set_title("④ 속도별 dq 전류 지령", fontsize=12)
ax.set_xlabel("회전속도 [rpm]")
ax.set_ylabel("전류 [A]")
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc="center right")

plt.tight_layout(rect=(0, 0, 1, 0.96))
plt.savefig("spms_foc_curves.png", dpi=140)
print("saved: spms_foc_curves.png")

# ---- 콘솔 요약 ------------------------------------------------------
print(f"\n극쌍수 P = {P},  Kt = 1.5*{P}*{LAMBDA_PM} = {KT:.4f} Nm/A")
print(f"최대 토크 (id=0, iq={IS_MAX:.0f}A) = {KT*IS_MAX:.3f} Nm")
if len(base_idx):
    print(f"기저속도 ≈ {rpm_base:,.0f} rpm")
print("\nrpm\tid[A]\tiq[A]\tTe[Nm]\tP[W]")
for r in (0, 2000, 4000, 6000, 8000, 10000, 12000, 14000):
    we = 2 * np.pi * r / 60 * P
    d, q, t = solve_operating_point(we)
    if np.isnan(t):
        print(f"{r}\t-\t-\t-\t-")
    else:
        print(f"{r}\t{d:.2f}\t{q:.2f}\t{t:.3f}\t{t*2*np.pi*r/60:.1f}")
