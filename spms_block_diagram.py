# -*- coding: utf-8 -*-
"""
SPMS FOC 제어 블록도 비교 : λpm 토크 보상 없음 vs 있음
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import font_manager

try:
    _fp = "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"
    font_manager.fontManager.addfont(_fp)
    plt.rcParams["font.family"] = font_manager.FontProperties(fname=_fp).get_name()
except Exception:
    pass
plt.rcParams["axes.unicode_minus"] = False

BASE = dict(fc="#dce9f7", ec="#2c5f8d")      # 공통 블록
NEW = dict(fc="#ffe0b8", ec="#d2691e")       # 보상 시 추가되는 블록
H = 0.82


def box(ax, x, y, w, text, style=BASE, fs=10.5, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - H / 2), w, H,
                                boxstyle="round,pad=0.06,rounding_size=0.12",
                                lw=1.8, **style))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal")
    return (x, y)


def arrow(ax, p1, p2, color="#333333", ls="-", rad=0.0, lw=1.7):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=16,
                                 lw=lw, color=color, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}",
                                 shrinkA=0, shrinkB=0))


fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 8.6))
for ax in (axA, axB):
    ax.set_xlim(0, 12)
    ax.set_ylim(-0.4, 10.6)
    ax.axis("off")

# ====================================================================
# A. λpm 보상 없음
# ====================================================================
axA.set_title("A.  λpm 보상 없음  (Worst-case 설계)", fontsize=14, fontweight="bold", pad=12)
X = 6.0
ys = [9.6, 8.1, 6.6, 5.1, 3.6, 2.1]
a1 = box(axA, X, ys[0], 5.4, "토크 지령  Te*", fs=11, bold=True)
a2 = box(axA, X, ys[1], 5.4, "iq* = Te* / Kt      (id* = 0)\nKt = 상수 (20°C 값 고정)")
a3 = box(axA, X, ys[2], 5.4, "전류 PI 제어기 (d, q축)")
a4 = box(axA, X, ys[3], 5.4, "역Park + SVPWM")
a5 = box(axA, X, ys[4], 5.4, "인버터  →  모터")
a6 = box(axA, X, ys[5], 5.4, "엔코더 / 전류센서")
for p, q in zip([a1, a2, a3, a4, a5], [a2, a3, a4, a5, a6]):
    arrow(axA, (p[0], p[1] - H / 2), (q[0], q[1] + H / 2))
# 전류/각도 피드백
axA.plot([X - 2.7, X - 4.0, X - 4.0, X - 2.7],
         [ys[5], ys[5], ys[2], ys[2]], color="#2c5f8d", lw=1.7)
arrow(axA, (X - 4.0 + 0.9, ys[2]), (X - 2.7, ys[2]), color="#2c5f8d")
axA.text(X - 4.15, (ys[2] + ys[5]) / 2, "i_abc, θ", rotation=90,
         ha="right", va="center", fontsize=9.5, color="#2c5f8d")

axA.text(X, 0.25,
         "· 온도 오르면 토크가 그냥 부족해짐 (최대 -18%)\n"
         "· 속도/위치 루프가 있으면 그 루프가 자동 흡수\n"
         "· 열 폭주 위험 낮음 · 개발 간단",
         ha="center", va="center", fontsize=10.5,
         bbox=dict(fc="#eef4fa", ec="#2c5f8d", boxstyle="round,pad=0.5"))

# ====================================================================
# B. λpm 보상 있음
# ====================================================================
axB.set_title("B.  λpm 보상 있음  (적응형 설계)", fontsize=14, fontweight="bold", pad=12)
XB, XR = 4.3, 9.7
yb = [9.9, 8.55, 7.2, 5.85, 4.5, 3.15, 1.75]
b1 = box(axB, XB, yb[0], 5.0, "토크 지령  Te*", fs=11, bold=True)
b2 = box(axB, XB, yb[1], 5.0, "열 디레이팅 (토크 제한)", style=NEW, bold=True)
b3 = box(axB, XB, yb[2], 5.0, "iq* = Te* / Kt(T)\nKt = 실시간 갱신", style=NEW)
b4 = box(axB, XB, yb[3], 5.0, "전류 PI 제어기 (d, q축)")
b5 = box(axB, XB, yb[4], 5.0, "역Park + SVPWM")
b6 = box(axB, XB, yb[5], 5.0, "인버터  →  모터")
b7 = box(axB, XB, yb[6], 5.0, "엔코더 / 전류센서")
for p, q in zip([b1, b2, b3, b4, b5, b6], [b2, b3, b4, b5, b6, b7]):
    arrow(axB, (p[0], p[1] - H / 2), (q[0], q[1] + H / 2))
axB.plot([XB - 2.5, XB - 3.7, XB - 3.7, XB - 2.5],
         [yb[6], yb[6], yb[3], yb[3]], color="#2c5f8d", lw=1.7)
arrow(axB, (XB - 3.7 + 0.8, yb[3]), (XB - 2.5, yb[3]), color="#2c5f8d")
axB.text(XB - 3.85, (yb[3] + yb[6]) / 2, "i_abc, θ", rotation=90,
         ha="right", va="center", fontsize=9.5, color="#2c5f8d")

# --- 추가 블록 (오른쪽) ---
r1 = box(axB, XR, yb[1], 3.9, "권선 온도센서\n(NTC / 열모델)", style=NEW)
r2 = box(axB, XR, yb[2], 3.9, "Kt = 1.5 · P · λpm", style=NEW)
r3 = box(axB, XR, yb[3], 3.9, "λpm 추정\n(vq* - R·iq) / ωe", style=NEW)

arrow(axB, (r1[0] - 3.9 / 2, r1[1]), (XB + 5.0 / 2, r1[1]), color="#d2691e")
arrow(axB, (r2[0] - 3.9 / 2, r2[1]), (XB + 5.0 / 2, r2[1]), color="#d2691e")
arrow(axB, (r3[0], r3[1] + H / 2), (r2[0], r2[1] - H / 2), color="#d2691e")
arrow(axB, (XB + 5.0 / 2, yb[3] - 0.25), (r3[0] - 3.9 / 2, yb[3] - 0.25),
      color="#d2691e", ls="--")
axB.text((XB + 2.5 + XR - 1.95) / 2, yb[3] - 0.62, "vq*", ha="center",
         fontsize=9.5, color="#d2691e", fontweight="bold")
axB.plot([XR + 1.55, XR + 1.55], [r1[1] - H / 2, r3[1]], color="#d2691e",
         lw=1.4, ls=":")
arrow(axB, (XR + 1.55, r3[1] + 0.3), (XR + 1.55, r3[1]), color="#d2691e", lw=1.4)
axB.text(XR + 1.68, (r1[1] + r3[1]) / 2, "R(T)", fontsize=9.5, color="#d2691e",
         va="center")

axB.text(6.9, 0.2,
         "· 토크 정확 (오차 -2~5%)  · 속도루프 게인 일정\n"
         "· 대신 전류가 늘어 열 폭주 위험 ↑\n"
         "· 주황 블록 = 추가되는 블록 (디레이팅은 필수)",
         ha="center", va="center", fontsize=10.5,
         bbox=dict(fc="#fff3e3", ec="#d2691e", boxstyle="round,pad=0.5"))

fig.suptitle("λpm 토크 보상 유무에 따른 FOC 제어기 구조 비교   "
             "(전류 PI 게인은 양쪽 동일 — 플랜트가 Ls, R 뿐)",
             fontsize=14, fontweight="bold", y=0.985)
plt.tight_layout(rect=(0, 0, 1, 0.945))
plt.savefig("spms_block_diagram.png", dpi=140)
print("saved: spms_block_diagram.png")
