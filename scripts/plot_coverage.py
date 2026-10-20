#!/usr/bin/env python3
# --------------------------------------------------------------------------------------------
# Grafico de cobertura por base del genoma mitocondrial.
# Entrada: salida de "samtools depth -a" 
#
# Uso:
#   plot_coverage.py --sample barcode01 --out barcode01_coverage.png \
#       --summary-tsv coverage_summary.tsv \
#       minimap2=barcode01_minimap2_depth.tsv ngmlr=barcode01_ngmlr_depth.tsv
#
# Dependencias: python3, pandas, matplotlib
# --------------------------------------------------------------------------------------------
import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")  # sin entorno grafico (servidor / WSL)
import matplotlib.pyplot as plt
import pandas as pd

COLORS = {"minimap2": "#2b6cb0", "ngmlr": "#dd6b20"}
DEFAULT_COLOR = "#4a5568"


def load_depth(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["chr", "pos", "depth"])
    if df.empty:
        sys.exit(f"ERROR: {path} esta vacio (samtools depth no devolvio posiciones)")
    return df


def stats(df):
    d = df["depth"]
    return {
        "ref_len": int(df["pos"].max()),
        "mean": d.mean(),
        "median": d.median(),
        "min": int(d.min()),
        "max": int(d.max()),
        "pct_1x": (d >= 1).mean() * 100,
        "pct_10x": (d >= 10).mean() * 100,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True, help="PNG de salida (se genera tambien un PDF vectorial)")
    ap.add_argument("--summary-tsv", help="TSV al que se anade una fila por alineador")
    ap.add_argument("--log", action="store_true", help="eje Y en escala symlog")
    ap.add_argument("inputs", nargs="+", metavar="HERRAMIENTA=depth.tsv")
    args = ap.parse_args()

    panels = []
    for item in args.inputs:
        if "=" not in item:
            sys.exit(f"ERROR: '{item}' debe tener el formato herramienta=archivo")
        tool, path = item.split("=", 1)
        df = load_depth(path)
        panels.append((tool, df, stats(df)))

    fig, axes = plt.subplots(len(panels), 1, figsize=(12, 3.6 * len(panels)), sharex=True, squeeze=False)
    axes = axes[:, 0]

    for ax, (tool, df, st) in zip(axes, panels):
        color = COLORS.get(tool, DEFAULT_COLOR)
        ax.fill_between(df["pos"], df["depth"], step="mid", color=color, alpha=0.85, linewidth=0)
        ax.axhline(st["mean"], color="red", ls="--", lw=1, label=f"Media: {st['mean']:.1f}x")
        ax.set_ylabel("Cobertura (x)")
        ax.set_xlim(1, st["ref_len"])
        if args.log:
            ax.set_yscale("symlog", linthresh=10)
        ax.set_ylim(bottom=0)
        ax.set_title(
            f"{tool}  |  media {st['mean']:.1f}x  ·  mediana {st['median']:.0f}x  ·  máx {st['max']}x  ·  "
            f"≥1x: {st['pct_1x']:.1f}%  ·  ≥10x: {st['pct_10x']:.1f}%",
            fontsize=10,
            loc="left",
        )
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(alpha=0.25)

    axes[-1].set_xlabel("Posición en el genoma mitocondrial (pb)")
    fig.suptitle(f"{args.sample} · cobertura por base del mtDNA (todas las lecturas alineadas, sin filtros)", fontsize=12)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    fig.savefig(os.path.splitext(args.out)[0] + ".pdf")
    plt.close(fig)

    if args.summary_tsv:
        new_file = not os.path.exists(args.summary_tsv)
        with open(args.summary_tsv, "a") as fh:
            if new_file:
                fh.write("sample\ttool\tref_len\tmean_depth\tmedian_depth\tmin_depth\tmax_depth\tpct_ge1x\tpct_ge10x\n")
            for tool, _, st in panels:
                fh.write(
                    f"{args.sample}\t{tool}\t{st['ref_len']}\t{st['mean']:.2f}\t{st['median']:.1f}\t"
                    f"{st['min']}\t{st['max']}\t{st['pct_1x']:.2f}\t{st['pct_10x']:.2f}\n"
                )

    for tool, _, st in panels:
        print(
            f"[Cobertura - {args.sample} - {tool}] media={st['mean']:.1f}x mediana={st['median']:.0f}x "
            f"min={st['min']} max={st['max']} >=1x={st['pct_1x']:.1f}% >=10x={st['pct_10x']:.1f}%"
        )
    print(f"[Cobertura - {args.sample}] Gráfico: {args.out}")


if __name__ == "__main__":
    main()
# --------------------------------------------------------------------------------------------