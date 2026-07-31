#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def discover_amplicon_folders(base_dir: Path) -> list[Path]:
	"""Return immediate subfolders that contain at least one TSV file."""
	folders: list[Path] = []
	for child in sorted(base_dir.iterdir()):
		if child.is_dir() and any(child.glob("*.tsv")):
			folders.append(child)
	return folders


def read_lengths_from_folder(folder: Path) -> pd.Series:
	"""Concatenate all TSV files in a folder and return numeric lengths from column 2."""
	all_values: list[pd.Series] = []

	for tsv_path in sorted(folder.glob("*.tsv")):
		df = pd.read_csv(tsv_path, sep="\t", header=None, usecols=[1])
		lengths = pd.to_numeric(df.iloc[:, 0], errors="coerce").dropna()
		all_values.append(lengths)

	if not all_values:
		return pd.Series(dtype="int64")

	concatenated = pd.concat(all_values, ignore_index=True)
	return concatenated.astype("int64")


def compute_stats(folder_name: str, lengths: pd.Series) -> pd.Series:
	"""Compute describe() statistics and add the folder label."""
	desc = lengths.describe()
	desc["folder"] = folder_name
	return desc


def plot_kde_panels(folder_lengths: list[tuple[str, pd.Series]], output_png: Path) -> None:
	"""Plot one KDE panel per folder in a single-column layout."""
	panel_count = len(folder_lengths)
	fig, axes = plt.subplots(panel_count, 1, figsize=(11, max(2.5 * panel_count, 4)), sharex=True)

	if panel_count == 1:
		axes = [axes]

	all_lengths = pd.concat([lengths for _, lengths in folder_lengths], ignore_index=True)
	x_min = float(all_lengths.min())
	x_max = float(all_lengths.max())
	if x_min == x_max:
		x_pad = 1.0
	else:
		x_pad = (x_max - x_min) * 0.02
	shared_xlim = (x_min - x_pad, x_max + x_pad)

	for ax, (folder_name, lengths) in zip(axes, folder_lengths):
		if lengths.nunique() > 1:
			lengths.plot.kde(ax=ax, linewidth=2)
		else:
			# KDE is undefined for a constant vector.
			ax.axvline(lengths.iloc[0], linewidth=2)
			ax.text(
				0.02,
				0.95,
				"constant value: KDE unavailable",
				transform=ax.transAxes,
				va="top",
			)

		ax.set_title(folder_name)
		ax.set_xlabel("ASV length (bp)")
		ax.set_ylabel("Density")
		ax.set_xlim(shared_xlim)
		ax.tick_params(axis="x", which="both", labelbottom=True)
		for tick_label in ax.get_xticklabels():
			tick_label.set_visible(True)

	fig.tight_layout()
	fig.savefig(output_png, dpi=300)
	plt.close(fig)


def main() -> None:
	parser = argparse.ArgumentParser(
		description=(
			"For each immediate subfolder, concatenate TSV files, compute length distribution "
			"statistics from column 2, and create one combined stats table plus KDE panels."
		)
	)
	parser.add_argument(
		"--base-dir",
		type=Path,
		default=Path(__file__).resolve().parent,
		help="Directory containing amplicon folders (default: script directory).",
	)
	parser.add_argument(
		"--stats-out",
		type=Path,
		default=Path(__file__).resolve().parent / "distribution_stats.tsv",
		help="Output TSV file for combined folder statistics.",
	)
	parser.add_argument(
		"--plot-out",
		type=Path,
		default=Path(__file__).resolve().parent / "kde_plots.png",
		help="Output PNG file for stacked KDE panels.",
	)
	args = parser.parse_args()

	base_dir = args.base_dir.resolve()
	folders = discover_amplicon_folders(base_dir)

	if not folders:
		raise SystemExit(f"No subfolders with TSV files found in: {base_dir}")

	rows: list[pd.Series] = []
	folder_lengths: list[tuple[str, pd.Series]] = []

	for folder in folders:
		lengths = read_lengths_from_folder(folder)
		if lengths.empty:
			continue
		rows.append(compute_stats(folder.name, lengths))
		folder_lengths.append((folder.name, lengths))

	if not rows:
		raise SystemExit("No numeric values found in column 2 across discovered folders.")

	stats_df = pd.DataFrame(rows)
	first_col = stats_df.pop("folder")
	stats_df.insert(0, "folder", first_col)
	stats_df.to_csv(args.stats_out, sep="\t", index=False)

	plot_kde_panels(folder_lengths, args.plot_out)

	print(f"Wrote stats table: {args.stats_out}")
	print(f"Wrote KDE panels : {args.plot_out}")


if __name__ == "__main__":
	main()
