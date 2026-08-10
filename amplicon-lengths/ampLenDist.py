#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd


ADAPTER_FORWARD = "ACACTGACGACATGGTTCTACA"
ADAPTER_REVERSE = "TACGGTAGCAGAGACTTGGTCT"
ADAPTER_TOTAL_BP = len(ADAPTER_FORWARD) + len(ADAPTER_REVERSE)


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
	desc["median"] = desc["50%"]
	desc["folder"] = folder_name
	return desc


def primer_length_from_sequence(value: str) -> int | None:
	"""Return primer length for simple IUPAC-like strings; None when missing/invalid."""
	if pd.isna(value):
		return None
	sequence = re.sub(r"\s+", "", str(value).strip()).upper()
	if not sequence or "MISSING" in sequence or not re.fullmatch(r"[A-Z]+", sequence):
		return None
	return len(sequence)


def load_primer_offsets(assays_tsv: Path) -> dict[str, int]:
	"""Map assay_name to forward+reverse primer length (bp)."""
	assays_df = pd.read_csv(
		assays_tsv,
		sep="\t",
		usecols=["assay_name", "pcr_primer_forward", "pcr_primer_reverse"],
	)

	offsets: dict[str, int] = {}
	for row in assays_df.itertuples(index=False):
		forward_len = primer_length_from_sequence(row.pcr_primer_forward)
		reverse_len = primer_length_from_sequence(row.pcr_primer_reverse)
		offsets[row.assay_name] = (forward_len or 0) + (reverse_len or 0)

	return offsets


def plot_kde_panels(
	folder_lengths: list[tuple[str, pd.Series]],
	output_png: Path,
	primer_offsets: dict[str, int],
) -> None:
	"""Plot one KDE panel per folder in a single-column layout."""
	panel_count = len(folder_lengths)
	fig, axes = plt.subplots(panel_count, 1, figsize=(11, max(2.5 * panel_count, 4)), sharex=False)

	if panel_count == 1:
		axes = [axes]

	all_lengths_parts: list[pd.Series] = []
	for folder_name, lengths in folder_lengths:
		primer_bp = primer_offsets.get(folder_name, 0)
		all_lengths_parts.append(lengths)
		all_lengths_parts.append(lengths + primer_bp)
		all_lengths_parts.append(lengths + primer_bp + ADAPTER_TOTAL_BP)
	all_lengths = pd.concat(all_lengths_parts, ignore_index=True)
	x_min = float(all_lengths.min())
	x_max = float(all_lengths.max())
	x_lower = 50 * math.floor(x_min / 50)
	x_upper = 50 * math.ceil(x_max / 50)
	if x_lower == x_upper:
		x_lower -= 50
		x_upper += 50
	shared_xlim = (x_lower, x_upper)

	for ax, (folder_name, lengths) in zip(axes, folder_lengths):
		primer_bp = primer_offsets.get(folder_name, 0)
		distributions = [
			("ASV only", lengths),
			("ASV + primers", lengths + primer_bp),
			("ASV + primers + adapters", lengths + primer_bp + ADAPTER_TOTAL_BP),
		]

		for label, series in distributions:
			if series.nunique() > 1:
				series.plot.kde(ax=ax, linewidth=2, label=label)
			else:
				# KDE is undefined for a constant vector.
				ax.axvline(series.iloc[0], linewidth=2, label=f"{label} (constant)")

		ax.set_title(folder_name)
		ax.set_xlabel("ASV length (bp)")
		ax.set_ylabel("Density")
		ax.set_xlim(shared_xlim)
		ax.xaxis.set_major_locator(MultipleLocator(50))
		ax.legend(loc="upper right", fontsize=8)

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
	parser.add_argument(
		"--assays-tsv",
		type=Path,
		default=Path(__file__).resolve().parent.parent / "assays.tsv",
		help="Assay metadata TSV containing primer sequences.",
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

	folder_lengths.sort(key=lambda item: item[1].median())
	rows.sort(key=lambda row: row["median"])
	primer_offsets = load_primer_offsets(args.assays_tsv)

	stats_df = pd.DataFrame(rows)
	first_col = stats_df.pop("folder")
	stats_df.insert(0, "folder", first_col)
	median_col = stats_df.pop("median")
	stats_df.insert(stats_df.columns.get_loc("50%") + 1, "median", median_col)
	stats_df.to_csv(args.stats_out, sep="\t", index=False)

	plot_kde_panels(folder_lengths, args.plot_out, primer_offsets)

	print(f"Wrote stats table: {args.stats_out}")
	print(f"Wrote KDE panels : {args.plot_out}")


if __name__ == "__main__":
	main()
