from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class NativeJutilsPlotMixin:
    """PNG rendering for Native Splicing Pattern Explorer outputs."""

    def _plot_native_jutils_tables(self, tables: dict[str, pd.DataFrame], plots_dir: Path) -> list[Path]:
        generated: list[Path] = []
        candidate_summary = tables.get("candidate_gene_as_summary", pd.DataFrame())
        if not candidate_summary.empty:
            generated.append(
                self._plot_candidate_gene_as_heatmap(
                    candidate_summary,
                    plots_dir / "native_candidate_gene_as_summary_heatmap.png",
                )
            )
        generated.extend(self._plot_native_per_comparison_psi_heatmaps(tables, plots_dir / "per_comparison"))
        generated.extend(self._plot_native_per_comparison_psi_pcas(tables, plots_dir / "per_comparison"))
        dpsi_matrix = tables.get("dpsi_matrix", pd.DataFrame())
        if not dpsi_matrix.empty:
            generated.append(
                self._plot_matrix_heatmap(
                    dpsi_matrix,
                    "Native AS dPSI pattern heatmap",
                    plots_dir / "native_as_dpsi_pattern_heatmap.png",
                    center_zero=True,
                    cluster_rows=True,
                )
            )
        psi_matrix = tables.get("sample_level_psi_matrix", pd.DataFrame())
        if not psi_matrix.empty:
            generated.append(
                self._plot_matrix_heatmap(
                    psi_matrix,
                    "Native sample-level PSI heatmap",
                    plots_dir / "native_sample_level_psi_heatmap.png",
                    center_zero=False,
                    cluster_rows=True,
                )
            )
        psi_zscore_matrix = tables.get("sample_level_psi_zscore_matrix", pd.DataFrame())
        if not psi_zscore_matrix.empty:
            generated.append(
                self._plot_matrix_heatmap(
                    psi_zscore_matrix,
                    "Sample-level PSI heatmap (row Z-score)",
                    plots_dir / "native_sample_level_psi_zscore_heatmap.png",
                    center_zero=True,
                    cluster_rows=True,
                    cluster_columns=True,
                )
            )
        similarity = tables.get("pairwise_as_similarity", pd.DataFrame())
        if not similarity.empty:
            generated.append(self._plot_similarity_heatmap(similarity, "weighted_cosine_similarity", "AS weighted cosine similarity", plots_dir / "native_as_similarity_heatmap.png"))
            generated.append(self._plot_similarity_heatmap(similarity, "direction_similarity", "AS direction agreement", plots_dir / "native_as_direction_similarity_heatmap.png"))
        sample_similarity = tables.get("sample_psi_spearman_matrix", pd.DataFrame())
        sample_overlap = tables.get("sample_psi_shared_event_count_matrix", pd.DataFrame())
        if not sample_similarity.empty:
            generated.append(
                self._plot_square_matrix_heatmap(
                    sample_similarity,
                    "Sample PSI similarity (Spearman)",
                    plots_dir / "native_sample_psi_similarity_heatmap.png",
                    value_label="Spearman rho",
                    minimum=-1.0,
                    maximum=1.0,
                    cmap="RdBu_r",
                )
            )
        if not sample_overlap.empty:
            maximum = pd.to_numeric(
                sample_overlap.drop(columns=[sample_overlap.columns[0]]).stack(), errors="coerce"
            ).max()
            generated.append(
                self._plot_square_matrix_heatmap(
                    sample_overlap,
                    "Shared measured AS events per sample pair",
                    plots_dir / "native_sample_psi_shared_event_counts.png",
                    value_label="shared events",
                    minimum=0.0,
                    maximum=max(float(maximum) if pd.notna(maximum) else 1.0, 1.0),
                    cmap="YlGnBu",
                    integer_labels=True,
                )
            )
        sample_pca = tables.get("sample_psi_pca", pd.DataFrame())
        if not sample_pca.empty:
            generated.append(
                self._plot_native_pca(
                    sample_pca,
                    "sample_id",
                    "Native sample-level PSI PCA",
                    plots_dir / "native_sample_psi_pca.png",
                    color_column="sample_group",
                )
            )
        pattern_clusters = tables.get("psi_pattern_clusters", pd.DataFrame())
        if not pattern_clusters.empty:
            generated.append(
                self._plot_psi_pattern_cluster_heatmap(
                    pattern_clusters,
                    plots_dir / "native_sample_psi_pattern_clusters.png",
                )
            )
        replicate_separation = tables.get("replicate_group_separation", pd.DataFrame())
        if not replicate_separation.empty:
            generated.append(
                self._plot_replicate_group_separation(
                    replicate_separation,
                    plots_dir / "native_replicate_group_separation.png",
                )
            )
        return [path for path in generated if path.exists()]

    def _plot_square_matrix_heatmap(
        self,
        frame: pd.DataFrame,
        title: str,
        path: Path,
        *,
        value_label: str,
        minimum: float,
        maximum: float,
        cmap: str,
        integer_labels: bool = False,
    ) -> Path:
        if frame.empty or frame.shape[1] <= 1:
            return path
        labels = frame.iloc[:, 0].astype(str).tolist()
        numeric = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        numeric.columns = numeric.columns.astype(str)
        numeric = numeric.reindex(columns=labels)
        data = numeric.to_numpy(dtype=float)
        size = max(6.5, min(16, 3.5 + len(labels) * 0.75))
        fig, ax = plt.subplots(figsize=(size, size))
        image = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, vmin=minimum, vmax=maximum)
        ax.set_title(title, fontweight="bold", pad=14)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, ha="center", fontsize=8)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        if len(labels) <= 18:
            for row_index in range(len(labels)):
                for column_index in range(len(labels)):
                    value = data[row_index, column_index]
                    if np.isfinite(value):
                        label = f"{int(round(value))}" if integer_labels else f"{value:.2f}"
                        ax.text(column_index, row_index, label, ha="center", va="center", fontsize=7)
        fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03, label=value_label)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return path

    def _plot_psi_pattern_cluster_heatmap(self, frame: pd.DataFrame, path: Path) -> Path:
        if frame.empty or "psi_pattern_cluster" not in frame.columns:
            return path
        metadata = {
            "event_matching_key",
            "gene_symbol",
            "event_type",
            "event_id",
            "psi_pattern_cluster",
        }
        value_columns = [
            column
            for column in frame.columns
            if column not in metadata and pd.api.types.is_numeric_dtype(frame[column])
        ]
        if not value_columns:
            return path
        working = frame.sort_values(["psi_pattern_cluster", "event_type", "gene_symbol"], na_position="last").head(120)
        data = working[value_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        fig, ax = plt.subplots(
            figsize=(max(9, min(18, 3 + len(value_columns) * 0.75)), max(6, min(18, 2 + len(working) * 0.16)))
        )
        image = ax.imshow(np.ma.masked_invalid(data), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        ax.set_title("Sample-level PSI trajectory clusters", fontweight="bold", pad=14)
        ax.set_xticks(range(len(value_columns)))
        ax.set_xticklabels(value_columns, rotation=90, fontsize=7)
        labels = [
            f"C{row.psi_pattern_cluster} | {str(getattr(row, 'gene_symbol', ''))[:18]} | {getattr(row, 'event_type', '')}"
            for row in working.itertuples(index=False)
        ]
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=6)
        fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02, label="row Z-score PSI")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return path

    def _plot_replicate_group_separation(self, frame: pd.DataFrame, path: Path) -> Path:
        working = frame.copy()
        if working.empty or "comparison_id" not in working.columns:
            return path
        within = pd.to_numeric(working["median_within_group_spearman"], errors="coerce")
        between = pd.to_numeric(working["median_between_group_spearman"], errors="coerce")
        x = np.arange(len(working))
        fig, ax = plt.subplots(figsize=(max(8, len(working) * 1.35), 5.5))
        width = 0.36
        ax.bar(x - width / 2, within, width, label="within-group replicates", color="#20a39e")
        ax.bar(x + width / 2, between, width, label="between groups", color="#f19c79")
        ax.axhline(0, color="#6d7b86", linewidth=0.8)
        ax.set_ylim(-1, 1)
        ax.set_ylabel("Median Spearman PSI similarity")
        ax.set_title("Replicate/group PSI similarity QC", fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(working["comparison_id"].astype(str), rotation=35, ha="right")
        ax.legend(frameon=False, loc="upper right")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return path

    def _plot_candidate_gene_as_heatmap(self, summary: pd.DataFrame, path: Path) -> Path:
        required = {"gene_symbol", "comparison_display_name", "best_dPSI", "n_significant_as_events"}
        if summary.empty or not required.issubset(summary.columns):
            return path
        working = summary.copy()
        working["best_dPSI_numeric"] = pd.to_numeric(working["best_dPSI"], errors="coerce")
        working["n_significant_numeric"] = pd.to_numeric(working["n_significant_as_events"], errors="coerce").fillna(0)
        pivot = working.pivot_table(
            index="gene_symbol",
            columns="comparison_display_name",
            values="best_dPSI_numeric",
            aggfunc="first",
            dropna=False,
        )
        counts = working.pivot_table(
            index="gene_symbol",
            columns="comparison_display_name",
            values="n_significant_numeric",
            aggfunc="max",
            fill_value=0,
            dropna=False,
        ).reindex(index=pivot.index, columns=pivot.columns, fill_value=0)
        order = (
            pd.DataFrame(
                {
                    "n_total_significant": counts.sum(axis=1),
                    "max_abs_dPSI": pivot.abs().max(axis=1).fillna(0),
                    "gene": pivot.index.astype(str),
                },
                index=pivot.index,
            )
            .sort_values(["n_total_significant", "max_abs_dPSI", "gene"], ascending=[False, False, True])
            .index
        )
        pivot = pivot.loc[order].head(80)
        counts = counts.loc[pivot.index]
        data = pivot.to_numpy(dtype=float)
        masked = np.ma.masked_invalid(data)
        vmax = float(np.nanmax(np.abs(data))) if np.isfinite(data).any() else 0.2
        vmax = max(vmax, 0.2)
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#d8d8d8")
        width = max(8, min(18, 2.8 + pivot.shape[1] * 1.4))
        height = max(5.5, min(18, 2.2 + pivot.shape[0] * 0.34))
        fig, ax = plt.subplots(figsize=(width, height))
        image = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
        ax.set_title("Candidate gene AS summary heatmap\nselected genes retained; gray = no significant AS under active cutoffs")
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels([str(col) for col in pivot.columns], rotation=35, ha="right", fontsize=8)
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels([str(item)[:45] for item in pivot.index], fontsize=8)
        for row_idx in range(pivot.shape[0]):
            for col_idx in range(pivot.shape[1]):
                count_value = int(counts.iloc[row_idx, col_idx]) if pd.notna(counts.iloc[row_idx, col_idx]) else 0
                if count_value > 0:
                    ax.text(col_idx, row_idx, f"n={count_value}", ha="center", va="center", fontsize=6, color="black")
        fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02, label="dominant significant AS dPSI")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return path

    def _plot_native_per_comparison_psi_pcas(self, tables: dict[str, pd.DataFrame], plots_dir: Path) -> list[Path]:
        frame = tables.get("per_comparison_sample_psi_pca", pd.DataFrame())
        if frame.empty or "comparison_id" not in frame.columns:
            return []
        generated: list[Path] = []
        for comparison_id, subset in frame.groupby("comparison_id", sort=False):
            display_name = (
                str(subset["comparison_display_name"].dropna().iloc[0])
                if "comparison_display_name" in subset.columns and subset["comparison_display_name"].notna().any()
                else str(comparison_id)
            )
            path = plots_dir / f"{self._safe_output_name(display_name)}_sample_psi_pca.png"
            generated.append(
                self._plot_native_pca(
                    subset,
                    "sample_id",
                    f"Native sample PSI PCA: {display_name}",
                    path,
                    color_column="sample_group",
                )
            )
        return generated

    def _plot_native_per_comparison_psi_heatmaps(self, tables: dict[str, pd.DataFrame], plots_dir: Path) -> list[Path]:
        psi_long = tables.get("sample_level_psi_long", pd.DataFrame())
        events = tables.get("unified_as_event_table", pd.DataFrame())
        if psi_long.empty or "comparison_id" not in psi_long.columns:
            return []
        generated: list[Path] = []
        event_meta = pd.DataFrame()
        if not events.empty and {"comparison_id", "event_matching_key"}.issubset(events.columns):
            meta_cols = [
                column
                for column in [
                    "comparison_id",
                    "comparison_display_name",
                    "event_matching_key",
                    "gene_symbol",
                    "event_type",
                    "FDR",
                    "dPSI",
                    "abs_dPSI",
                    "native_jutils_selection_score",
                ]
                if column in events.columns
            ]
            event_meta = events[meta_cols].drop_duplicates(["comparison_id", "event_matching_key"])
        for comparison_id, subset in psi_long.groupby("comparison_id", sort=False):
            pivot = subset.pivot_table(
                index="event_matching_key",
                columns="sample_id",
                values="PSI",
                aggfunc="mean",
            )
            if pivot.empty:
                continue
            if not event_meta.empty:
                meta = event_meta.loc[event_meta["comparison_id"].astype(str) == str(comparison_id)].copy()
                if not meta.empty:
                    sort_columns = [column for column in ["native_jutils_selection_score", "abs_dPSI", "FDR"] if column in meta.columns]
                    sort_ascending = [False if column != "FDR" else True for column in sort_columns]
                    order_source = meta.sort_values(sort_columns, ascending=sort_ascending, na_position="last") if sort_columns else meta
                    order = order_source["event_matching_key"].astype(str).tolist()
                    pivot = pivot.loc[[key for key in order if key in pivot.index]]
                    label_lookup = {
                        str(row.event_matching_key): f"{str(row.gene_symbol)[:18]} | {str(row.event_type)}"
                        for row in meta[["event_matching_key", "gene_symbol", "event_type"]].itertuples(index=False)
                    }
                    pivot.index = [label_lookup.get(str(index), str(index)[:50]) for index in pivot.index]
            display_name = str(subset["comparison_display_name"].dropna().iloc[0]) if "comparison_display_name" in subset.columns and subset["comparison_display_name"].notna().any() else str(comparison_id)
            path = plots_dir / f"{self._safe_output_name(display_name)}_sample_psi_heatmap.png"
            generated.append(
                self._plot_matrix_heatmap(
                    pivot.reset_index().rename(columns={"index": "event"}),
                    f"Native sample PSI heatmap: {display_name}",
                    path,
                    center_zero=False,
                    cluster_rows=True,
                )
            )
            clustered_path = plots_dir / f"{self._safe_output_name(display_name)}_sample_psi_heatmap_clustered.png"
            generated.append(
                self._plot_matrix_heatmap(
                    pivot.reset_index().rename(columns={"index": "event"}),
                    f"Native clustered sample PSI heatmap: {display_name}",
                    clustered_path,
                    center_zero=False,
                    cluster_rows=True,
                    cluster_columns=True,
                )
            )
        return generated

    def _plot_matrix_heatmap(
        self,
        frame: pd.DataFrame,
        title: str,
        path: Path,
        *,
        center_zero: bool,
        cluster_rows: bool = True,
        cluster_columns: bool = False,
    ) -> Path:
        working = frame.copy()
        if working.empty or working.shape[1] <= 1:
            return path
        index_col = working.columns[0]
        working = working.set_index(index_col)
        numeric = working.apply(pd.to_numeric, errors="coerce")
        row_score = numeric.abs().max(axis=1).fillna(0)
        numeric = numeric.loc[row_score.sort_values(ascending=False).index].head(120)
        if cluster_rows and len(numeric) > 2:
            numeric = self._cluster_order_matrix_rows(numeric)
        if cluster_columns and numeric.shape[1] > 2:
            numeric = self._cluster_order_matrix_rows(numeric.T).T
        numeric = numeric.head(80)
        data = numeric.to_numpy(dtype=float)
        masked = np.ma.masked_invalid(data)
        width = max(8, min(18, 2 + numeric.shape[1] * 1.2))
        height = max(5, min(18, 2 + numeric.shape[0] * 0.22))
        fig, ax = plt.subplots(figsize=(width, height))
        if center_zero:
            vmax = float(np.nanmax(np.abs(data))) if np.isfinite(data).any() else 1.0
            vmax = max(vmax, 0.2)
            image = ax.imshow(masked, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        else:
            image = ax.imshow(masked, aspect="auto", cmap="magma", vmin=0, vmax=1)
        ax.set_title(title)
        ax.set_xticks(range(numeric.shape[1]))
        ax.set_xticklabels([str(col) for col in numeric.columns], rotation=45, ha="right", fontsize=8)
        y_labels = [str(item)[:70] for item in numeric.index]
        ax.set_yticks(range(numeric.shape[0]))
        ax.set_yticklabels(y_labels, fontsize=6)
        fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return path

    def _cluster_order_matrix_rows(self, numeric: pd.DataFrame) -> pd.DataFrame:
        working = numeric.apply(pd.to_numeric, errors="coerce")
        filled = working.fillna(working.mean(axis=0)).fillna(0.0)
        if len(filled) <= 2:
            return working
        try:
            from scipy.cluster.hierarchy import leaves_list, linkage
            from scipy.spatial.distance import pdist

            distances = pdist(filled.to_numpy(dtype=float), metric="cityblock")
            if len(distances) and np.isfinite(distances).all():
                ordered_index = filled.index[leaves_list(linkage(distances, method="weighted"))]
                return working.loc[ordered_index]
        except Exception:
            pass
        # Stable fallback when SciPy is unavailable or the matrix is degenerate.
        order_frame = pd.DataFrame(
            {
                "row": filled.index.astype(str),
                "mean_abs_signal": filled.abs().mean(axis=1).to_numpy(dtype=float),
                "max_abs_signal": filled.abs().max(axis=1).to_numpy(dtype=float),
            },
            index=filled.index,
        )
        ordered_index = order_frame.sort_values(["max_abs_signal", "mean_abs_signal", "row"], ascending=[False, False, True]).index
        return working.loc[ordered_index]

    def _plot_similarity_heatmap(self, frame: pd.DataFrame, metric: str, title: str, path: Path) -> Path:
        labels = sorted(set(frame["comparison_A_display_name"].astype(str)) | set(frame["comparison_B_display_name"].astype(str)))
        matrix = pd.DataFrame(np.nan, index=labels, columns=labels, dtype=float)
        for label in labels:
            matrix.loc[label, label] = 1.0
        for _, row in frame.iterrows():
            left = str(row["comparison_A_display_name"])
            right = str(row["comparison_B_display_name"])
            value = pd.to_numeric(row.get(metric), errors="coerce")
            matrix.loc[left, right] = value
            matrix.loc[right, left] = value
        fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.2), max(5, len(labels) * 1.0)))
        image = ax.imshow(np.ma.masked_invalid(matrix.to_numpy(dtype=float)), cmap="viridis", vmin=-1 if "cosine" in metric else 0, vmax=1)
        ax.set_title(title)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        for row_idx in range(len(labels)):
            for col_idx in range(len(labels)):
                value = matrix.iloc[row_idx, col_idx]
                if pd.notna(value):
                    ax.text(col_idx, row_idx, f"{float(value):.2f}", ha="center", va="center", color="white" if float(value) < 0.45 else "black", fontsize=8)
        fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return path

    def _plot_native_pca(self, frame: pd.DataFrame, label_column: str, title: str, path: Path, *, color_column: str | None = None) -> Path:
        if frame.empty or "PC1" not in frame.columns or "PC2" not in frame.columns:
            return path
        fig, ax = plt.subplots(figsize=(8, 6))
        label_offsets: dict[object, tuple[int, int]] = {}
        if color_column and color_column in frame.columns:
            groups = frame[color_column].astype(str).fillna("")
            unique_groups = sorted(groups.unique().tolist())
            cmap = plt.get_cmap("tab10")
            for idx, group in enumerate(unique_groups):
                subset = frame.loc[groups == group]
                ax.scatter(subset["PC1"], subset["PC2"], s=70, label=group, color=cmap(idx % 10), edgecolor="black", linewidth=0.4)
                ordered_indices = subset.sort_values(["PC2", "PC1"]).index.tolist()
                offsets = [-12, 0, 12, -20, 20]
                for rank, row_index in enumerate(ordered_indices):
                    label_offsets[row_index] = (7, offsets[rank % len(offsets)])
        else:
            ax.scatter(frame["PC1"], frame["PC2"], s=85, color="#4C78A8", edgecolor="black", linewidth=0.5)
        for row_index, row in frame.iterrows():
            label = str(row.get(label_column, ""))[:38]
            offset = label_offsets.get(row_index, (7, 0))
            ax.annotate(
                label,
                (float(row["PC1"]), float(row["PC2"])),
                xytext=offset,
                textcoords="offset points",
                fontsize=8,
                va="center",
                arrowprops={"arrowstyle": "-", "color": "#7b8794", "linewidth": 0.45},
            )
        pc1_var = pd.to_numeric(frame.get("PC1_variance_explained"), errors="coerce").dropna()
        pc2_var = pd.to_numeric(frame.get("PC2_variance_explained"), errors="coerce").dropna()
        x_label = f"PC1 ({float(pc1_var.iloc[0]) * 100:.1f}% var)" if not pc1_var.empty else "PC1"
        y_label = f"PC2 ({float(pc2_var.iloc[0]) * 100:.1f}% var)" if not pc2_var.empty else "PC2"
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.axhline(0, color="#999999", linewidth=0.8)
        ax.axvline(0, color="#999999", linewidth=0.8)
        if color_column and color_column in frame.columns:
            ax.legend(title=color_column, bbox_to_anchor=(1.02, 1), loc="upper left")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return path
