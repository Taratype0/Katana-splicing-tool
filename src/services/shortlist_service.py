from __future__ import annotations

import pandas as pd


class ShortlistService:
    KEEP_CLASSES = {"A_only", "B_only", "opposite_direction", "shared_same_direction", "same_direction_large_delta"}

    def rank_candidates(self, program_events: pd.DataFrame) -> pd.DataFrame:
        if program_events.empty:
            return pd.DataFrame()
        frame = program_events.copy()
        frame = frame[frame["class"].isin(self.KEEP_CLASSES)].copy()
        if frame.empty:
            return frame
        if "comparison_A_name" not in frame.columns:
            frame["comparison_A_name"] = frame["comparison_A"]
        if "comparison_B_name" not in frame.columns:
            frame["comparison_B_name"] = frame["comparison_B"]
        frame["priority_line"] = (
            frame["comparison_A_name"].astype(str) + "__vs__" + frame["comparison_B_name"].astype(str)
        )
        frame["rank_score"] = (
            pd.to_numeric(frame.get("abs_delta_between_programs"), errors="coerce").fillna(0)
            + pd.to_numeric(frame.get("abs_dPSI_A"), errors="coerce").fillna(0)
            + pd.to_numeric(frame.get("abs_dPSI_B"), errors="coerce").fillna(0)
        )
        frame = frame.sort_values(["priority_line", "rank_score", "gene_symbol"], ascending=[True, False, True])
        return frame

    def split_shortlists(self, ranked: pd.DataFrame, top_n_dp: int = 12, top_n_ko: int = 12) -> tuple[pd.DataFrame, pd.DataFrame]:
        if ranked.empty:
            return pd.DataFrame(), pd.DataFrame()
        lines = [line for line in ranked["priority_line"].dropna().astype(str).drop_duplicates().tolist()]
        first = ranked[ranked["priority_line"] == lines[0]].copy() if len(lines) >= 1 else pd.DataFrame()
        second = ranked[ranked["priority_line"] == lines[1]].copy() if len(lines) >= 2 else pd.DataFrame()
        first = self._annotate_generic(first, rank_limit=top_n_dp)
        second = self._annotate_generic(second, rank_limit=top_n_ko)
        return first, second

    def build_all_shortlists(self, ranked: pd.DataFrame, rank_limit: int = 12) -> pd.DataFrame:
        if ranked.empty:
            return pd.DataFrame()
        lines = ranked["priority_line"].dropna().astype(str).drop_duplicates().tolist()
        frames = []
        for line in lines:
            subset = ranked[ranked["priority_line"] == line].copy()
            annotated = self._annotate_generic(subset, rank_limit=rank_limit)
            if not annotated.empty:
                frames.append(annotated)
        return pd.concat(frames, axis=0, ignore_index=True) if frames else pd.DataFrame()

    def build_cards_input(self, shortlist_dp: pd.DataFrame, shortlist_ko: pd.DataFrame) -> pd.DataFrame:
        frames = [frame for frame in (shortlist_dp, shortlist_ko) if frame is not None and not frame.empty]
        return pd.concat(frames, axis=0, ignore_index=True) if frames else pd.DataFrame()

    def _annotate_generic(self, frame: pd.DataFrame, rank_limit: int) -> pd.DataFrame:
        if frame.empty:
            return frame
        frame = frame.copy()
        comparison_a = str(frame["comparison_A_name"].iloc[0])
        comparison_b = str(frame["comparison_B_name"].iloc[0])
        group_name = self._safe_group_name(f"{comparison_a}_vs_{comparison_b}_shortlist")
        frame["biological_question"] = f"Compare splicing programs between {comparison_a} and {comparison_b}"
        frame["class_interpretation"] = frame["class"].map(
            {
                "A_only": f"Only_in_{self._safe_group_name(comparison_a)}",
                "B_only": f"Only_in_{self._safe_group_name(comparison_b)}",
                "shared_same_direction": "shared_same_direction",
                "opposite_direction": "opposite_direction",
            }
        ).fillna(frame["class"])
        frame["display_group"] = group_name
        frame["geneSymbol"] = frame.get("gene_symbol", frame.get("geneSymbol"))
        frame["GeneID"] = frame.get("gene_id", frame.get("GeneID"))
        frame["event_type"] = frame.get("event_type", frame.get("event_type_A"))
        frame["event_uid"] = frame.get("event_uid")
        frame["comparison_A"] = frame.get("comparison_A_name", frame.get("comparison_A"))
        frame["comparison_B"] = frame.get("comparison_B_name", frame.get("comparison_B"))
        frame = frame.drop_duplicates("geneSymbol").head(rank_limit).copy()
        return frame

    def _safe_group_name(self, value: str) -> str:
        return "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in value).strip("_")
