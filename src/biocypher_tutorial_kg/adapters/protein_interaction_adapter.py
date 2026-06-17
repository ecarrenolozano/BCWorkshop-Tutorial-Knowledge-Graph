"""
ProteinInteractionAdapter Adapter

This adapter handles protein interaction TSV data for BioCypher.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = {
    "source",
    "source_genesymbol",
    "ncbi_tax_id_source",
    "entity_type_source",
    "target",
    "target_genesymbol",
    "ncbi_tax_id_target",
    "entity_type_target",
    "type",
    "is_directed",
    "is_stimulation",
    "is_inhibition",
    "consensus_direction",
    "consensus_stimulation",
    "consensus_inhibition",
}

BOOLEAN_COLUMNS = {
    "is_directed",
    "is_stimulation",
    "is_inhibition",
    "consensus_direction",
    "consensus_stimulation",
    "consensus_inhibition",
}


class ProteinInteractionAdapter:
    """
    Adapter for protein interaction TSV data source.

    This adapter implements the BioCypher adapter interface and is intentionally
    scaffolded as a workshop exercise.
    """

    def __init__(self, data_source: str | Path, **kwargs):
        """
        Initialize the adapter.

        Args:
            data_source: Path to the TSV data source
            **kwargs: Additional configuration parameters
        """
        self.data_source = data_source
        self.config = kwargs
        logger.info(
            "Initialized ProteinInteractionAdapter with data source: %s", data_source
        )

    def get_nodes(self):
        """
        Extract nodes from the data source.

        Yields:
            Tuples of (node_id, node_label, properties) for each node.
        """
        logger.info("Extracting nodes from data source")

        # BioCypher node tuple format:
        # (node_id, node_label, properties)

        # ── REFERENCE ────────────────────────────────────────────────────────
        # Steps 1–2.4 are already implemented. Read them to understand the
        # data preparation before you implement the yield step (TODO 4).
        # ─────────────────────────────────────────────────────────────────────

        # Step 1: Read the interaction TSV into a DataFrame.
        df = pd.read_csv(self.data_source, sep="\t")
        self._validate_columns(df)

        # Step 2.1: Build the source-side Protein table.
        source_nodes = df[
            ["source", "source_genesymbol", "ncbi_tax_id_source", "entity_type_source"]
        ].copy()
        source_nodes = source_nodes.rename(
            columns={
                "source": "node_id",
                "source_genesymbol": "genesymbol",
                "ncbi_tax_id_source": "ncbi_tax_id",
                "entity_type_source": "entity_type",
            }
        )

        # Step 2.2: Build the target-side Protein table.
        target_nodes = df[
            ["target", "target_genesymbol", "ncbi_tax_id_target", "entity_type_target"]
        ].copy()
        target_nodes = target_nodes.rename(
            columns={
                "target": "node_id",
                "target_genesymbol": "genesymbol",
                "ncbi_tax_id_target": "ncbi_tax_id",
                "entity_type_target": "entity_type",
            }
        )

        # Step 2.3: Merge source-side and target-side proteins.
        proteins_df = pd.concat([source_nodes, target_nodes], ignore_index=True)

        # Step 2.4: Remove duplicate proteins by node_id.
        proteins_df = proteins_df.drop_duplicates(subset=["node_id"])

        # ── TODO 4 ───────────────────────────────────────────────────────────
        # Task: Yield one BioCypher node tuple per unique protein.
        # Why:  BioCypher expects this exact 3-element tuple shape:
        #       (node_id, node_label, properties)
        # Hint: Iterate over proteins_df rows. For each row yield:
        #         (str(row["node_id"]),
        #          "uniprot_protein",
        #          {"genesymbol": ..., "ncbi_tax_id": str(...), "entity_type": ...})
        # ─────────────────────────────────────────────────────────────────────
        # ---------- YOUR CODE STARTS HERE ----------
        for _, row in proteins_df.iterrows():
            yield (
                str(row["node_id"]),
                "uniprot_protein",
                {
                    "genesymbol": str(row["genesymbol"]),
                    "ncbi_tax_id": str(row["ncbi_tax_id"]),
                    "entity_type": str(row["entity_type"]),
                },
            )
        # ---------- YOUR CODE ENDS HERE ----------

    def get_edges(self):
        """
        Extract edges from the data source.

        Yields:
            Tuples of (edge_id, source_id, target_id, edge_label, properties)
            for each edge.
        """
        logger.info("Extracting edges from data source")

        # BioCypher edge tuple format:
        # (edge_id, source_id, target_id, edge_label, properties)

        # ── TODO 5a ──────────────────────────────────────────────────────────
        # Task: Read the interaction TSV into a DataFrame.
        # Why:  Each row represents one interaction and maps to one edge tuple.
        # Hint: Use pd.read_csv(self.data_source, sep="\t").
        # ─────────────────────────────────────────────────────────────────────
        # ---------- YOUR CODE STARTS HERE ----------
        df = pd.read_csv(self.data_source, sep="\t")
        # ---------- YOUR CODE ENDS HERE ----------

        # ── TODO 5e ──────────────────────────────────────────────────────────
        # Task: Yield one BioCypher edge tuple per interaction row.
        # Why:  BioCypher expects this exact 5-element tuple shape:
        #       (edge_id, source_id, target_id, edge_label, properties)
        # Hint: Use interaction_type as the edge_label. Optionally skip rows
        #       where source, target, or type is missing.
        # ─────────────────────────────────────────────────────────────────────
        # ---------- YOUR CODE STARTS HERE ----------
        for _, row in df.iterrows():
            source_id = str(row["source"])
            target_id = str(row["target"])
            interaction_type = str(row["type"])
            edge_id = f"{source_id}_{target_id}_{interaction_type}"
            properties = {
                column: self._parse_bool(row[column]) for column in BOOLEAN_COLUMNS
            }
            yield (edge_id, source_id, target_id, interaction_type, properties)

        # ---------- YOUR CODE ENDS HERE ----------

    def get_metadata(self) -> dict[str, Any]:
        """
        Get metadata about the data source.

        Returns:
            Dictionary containing metadata
        """
        return {
            "name": "ProteinInteractionAdapter",
            "data_source": str(self.data_source),
            "data_type": "tsv",
            "version": "0.1.0",
            "adapter_class": "ProteinInteractionAdapter",
        }

    def validate_data_source(self) -> bool:
        """
        Validate that the TSV data source is accessible and properly formatted.

        Returns:
            True if data source is valid, False otherwise
        """
        try:
            data_path = Path(self.data_source)
            if not data_path.exists() or not data_path.is_file():
                return False

            # Read one row to validate TSV format and required columns.
            df = pd.read_csv(data_path, sep="\t", nrows=1)
            self._validate_columns(df)
            return True

        except (OSError, ValueError, pd.errors.ParserError) as exc:
            logger.exception("Data source validation failed: %s", exc)
            return False

    @staticmethod
    def _validate_columns(df: pd.DataFrame) -> None:
        missing_columns = REQUIRED_COLUMNS.difference(df.columns)
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"Missing required TSV columns: {missing}")

    @staticmethod
    def _parse_bool(value: Any) -> bool:
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)):
            return bool(value)
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "t", "yes", "y"}
        return bool(value)
