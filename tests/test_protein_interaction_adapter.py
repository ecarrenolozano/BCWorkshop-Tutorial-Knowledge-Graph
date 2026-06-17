"""Tests for the protein interaction BioCypher adapter."""

from pathlib import Path

import pandas as pd

from biocypher_tutorial_kg.adapters.protein_interaction_adapter import (
    ProteinInteractionAdapter,
)

REQUIRED_COLUMNS = [
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
]


def write_interaction_tsv(tmp_path: Path, rows: list[dict]) -> Path:
    data_path = tmp_path / "protein_interactions.tsv"
    pd.DataFrame(rows, columns=REQUIRED_COLUMNS).to_csv(
        data_path,
        sep="\t",
        index=False,
    )
    return data_path


def sample_interactions() -> list[dict]:
    return [
        {
            "source": "P00533",
            "source_genesymbol": "EGFR",
            "ncbi_tax_id_source": 9606,
            "entity_type_source": "protein",
            "target": "P04626",
            "target_genesymbol": "ERBB2",
            "ncbi_tax_id_target": 9606,
            "entity_type_target": "protein",
            "type": "activation",
            "is_directed": 1,
            "is_stimulation": 1,
            "is_inhibition": 0,
            "consensus_direction": 1,
            "consensus_stimulation": 1,
            "consensus_inhibition": 0,
        },
        {
            "source": "P00533",
            "source_genesymbol": "EGFR",
            "ncbi_tax_id_source": 9606,
            "entity_type_source": "protein",
            "target": "Q9Y243",
            "target_genesymbol": "AKT3",
            "ncbi_tax_id_target": 9606,
            "entity_type_target": "protein",
            "type": "inhibition",
            "is_directed": 1,
            "is_stimulation": 0,
            "is_inhibition": 1,
            "consensus_direction": 1,
            "consensus_stimulation": 0,
            "consensus_inhibition": 1,
        },
    ]


def test_adapter_initialization_preserves_source_and_config(tmp_path):
    data_path = tmp_path / "protein_interactions.tsv"
    config = {"batch_size": 500}

    adapter = ProteinInteractionAdapter(data_path, **config)

    assert adapter.data_source == data_path
    assert adapter.config == config


def test_metadata_describes_tsv_adapter(tmp_path):
    data_path = tmp_path / "protein_interactions.tsv"
    adapter = ProteinInteractionAdapter(data_path)

    metadata = adapter.get_metadata()

    assert metadata == {
        "name": "ProteinInteractionAdapter",
        "data_source": str(data_path),
        "data_type": "tsv",
        "version": "0.1.0",
        "adapter_class": "ProteinInteractionAdapter",
    }


def test_get_nodes_yields_unique_uniprot_protein_nodes(tmp_path):
    data_path = write_interaction_tsv(tmp_path, sample_interactions())
    adapter = ProteinInteractionAdapter(data_path)

    nodes = list(adapter.get_nodes())

    assert nodes == [
        (
            "P00533",
            "uniprot_protein",
            {
                "genesymbol": "EGFR",
                "ncbi_tax_id": "9606",
                "entity_type": "protein",
            },
        ),
        (
            "P04626",
            "uniprot_protein",
            {
                "genesymbol": "ERBB2",
                "ncbi_tax_id": "9606",
                "entity_type": "protein",
            },
        ),
        (
            "Q9Y243",
            "uniprot_protein",
            {
                "genesymbol": "AKT3",
                "ncbi_tax_id": "9606",
                "entity_type": "protein",
            },
        ),
    ]


def test_get_edges_yields_biocypher_edge_tuples(tmp_path):
    data_path = write_interaction_tsv(tmp_path, sample_interactions())
    adapter = ProteinInteractionAdapter(data_path)

    edges = list(adapter.get_edges())

    assert edges == [
        (
            "P00533_P04626_activation",
            "P00533",
            "P04626",
            "activation",
            {
                "is_directed": True,
                "is_stimulation": True,
                "is_inhibition": False,
                "consensus_direction": True,
                "consensus_stimulation": True,
                "consensus_inhibition": False,
            },
        ),
        (
            "P00533_Q9Y243_inhibition",
            "P00533",
            "Q9Y243",
            "inhibition",
            {
                "is_directed": True,
                "is_stimulation": False,
                "is_inhibition": True,
                "consensus_direction": True,
                "consensus_stimulation": False,
                "consensus_inhibition": True,
            },
        ),
    ]


def test_validate_data_source_requires_existing_tsv_with_interaction_columns(tmp_path):
    valid_path = write_interaction_tsv(tmp_path, sample_interactions())
    invalid_path = tmp_path / "invalid.tsv"
    pd.DataFrame([{"source": "P00533"}]).to_csv(
        invalid_path,
        sep="\t",
        index=False,
    )

    assert ProteinInteractionAdapter(valid_path).validate_data_source() is True
    assert ProteinInteractionAdapter(invalid_path).validate_data_source() is False
    assert (
        ProteinInteractionAdapter(tmp_path / "missing.tsv").validate_data_source()
        is False
    )
