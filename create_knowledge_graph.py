#!/usr/bin/env python3
"""
BioCypherWorkshop Tutorial Knowledge Graph - A beginner-friendly BioCypher project for learning how to build and extend a small knowledge graph from multiple data sources.

This script creates a knowledge graph using BioCypher and the ProteinInteractionAdapter.
"""

import logging
from pathlib import Path

from biocypher import BioCypher
from biocypher_tutorial_kg.adapters.protein_interaction_adapter import ProteinInteractionAdapter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Main function to create the knowledge graph."""
    logger.info("Starting BioCypherWorkshop Tutorial Knowledge Graph knowledge graph creation")
    
    # Initialize BioCypher
    bc = BioCypher(
        biocypher_config_path="config/biocypher_config.yaml",
        schema_config_path="config/schema_config.yaml"
    )
    
    # Initialize the adapter
    # TODO: Configure your CSV data source path here
    data_source = "data/your_data.csv"  # Update this with your actual CSV file path
    
    adapter = ProteinInteractionAdapter(
        data_source=data_source,
        # Add any additional configuration parameters here
    )
    
    # Create the knowledge graph
    logger.info("Creating knowledge graph...")
    bc.write_nodes(adapter.get_nodes())
    bc.write_edges(adapter.get_edges())
    
    logger.info("Knowledge graph creation completed successfully!")

    # Create final summary
    bc.summary()


if __name__ == "__main__":
    main()
