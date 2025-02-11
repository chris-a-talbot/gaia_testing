#!/usr/bin/env python3

import sys
import tskit
import pandas as pd
import json
from pathlib import Path
import logging
from dataclasses import dataclass
from typing import Dict, Set, List, Optional
import numpy as np
from numpy import integer as np_integer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('validation_report.log')
    ]
)
logger = logging.getLogger(__name__)


@dataclass
class ValidationError:
    error_type: str
    description: str
    subset_name: Optional[str] = None
    details: Optional[Dict] = None


class DataValidator:
    def __init__(self, tree_name: str):
        """Initialize validator with paths and load base data"""
        self.tree_name = tree_name
        self.errors: List[ValidationError] = []

        # Define paths
        self.base_tree_path = Path(f"./trees/{tree_name}.trees")
        self.subsets_dir = Path(f"./trees/subsets/{tree_name}")
        self.base_locations_path = Path(f"./sample_locations/{tree_name}_locations.csv")
        self.subset_locations_dir = Path(f"./sample_locations/subsets/{tree_name}")

        # Load base data
        self.load_base_data()

    def load_base_data(self):
        """Load the base tree sequence and locations data"""
        try:
            self.base_ts = tskit.load(str(self.base_tree_path))
            logger.info(f"Loaded base tree sequence: {self.base_tree_path}")

            self.base_locations = pd.read_csv(self.base_locations_path)
            logger.info(f"Loaded base locations data: {self.base_locations_path}")

            # Get base sample nodes
            self.base_sample_nodes = set(self.base_ts.samples())
            logger.info(f"Found {len(self.base_sample_nodes)} sample nodes in base tree")

        except Exception as e:
            raise RuntimeError(f"Failed to load base data: {str(e)}")

    def validate_base_data_integrity(self):
        """Validate the integrity of the base data"""
        # Check base locations data columns
        required_cols = {'node_id', 'x', 'y', 'z'}
        missing_cols = required_cols - set(self.base_locations.columns)
        if missing_cols:
            self.errors.append(ValidationError(
                "missing_columns",
                f"Base locations file missing required columns: {missing_cols}"
            ))

        # Get sample nodes from location data
        location_node_ids = set(self.base_locations['node_id'])

        # Check for missing samples
        missing_samples = self.base_sample_nodes - location_node_ids
        if missing_samples:
            self.errors.append(ValidationError(
                "missing_sample_locations",
                f"Found {len(missing_samples)} sample nodes in tree without location data",
                details={'missing_nodes': sorted(list(missing_samples))}
            ))

        # Check for extra samples
        extra_samples = location_node_ids - self.base_sample_nodes
        if extra_samples:
            self.errors.append(ValidationError(
                "extra_sample_locations",
                f"Found {len(extra_samples)} location entries for non-sample nodes",
                details={'extra_nodes': sorted(list(extra_samples))}
            ))

    def is_ancient(self, time_value: float) -> bool:
        """Determine if a sample is ancient based on its time value"""
        return time_value not in [0.0, 1.0]

    def validate_subset(self, subset_path: Path):
        """Validate a single subset tree sequence and its corresponding location data"""
        subset_name = subset_path.stem
        logger.info(f"Validating subset: {subset_name}")

        try:
            # Load subset tree and metadata
            subset_ts = tskit.load(str(subset_path))
            metadata = self.load_subset_metadata(subset_path)

            # Get sample nodes from subset
            subset_sample_nodes = set(subset_ts.samples())

            # Load corresponding location data
            location_path = self.subset_locations_dir / f"{subset_name}_locations.csv"
            if not location_path.exists():
                self.errors.append(ValidationError(
                    "missing_location_file",
                    f"No location file found for subset {subset_name}",
                    subset_name
                ))
                return

            subset_locations = pd.read_csv(location_path)

            # Check required columns
            required_cols = {'node_id', 'x', 'y', 'z'}
            missing_cols = required_cols - set(subset_locations.columns)
            if missing_cols:
                self.errors.append(ValidationError(
                    "missing_columns",
                    f"Subset locations file missing required columns: {missing_cols}",
                    subset_name
                ))
                return

            # Validate sampling scheme compliance
            self.validate_sampling_scheme(subset_ts, subset_locations, metadata, subset_name)

            # Validate data correspondence
            self.validate_data_correspondence(subset_sample_nodes, subset_locations, subset_name)

        except Exception as e:
            logger.error(f"Error processing subset {subset_name}: {str(e)}", exc_info=True)
            self.errors.append(ValidationError(
                "processing_error",
                f"Error processing subset {subset_name}: {str(e)}",
                subset_name
            ))

    def validate_sampling_scheme(self, subset_ts: tskit.TreeSequence,
                                 subset_locations: pd.DataFrame,
                                 metadata: Dict,
                                 subset_name: str):
        """Validate that the subset matches its intended sampling scheme"""

        # Check number of samples
        if 'n_nodes' in metadata:
            expected_samples = metadata['n_nodes']
            actual_samples = len(subset_ts.samples())
            if actual_samples != expected_samples:
                self.errors.append(ValidationError(
                    "sample_count_mismatch",
                    f"Expected {expected_samples} sample nodes, found {actual_samples}",
                    subset_name
                ))

    def validate_data_correspondence(self, subset_sample_nodes: Set[int],
                                     subset_locations: pd.DataFrame,
                                     subset_name: str):
        """Validate correspondence between tree sequence samples and location data"""
        # Get node IDs from locations
        location_node_ids = set(subset_locations['node_id'])

        # Validate node counts match
        if len(subset_sample_nodes) != len(location_node_ids):
            self.errors.append(ValidationError(
                "node_count_mismatch",
                f"Tree has {len(subset_sample_nodes)} sample nodes but locations file has {len(location_node_ids)} entries",
                subset_name
            ))

    def load_subset_metadata(self, subset_path: Path) -> Dict:
        """Load metadata for a subset tree sequence"""
        meta_path = subset_path.parent / f"{subset_path.stem}_meta.json"
        try:
            with open(meta_path) as f:
                return json.load(f)
        except Exception as e:
            self.errors.append(ValidationError(
                "metadata_error",
                f"Failed to load metadata for {subset_path.stem}: {str(e)}",
                subset_path.stem
            ))
            return {}

    def validate_all(self):
        """Run all validations and generate report"""
        logger.info("Starting validation process")

        # Validate base data
        self.validate_base_data_integrity()

        # Validate each subset
        subset_files = list(self.subsets_dir.glob("*.trees"))
        logger.info(f"Found {len(subset_files)} subset files to validate")

        for subset_file in subset_files:
            self.validate_subset(subset_file)

        # Generate validation report
        self.generate_report()

    def generate_report(self):
        """Generate a detailed validation report"""
        report_path = Path(f"validation_report_{self.tree_name}.json")

        report = {
            'tree_name': self.tree_name,
            'timestamp': pd.Timestamp.now().isoformat(),
            'total_subsets_checked': len(list(self.subsets_dir.glob("*.trees"))),
            'total_errors': len(self.errors),
            'errors_by_type': {},
            'error_details': []
        }

        # Group errors by type
        for error in self.errors:
            if error.error_type not in report['errors_by_type']:
                report['errors_by_type'][error.error_type] = 0
            report['errors_by_type'][error.error_type] += 1

            error_detail = {
                'type': error.error_type,
                'description': error.description,
                'subset': error.subset_name
            }
            if error.details:
                # Convert numpy types to native Python types
                details = {}
                for key, value in error.details.items():
                    if isinstance(value, list):
                        details[key] = [int(x) if isinstance(x, np_integer) else x for x in value]
                    else:
                        details[key] = int(value) if isinstance(value, np_integer) else value
                error_detail['details'] = details

            report['error_details'].append(error_detail)

        # Save report
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)

        logger.info(f"Validation complete. Report saved to {report_path}")

        # Log summary
        if self.errors:
            logger.error(f"Found {len(self.errors)} validation errors")
            for error_type, count in report['errors_by_type'].items():
                logger.error(f"  {error_type}: {count} errors")
        else:
            logger.info("No validation errors found")


def main():
    if len(sys.argv) != 2:
        print("Usage: python validate_data.py <tree_name>")
        print("Example: python validate_data.py tree-S0.3-R1")
        sys.exit(1)

    tree_name = sys.argv[1]

    try:
        validator = DataValidator(tree_name)
        validator.validate_all()
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()