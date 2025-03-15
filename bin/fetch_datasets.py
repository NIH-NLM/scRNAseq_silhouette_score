import os
import sys
import json
import requests

# Constants
DATASETS_API_URL = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
SMALLEST_DATASET_ID = "0895c838-e550-48a3-a777-dbcd35d30272"

def fetch_datasets(collections_json, test_mode):
    """Fetch datasets from collections_info.json or API (test mode)."""

    if test_mode:
        # Fetch only the smallest dataset for testing
        print(f"Fetching test dataset: {SMALLEST_DATASET_ID}")
        url = f"{DATASETS_API_URL}/{SMALLEST_DATASET_ID}/versions"

        try:
            response = requests.get(url)
            response.raise_for_status()
            datasets = response.json()
        except requests.exceptions.RequestException as e:
            print(f"ERROR: Failed to fetch dataset {SMALLEST_DATASET_ID} - {e}")
            return

        datasets_output_filename = f"datasets_{SMALLEST_DATASET_ID}.json"
        with open(datasets_output_filename, "w") as f:
            json.dump(datasets, f, indent=2)
        print(f"Saved test dataset: {datasets_output_filename}")

    else:
        # Fetch all datasets from the provided JSON file
        try:
            with open(collections_json, "r") as f:
                collections_data = json.load(f)
        except (json.JSONDecodeError, FileNotFoundError) as e:
            print(f"ERROR: Unable to read JSON file {collections_json} - {e}")
            return

        if "datasets" not in collections_data:
            print("ERROR: JSON file does not contain a 'datasets' field.")
            return

        datasets = collections_data["datasets"]

        for dataset in datasets:
            dataset_id = dataset.get("dataset_id")

            if not dataset_id:
                print("WARNING: Skipping dataset without dataset_id")
                continue

            dataset_filename = f"datasets_{dataset_id}.json"
            with open(dataset_filename, "w") as f:
                json.dump(dataset, f, indent=2)
            print(f"Saved dataset: {dataset_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("ERROR: Usage: python fetch_datasets.py <collections_json> <test_mode>")
        sys.exit(1)

    collections_json = sys.argv[1]
    test_mode = sys.argv[2].lower() == "true"

    fetch_datasets(collections_json, test_mode)

