import gdown
import hashlib
import os

import requests

PROTEOMICS_CORRELATION_URL = 'https://drive.google.com/uc?id=1_fIcLz6aUNHUtrQV6CpFo2ccToRlFLDw'
GSDC1_TTEST = 'https://drive.google.com/file/d/1bdFlgMZA0JfREIaBqrJTRd-zf3bdd6nH'
DEPMAP_TTEST = 'https://drive.google.com/file/d/1sAC4mm8xfAfGvSV-wCvB5yuyxaUjvJRt'
DEPMAP_GDS1_TTEST_TOP_HITS_MERGED = 'https://drive.google.com/file/d/1g3GOx_bu2nAIt9pj0RBXl9Uz-Xphf0Zn'
DEPMAP_GDSC1_TTEST_ALL_MERGED = 'https://drive.google.com/file/d/104VMExALrb3XjGkZK_aAKdtnzqAWmFNl'

DEPMAP_MUTATION_DATA = 'https://drive.google.com/file/d/1iedYFEZoDZxIrBZys8LXAmzD79DMIKsA'
# proteomics data taken from https://depmap.sanger.ac.uk/documentation/datasets/proteomics/
PROTEOMICS_MATRIX_AVERAGE_ZSCORE = 'https://drive.google.com/file/d/10Vonswa_Cp5zAWsqvs7k0rNLVbM100v7'
GDSC1_DOSE_RESPONSE = 'https://drive.google.com/file/d/1XyJyM2ToI1U1I29ZwTYL2WSkNfmg275f'
DEPMAP_MODEL_MAPPING='https://drive.google.com/file/d/1oAcdxGsHXq9KZ4U7I1SeYOflLkhxWWKE'

PRISM_DRUG_REPURPOSING = 'https://drive.google.com/file/d/1DcvKjdJlKE6zy1vHdw-8iHPn3cxF_QvY'

DEPMAP_EXPRESSION = "https://drive.google.com/file/d/1c-mgX0K_-OOK3MmR93G3og1BLHTOkxQa"

def download_google_drive_url(url):
    """
    Downloads a CSV file from Google Drive only if it doesn't exist locally.
    
    Args:
        url (str): Google Drive URL of the CSV file
        
    Returns:
        str: Full path to the downloaded or existing file
    """
    # Create data directory if it doesn't exist
    data_dir = os.path.expanduser("~/workspace/bioinformatics_analysis/data")
    os.makedirs(data_dir, exist_ok=True)
    
    # Create hash of URL for filename
    url_hash = hashlib.md5(url.encode()).hexdigest()
    output_path = os.path.join(data_dir, f"{url_hash}.csv")
    
    # Convert file/d/ URLs to uc?id= format
    if '/file/d/' in url:
        file_id = url.split('/file/d/')[1].split('/')[0]
        url = f'https://drive.google.com/uc?id={file_id}'
    
    # Download only if file doesn't exist
    if not os.path.exists(output_path):
        print(f"Downloading file from {url}")
        try:
            gdown.download(url, output_path, quiet=False, fuzzy=True)
            # Verify if downloaded file is CSV
            with open(output_path, 'r', encoding='utf-8') as f:
                first_line = f.readline()
                if '<html' in first_line.lower():
                    os.remove(output_path)
                    raise ValueError("Downloaded file appears to be HTML instead of CSV")
        except Exception as e:
            print(f"Error downloading file: {e}")
            if os.path.exists(output_path):
                os.remove(output_path)
            return None
    else:
        print(f"File already exists at {output_path}")
    
    return output_path


def download_url(url):
    # Define the local filename to save the CSV
    data_dir = os.path.expanduser("~/workspace/bioinformatics_analysis/data")
    os.makedirs(data_dir, exist_ok=True)

    # Create hash of URL for filename
    url_hash = hashlib.md5(url.encode()).hexdigest()
    output_path = os.path.join(data_dir, f"{url_hash}.csv")

    try:
        # Send a GET request to the URL
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)

        # Open the local file in binary write mode
        with open(output_path, 'wb') as f:
            # Iterate over the response content in chunks
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        print(f"File '{output_path}' downloaded successfully.")

    except requests.exceptions.RequestException as e:
        print(f"Error downloading the file: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return output_path
