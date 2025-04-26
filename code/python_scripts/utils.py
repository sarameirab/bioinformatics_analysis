import requests


def map_uniprot_to_gene(protein_ids):
    """
    Map UniProt protein IDs to gene names using UniProt's mapping service

    Args:
        protein_ids (list): List of UniProt protein IDs

    Returns:
        dict: Mapping of protein IDs to gene names
    """
    # Join protein IDs with spaces
    ids = ' '.join(protein_ids)

    # Configure the UniProt API request
    url = 'https://rest.uniprot.org/idmapping/run'
    params = {
        'from': 'UniProtKB_AC-ID',
        'to': 'Gene_Name',
        'ids': ids
    }

    try:
        # Submit the mapping job
        response = requests.post(url, data=params)
        response.raise_for_status()  # Raise exception for bad status codes
        job_id = response.json()['jobId']

        # Get the results
        result_url = f'https://rest.uniprot.org/idmapping/stream/{job_id}'
        results = requests.get(result_url)
        results.raise_for_status()

        # Parse the results - updated format
        mapping = {}
        result_data = results.json()

        if 'results' in result_data:
            for entry in result_data['results']:
                mapping[entry['from']] = entry.get('to', '')
        else:
            # Handle new API format
            for entry in result_data:
                from_id = entry.get('from', '')
                to_name = entry.get('to', {}).get('geneName', '')
                if from_id and to_name:
                    mapping[from_id] = to_name

        return mapping

    except requests.exceptions.RequestException as e:
        print(f"Error during API request: {e}")
        return {}
    except KeyError as e:
        print(f"Error parsing API response: {e}")
        print("API Response:", results.text[:200])  # Print first 200 chars of response
        return {}
    except Exception as e:
        print(f"Unexpected error: {e}")
        return {}