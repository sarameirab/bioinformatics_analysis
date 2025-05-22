import os
import shutil
import gdown


def download_gdrive_folder(folder_url, destination_path):
    """
    Downloads a Google Drive folder to the specified destination path.

    Args:
        folder_url (str): URL of the Google Drive folder
        destination_path (str): Local path where the folder should be downloaded
    """
    # Expand the ~ to full path
    destination_path = os.path.expanduser(destination_path)

    # Remove the directory if it already exists
    if os.path.exists(destination_path):
        print(f"Removing existing directory: {destination_path}")
        shutil.rmtree(destination_path)

    # Create the parent directory if it doesn't exist
    os.makedirs(os.path.dirname(destination_path), exist_ok=True)

    # Download the folder
    print(f"Downloading Google Drive folder to: {destination_path}")
    gdown.download_folder(url=folder_url,
                          output=destination_path,
                          quiet=False)

    print("Download complete!")


if __name__ == "__main__":
    # Your Google Drive folder URL
    folder_url = "https://drive.google.com/drive/folders/1fR-ecJqCcfijP6LaiYw2_4QfRYEIs7-K?usp=drive_link"

    # Destination path
    dest_path = "~/workspace/ai_agent_lab_data"

    download_gdrive_folder(folder_url, dest_path)