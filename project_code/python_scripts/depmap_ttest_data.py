import pandas as pd
from urls import *
from utils import *

path = download_google_drive_url(DEPMAP_GDS1_TTEST_TOP_HITS_MERGED)
df = pd.read_csv(path)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', None)
print(df.head())
print(f"DataFrame dimensions (rows, columns): {df.shape}")