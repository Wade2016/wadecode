import pandas as pd
import sys
data=pd.read_csv(sys.argv[1])
print(data.mean())

