import itertools as it
import pandas as pd
params = {
          'pair': ('y', 'n'),
          'cnv': ('y','n'),
          'ffpe': ('y')
          }

df = pd.DataFrame([row for row in it.product(*params.values())], 
                   columns = params.keys())
df.query('pair == "y" | cnv == "n"').to_csv('config/params.csv', index=False)
