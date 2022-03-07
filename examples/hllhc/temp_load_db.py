import io
import numpy as np
import pandas as pd




def temp_load_colldb(filename):
   with open(filename, "r") as infile:
       coll_data_string = ""
       family_settings = {}
       family_types = {}
       onesided = {}

       for l_no, line in enumerate(infile):
           if line.startswith("#"):
               continue # Comment

           sline = line.split()
           if len(sline) < 6:
               if sline[0].lower() == "nsig_fam":
                   family_settings[sline[1]] = float(sline[2])
                   family_types[sline[1]] = sline[3]
               elif sline[0].lower() == "onesided":
                   onesided[sline[1]] = int(sline[2])
               elif sline[0].lower() == "settings":
                   pass # Acknowledge and ignore this line
               else:
                   print(f"Unknown setting {line}")
           else:
               coll_data_string += line

   names = ["name", "opening", "material", "length", "angle", "offset"]

   df = pd.read_csv(io.StringIO(coll_data_string), delim_whitespace=True,
                    index_col=False, names=names)

   df["angle"] = np.deg2rad(df["angle"]) # convert to radians
   df["name"] = df["name"].str.lower() # Make the names lowercase for easy processing
   df["nsigma"] = df["opening"].apply(lambda s: family_settings.get(s, s))
   df["type"] = df["opening"].apply(lambda s: family_types.get(s, "UNKNOWN"))
   df["side"] = df["name"].apply(lambda s: onesided.get(s, 0))
   df = df.set_index("name").T

   return df.to_dict()

