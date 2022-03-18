import io

import numpy as np
import pandas as pd
from .beam_elements import Collimator

_parameters_to_be_extracted_from_twiss = (
            'x y px py betx bety alfx alfy gamx gamy dx dpx dy dpy mux muy'.split())
_locations = ['at_center_active_part', 'at_start_active_part', 'at_start_element']


class CollimatorManager:
    def __init__(self, coll_db_txt_file, nemitt_ref_x, nemitt_ref_y):

        colldb = _load_colldb(coll_db_txt_file)

        # Assumes marker in the the line is at the center of the active length
        inputcolldf = pd.DataFrame.from_dict(colldb).transpose()\
                            .rename(columns={'length': 'active_length'})
        inputcolldf['inactive_length_at_start'] = 0 # To be used for fluka
        inputcolldf['inactive_length_at_end'] = 0   # To be used for fluka

        # Prepare data structure for optics and orbit (initialized with NaN for now)
        temp_dfs = []
        for ll in _locations:
            cols = pd.MultiIndex.from_tuples(
            [(ll, nn) for nn in ['s'] + _parameters_to_be_extracted_from_twiss])
            newdf = pd.DataFrame(index=inputcolldf.index, columns=cols)
            temp_dfs.append(newdf)

        colldf = temp_dfs[0]
        for dd in temp_dfs[1:]:
            colldf = colldf.join(dd)

        # Get collimator info from input
        for cc in inputcolldf.columns[::-1]: #Try to get a nice ordering...
            colldf.insert(0, column=cc, value=inputcolldf[cc])

        assert np.all(colldf['offset']==0), "Halfgap offset not implemented"

        colldf['angle_rad'] = colldf['angle']
        colldf['angle_deg'] = colldf['angle']*180/np.pi
        colldf.drop('angle', axis=1, inplace=True, level=0)

        colldf['length'] = (colldf['active_length']
                            + colldf['inactive_length_at_start']
                            + colldf['inactive_length_at_end'])
        colldf['name'] = colldf.index

        self.table = colldf
        self.line = None
        self.nemitt_ref_x = nemitt_ref_x
        self.nemitt_ref_y = nemitt_ref_y

    @property
    def collimator_names(self):
        return list(self.table.index.values)

    def install_collimators_in_line(self, line, marker_positions='center_active'):

        assert self.line is None, (
            "A line is already associated to this collimator manager")

        if marker_positions != 'center_active':
            raise NotImplementedError
            
        for kk in self.collimator_names:
            assert kk in line.element_names

        self.line = line

        colldf=self.table

        colldf['at_center_active_part', 's'] = line.get_s_position(colldf['name'].to_list())
        colldf['at_start_active_part', 's'] = (
            colldf['at_center_active_part', 's'] - colldf['active_length']/2)
        colldf['at_start_element', 's'] = (
            colldf['at_start_active_part', 's'] - colldf['inactive_length_at_start'])

        for nn in colldf.index.values:
            print(f'Install {nn}')
            newcoll = Collimator(
                    inactive_length_at_start=colldf.loc[nn, 'inactive_length_at_start'],
                    inactive_length_at_end=colldf.loc[nn, 'inactive_length_at_end'],
                    active_length=colldf.loc[nn, 'active_length'],
                    n_slices=1, angle=colldf.loc[nn, 'angle_deg'].values[0],
                    a_min=-1, a_max=1,
                    b_min=-1, b_max=1
                    )
            s_insert = colldf['at_start_element', 's'][nn]
            line.insert_element(element=newcoll, name=nn, at_s=s_insert)

    def compute_optics_at_collimators(self):

        assert self.line.tracker is not None, (
            "Please build tracker before calling this method"
        )

        tracker = self.line.tracker
        colldf = self.table

        s_twiss = []
        for ll in _locations:
            s_twiss.extend(colldf[ll]['s'].to_list())

        n_coll = len(colldf)
        assert len(s_twiss) == len(_locations) * n_coll

        tw = tracker.twiss(at_s=s_twiss)

        # Extract optics info
        for ill, ll in enumerate(_locations):
            for nn in _parameters_to_be_extracted_from_twiss:
                colldf[ll, nn] = tw[nn][ill*n_coll:(ill+1)*n_coll]

        beta0_gamma0 = (tracker.particle_ref._xobject.beta0[0]
                * tracker.particle_ref._xobject.gamma0[0])
        for ll in _locations:
            colldf[ll, 'sigmax'] = np.sqrt(colldf[ll, 'betx']*self.nemitt_ref_x/beta0_gamma0)
            colldf[ll, 'sigmay'] = np.sqrt(colldf[ll, 'bety']*self.nemitt_ref_y/beta0_gamma0)

    def set_collimator_openings(self, collimator_names='all'):

        if collimator_names != 'all':
            raise NotImplementedError

        assert self.line.tracker is not None, (
            "Please build tracker before calling this method"
        )

        self.compute_optics_at_collimators()

        colldf = self.table

        # Compute halfgap
        colldf['halfgap_m'] = colldf['nsigma'].values * np.sqrt(
            (colldf['at_center_active_part', 'sigmax']*np.cos(np.float_(colldf['angle_rad'].values)))**2
            + (colldf['at_center_active_part', 'sigmay']*np.sin(np.float_(colldf['angle_rad'].values)))**2)

        line = self.line

        # Configure collimators
        for nn in self.collimator_names:
            line[nn].dx = colldf['at_center_active_part', 'x'][nn]
            line[nn].dy = colldf['at_center_active_part', 'y'][nn]
            line[nn].a_max = colldf['halfgap_m'][nn]
            line[nn].a_min = -colldf['halfgap_m'][nn]


def _load_colldb(filename):
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
