from rootpy.io import root_open
from rootpy import asrootpy


directory = '~/msc/cern_stay/results/single_track_eff_golden_z_pt/'
f = root_open(directory + 'eff.root', 'read')
f_pt_z_only = root_open(directory + 'eff_z_pt_only.root', 'recreate')
mc = f.raw.tracks_mc
re = f.raw.tracks_re

# 3: p_T; 2: z_vtx
h = (asrootpy(re.Projection(3,2)) / asrootpy(mc.Projection(3,2)))  # y, x
h.name = 'single_track_eff'
f_pt_z_only.cd()
h.Write()
f.close()
f_pt_z_only.close()

