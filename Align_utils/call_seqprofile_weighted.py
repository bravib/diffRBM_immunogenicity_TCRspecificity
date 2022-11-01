rootf= path_to_folder + '/Example_diffRBM/Align_utils'

import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(path_c)
eng.call_seqprofile_weighted(21)
eng.quit()
