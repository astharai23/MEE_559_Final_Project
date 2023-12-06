import numpy as np
import h5py


hf = h5py.File("input_data.mat", "r")
data = hf['ans'][()]

data.astype(np.float64).tofile("data.bin")