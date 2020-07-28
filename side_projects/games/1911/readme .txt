pickle and than zlib is the best :)
didn't really tried hdf5

z.write(zlib.compress(pickle.dumps(chr_dict, pickle.HIGHEST_PROTOCOL),9))