import struct
import numpy as np

if __name__== "__main__":
    nx = 20
    ny = 20
    fomTotDofs = nx * ny * 2

    D = np.fromfile("burgers_outflow2d_solution.bin")
    nt = int(np.size(D) / fomTotDofs)
    D = np.reshape(D, (nt, fomTotDofs))
    D = D[-1, :]
    D = np.reshape(D, (nx * ny, 2))
    u = D[:,0]
    np.savetxt("u.txt", u)

    goldD = np.loadtxt("u_gold.txt")
    assert np.allclose(u.shape, goldD.shape)
    assert np.isnan(u).all() == False
    assert np.isnan(goldD).all() == False
    assert np.allclose(u, goldD, rtol=1e-10, atol=1e-12)

    # check runtime file
    f = open('runtime.bin', 'rb')
    contents = f.read()
    assert len(contents) == 16
    niters = struct.unpack('Q', contents[:8])[0]
    assert niters == 1
    timeval = struct.unpack('d', contents[8:16])[0]
    assert timeval > 0.0
