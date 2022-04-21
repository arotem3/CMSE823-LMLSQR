import numpy as np
import matplotlib.pyplot as plt

times_chol = np.loadtxt('timing_chol.txt')
times_qr = np.loadtxt('timing_qr.txt')
times_svd = np.loadtxt('timing_svd.txt')

ms = np.loadtxt('ms.txt', dtype=int)
ratios = np.loadtxt('ratios.txt')

for (m, tchol, tqr, tsvd) in zip(ms, times_chol, times_qr, times_svd):
    plt.clf()
    plt.semilogy(ratios, tchol, 'o-', label='Cholesky')
    plt.semilogy(ratios, tsvd, 'd-', label='SVD')
    plt.semilogy(ratios, tqr, 's-', label='QR')
    plt.xlabel('n/m')
    plt.ylabel('runtime (s)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'm{m}.png')

