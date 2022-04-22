import numpy as np
import matplotlib.pyplot as plt

ms = np.loadtxt('ms.txt', dtype=int)
ratios = np.loadtxt('ratios.txt')

for k in [2, 5, 10]:
    times_chol = np.loadtxt(f'timing_chol_{k}.txt')
    times_qr = np.loadtxt(f'timing_qr_{k}.txt')
    times_svd = np.loadtxt(f'timing_svd_{k}.txt')
    times_qrf = np.loadtxt(f'timing_qrf_{k}.txt')

    for (m, tchol, tqr, tqrf, tsvd) in zip(ms, times_chol, times_qr, times_qrf, times_svd):
        plt.clf()
        plt.semilogy(ratios, tchol, 'o-', label='Cholesky')
        plt.semilogy(ratios, tsvd, 'd-', label='SVD')
        plt.semilogy(ratios, tqrf, '*-', label='naive QR')
        plt.semilogy(ratios, tqr, 's-', label='updated QR')
        plt.xlabel('n/m')
        plt.ylabel('runtime (s)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'm{m}_k{k}.png')

tchol = np.loadtxt('lmsn_chol.txt')
tsvd = np.loadtxt('lmsn_svd.txt')
tqr = np.loadtxt('lmsn_qr.txt')
tqrf = np.loadtxt('lmsn_qrf.txt')
ns = np.loadtxt('lmsn_n.txt', dtype=int)

plt.clf()
plt.semilogy(ns, tchol, 'o-', label='Cholesky')
plt.semilogy(ns, tsvd, 'd-', label='SVD')
plt.semilogy(ns, tqrf, '*-', label='naive QR')
plt.semilogy(ns, tqr, 's-', label='updated QR')
plt.xlabel('n')
plt.ylabel('runtime (s)')
plt.legend()
plt.tight_layout()
plt.savefig("lmsn.png")