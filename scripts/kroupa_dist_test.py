import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
# masses = np.linspace(5,1000,150)
# masses_sec = np.linspace(5,1000,150)
# N = 100000

# kroupa_imf = masses**(-2.3)
# kroupa_imf_sec = []
# for m in masses_sec:
#     if m <= 0.08 :
#         kroupa_imf_sec.append(m**(-0.3))
#     if 0.08 < m < 0.5:
#         kroupa_imf_sec.append(m**(-1.3))
#     if m >= 0.5:
#         kroupa_imf_sec.append(m**(-2.3))

# kroupa_imf_sec = np.asarray(kroupa_imf_sec)
# eta = N / np.sum(kroupa_imf)
# eta_sec = N / np.sum(kroupa_imf_sec)

# imf = eta * kroupa_imf
# imf_sec = eta_sec * kroupa_imf_sec
# plt.plot(masses, imf)
# plt.show()
# print(np.sum(imf[103:]))
# print(np.sum(imf))
# print(np.sum(imf_sec[103:]))
# print(np.sum(imf_sec))

# plt.bar(np.linspace(0, N,N), height=imf)
# plt.show()
# imf_sec_counts, bins =  np.histogram(imf_sec)
# plt.hist(bins[:-1], bins, weights=imf_sec_counts)
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt

# masses = np.linspace(5, 1000, 150)  # 0.01 to 100 M☉

# def kroupa_imf(m):
#     return np.piecewise(m,
#         [m <= 0.08, (m > 0.08) & (m < 0.5), m >= 0.5],
#         [lambda m: m**(-0.3), lambda m: m**(-1.3), lambda m: m**(-2.3)]
#     )

# imf = kroupa_imf(masses)
# imf /= np.sum(imf)  # Normalize (or integrate properly for continuous)

# plt.figure()
# plt.plot(masses, imf)
# # plt.xscale('log')
# # plt.yscale('log')
# plt.xlabel('Mass [$M_\odot$]')
# plt.ylabel('IMF (normalized)')
# plt.title('Kroupa IMF')
# plt.grid(True)
# plt.show()


# Define IMF and normalized PDF again for reference
def kroupa_imf_continuous(m):
    if m <= 0.08:
        return m**(-0.3)
    elif m < 0.5:
        return m**(-1.3)
    else:
        return m**(-2.3)

m_min, m_max = 5, 150
masses = np.logspace(np.log10(m_min), np.log10(m_max), 1000)
normalization, _ = quad(kroupa_imf_continuous, m_min, m_max)

def kroupa_pdf(m):
    return kroupa_imf_continuous(m) / normalization

# Calculate number of stars between 25 and 150 M_sun
m1a, m2a = 25,60
fractiona, _ = quad(kroupa_pdf, m1a, m2a)

N_total = 100000
N_in_range_1 = fractiona * N_total
print(f"Expected number of stars between {m1a}-{m2a} M☉ = {N_in_range_1  :.2f}, percentage = {fractiona*100:.2f}%")

m1b, m2b = 12,150
fractionb, _ = quad(kroupa_pdf, m1b, m2b)
N_in_range_2 = fractionb * N_total
print(f"Expected number of stars between {m1b}-{m2b} M☉ = {N_in_range_2  :.2f}, percentage = {fractionb*100:.2f}%")

fractionc, _ = quad(kroupa_pdf, m_min, m_max)
print(f"Expected number of stars between {m_min}-{m_max} M☉ = {N_total:.2f}, percentage = {fractionc*100:.2f}%")

pdf_values = np.array([kroupa_pdf(m) for m in masses])

plt.figure(figsize=(3.5,3.5))
plt.bar(masses, pdf_values, label='Normalized Kroupa IMF', color='gray', width=np.diff(masses, prepend=masses[0]), alpha=0.5)
plt.fill_between([m for m in masses if m1a <= m <= m2a], 0,np.max(pdf_values),alpha=0.3, label='25-60 M☉ Range', hatch='xx', color='blue')
plt.fill_between([m for m in masses if m1b <= m <= m2b], 0,np.max(pdf_values),alpha=0.3, color='red', label='20-150 M☉ Range')
plt.text(16, 0.1, fr'$\Delta$m(25-60) = {fractiona*100:.2f}%')
plt.text(16, 0.05, fr'$\Delta$m(20-150) = {fractionb*100:.2f}%')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Mass [$M_\odot$]')
plt.ylabel('Probability Density [1/$M_\odot$]')
plt.legend()
plt.title('Kroupa IMF')
plt.grid(True, which='both', ls='--')
plt.tight_layout()
plt.savefig('Plots/kroupa_imf_test.pdf', bbox_inches='tight')
plt.close()

# Generate grid for masses and compute PDF
masses = np.logspace(np.log10(m_min), np.log10(m_max), 10000)
pdf_values = np.array([kroupa_pdf(m) for m in masses])

cdf = np.cumsum(pdf_values * np.gradient(masses))
cdf = cdf - cdf[0]   # force start at 0 exactly
cdf /= cdf[-1] 
# Create inverse CDF interpolator
inverse_cdf = interp1d(cdf, masses, kind='linear') 

N_samples = 100000
u = np.random.rand(N_samples)
sampled_masses = inverse_cdf(u)

# Plot histogram of sampled masses
plt.figure(figsize=(3.5,3.5))
plt.hist(sampled_masses, bins=np.logspace(np.log10(m_min), np.log10(m_max), 50),
         density=True, alpha=0.7, color='skyblue', label='Sampled Distribution')
plt.bar(masses, pdf_values, label='Normalized Kroupa IMF', color='gray', width=np.diff(masses, prepend=masses[0]), alpha=0.3)
# plt.fill_between([m for m in masses if m1a <= m <= m2a], 0,np.max(pdf_values),alpha=0.3, label='25-60 M☉ Range', hatch='xx', color='blue')
# plt.fill_between([m for m in masses if m1b <= m <= m2b], 0,np.max(pdf_values),alpha=0.3, color='red', label='20-150 M☉ Range')
plt.text(8, 0.02, fr'N($\Delta$m({m1a}-{m2a})) M☉: {np.sum((sampled_masses >= m1a) & (sampled_masses <= m2a))}')
plt.text(8, 0.01, fr'N($\Delta$m({m1b}-{m2b})) M☉: {np.sum((sampled_masses >= m1b) & (sampled_masses <= m2b))}')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Mass [$M_\odot$]')
plt.ylabel('Probability Density [1/$M_\odot$]')
plt.legend()
plt.title(f'Sampled from Kroupa IMF (N:{N_samples})')
plt.grid(True, which='both', ls='--')
plt.tight_layout()
plt.savefig('Plots/kroupa_imf_sampled.pdf', bbox_inches='tight')
plt.close()
# Print the first 10 sampled masses
print("First 10 sampled masses:", sampled_masses[:10])
# Print the number of stars in the sampled range
print(f"Number of stars in the range {m1a}-{m2a} M☉: {np.sum((sampled_masses >= m1a) & (sampled_masses <= m2a))}")
# Print the number of stars in the sampled range
print(f"Number of stars in the range {m1b}-{m2b} M☉: {np.sum((sampled_masses >= m1b) & (sampled_masses <= m2b))}")
# Print the total number of sampled stars
print(f"Total number of sampled stars: {len(sampled_masses)}")
# Print the total number of stars in the sampled range
print(r'For interval $[25,60]$ M☉: f$_{BH}$ ', fractiona/fractionc)
print(r'For interval $[20,150]$ M☉: f$_{BH}$ ', fractionb/fractionc)

