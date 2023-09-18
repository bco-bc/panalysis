"""Some very simple probability plots
"""

import math
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # Inclination angle (theta) with positive z-axis, no external field.
    bin_size = 0.001
    thetas = np.arange(start=bin_size, stop=math.pi - bin_size, step=bin_size)
    p_theta = np.zeros(shape=len(thetas))
    prob_theta = np.zeros(shape=len(thetas))
    one_over_pi = np.zeros(shape=len(thetas))
    p_cos_theta = np.zeros(shape=len(thetas))
    areas = np.zeros(shape=len(thetas))
    cos_thetas = np.zeros(shape=len(thetas))
    ave_theta = 0.0
    for k in np.arange(0, len(thetas)):
        angle = thetas[k]
        area_k = 4 * math.pi * (1.0 - math.cos(angle))
        area_k_next = 4 * math.pi * (1.0 - math.cos(angle + bin_size))
        areas[k] = math.fabs(area_k - area_k_next)
        p_theta[k] = 0.5 * math.sin(angle)
        prob_theta[k] = p_theta[k] * bin_size
        one_over_pi[k] = 1.0 / math.pi
        cos_theta = math.cos(angle)
        ave_theta += angle * prob_theta[k]
        l = len(thetas) - k - 1
        cos_thetas[l] = cos_theta
        p_cos_theta[l] = p_theta[k] / math.sqrt(1.0 - cos_theta * cos_theta)
        print(f'angle. p_theta, prob_theta, cos_theta, p_cos_theta area: '
              f'{angle, p_theta[k], prob_theta[k], cos_theta, p_cos_theta[k], areas[k]}')
    print(f'Average theta (should be 1/2*math.pi): {ave_theta}')

    print()
    print(f'thetas: {thetas}')
    print(f'p_theta: {p_theta}')
    print(f'cos_thetas: {cos_thetas}')
    print(f'p_cos_theta: {p_cos_theta}')

    figure = plt.figure()

    plt.subplot(3, 1, 1)
    plt.plot(cos_thetas, p_cos_theta, color='red', label=r'p($\theta$)')
    plt.ylabel(r'p(cos($\theta$))')
    plt.xlabel(r'cos($\theta$)')

    plt.subplot(3, 1, 2)
    plt.plot(thetas, p_theta, color='red', label=r'p($\theta$)')
    plt.plot(thetas, prob_theta, color='blue', label=r'p($\theta$)d$\theta$')
    plt.xlabel(r'$\theta$')
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(thetas, areas, color='green')
    plt.ylabel(r'A($\theta$)')
    plt.xlabel(r'$\theta$')

    figure.tight_layout(pad=1.0)
    plt.show()
