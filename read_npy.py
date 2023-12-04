import os

import numpy
import matplotlib.pyplot as P


def get_slope(element):
    path = os.getcwd()
    xc = numpy.load(path + f'/{element}_global_1_xc.npy')
    mu = numpy.load(path + f'/{element}_global_1_mu.npy')
    slope = (mu[-1] - mu[0]) / (xc[-1] - xc[0])
    print(slope)
    # image = numpy.load(args.input + '_image.npy')
    # xc = numpy.load(args.input + '_xc.npy')
    # yc = numpy.load(args.input + '_yc.npy')
    # mu = numpy.load(args.input + '_mu.npy')


elements = ['Al2O3', 'CaO', 'FeOT', 'K2O', 'Na2O', 'P2O5', 'SiO2', 'TiO2', 'λ0', 'λ1', 'λ2']

for k in elements:
    get_slope(k)
