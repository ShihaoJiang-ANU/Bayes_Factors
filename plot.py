
import os
import argparse

import numpy
import matplotlib.pyplot as P

if __name__ == '__main__':
    def format(x):
        for i in [0, 1, 2, 3, 5]:
            x = x.replace(f'{i}', f'$_{i}$')
        return x

    
    P.rcParams['text.usetex'] = False
    
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input prefix')

    parser.add_argument('-d', '--data', type = str, default = None, help = 'Data points')

    parser.add_argument('-b', '--break-point', action = 'store_true', default = False, help = 'Show mean break point')
    parser.add_argument('-l', '--label', type = str, default = None, help = 'Y axis label')

    parser.add_argument('-p', '--pdf', type = str, default = None, help = 'PDF output')
    
    args = parser.parse_args()

    logz = numpy.load(args.input + '_logz.npy')
    image = numpy.load(args.input + '_image.npy')
    xc = numpy.load(args.input + '_xc.npy')
    yc = numpy.load(args.input + '_yc.npy')
    mu = numpy.load(args.input + '_mu.npy')

    fig, ax = P.subplots()
    fig.set_tight_layout(True)
    fig.set_size_inches((8, 2))

    ax.set_title(r'$\log_{10}p(\mathbf{d})$ = %.3f' % (logz[-1] * numpy.log10(numpy.exp(1))))

    dx = xc[1] - xc[0]

    xmin = xc[0] - dx/2.0
    xmax = xc[-1] + dx/2.0

    dy = yc[1] - yc[0]

    ymin = yc[0] - dy/2.0
    ymax = yc[-1] + dy/2.0
    
    ax.imshow(numpy.log10(image),
              extent = [xmin, xmax, ymin, ymax], origin = 'lower',
              cmap = 'Blues',
              aspect = 'auto')
    ax.plot(xc, mu, 'y-')

    name = args.input.split('_')[-1]
    if 'P' in name:
        name = name.split('P')[0]
    subtitle = {'0': '(a)', '1': '(b)', '2': '(c)'}
    ax.text(0.01, 0.85, subtitle[name], transform=ax.transAxes)

    if not (args.data is None):

        data = numpy.loadtxt(args.data)
        ax.scatter(data[:, 0], data[:, 1], marker = '+', linewidth = 0.5, color = 'black')

    ax.axvline(75.0, linestyle = '--', color = 'red', linewidth = 0.5)
    ax.axvline(85.0, linestyle = '--', color = 'red', linewidth = 0.5)

    ax.set_xlabel('Depth (km)')
    if not (args.label is None):
        ax.set_ylabel(format(args.label))

    figb = None
    bx = None
    if args.break_point:
        samples = numpy.load(args.input + '_samples.npy')
        weights = numpy.load(args.input + '_weights.npy')

        mean_break = numpy.average(samples[:, 1], weights = weights)

        figb, bx = P.subplots()
        figb.set_tight_layout(True)
        figb.set_size_inches((8, 2))

        # bx.hist(samples[:, 1], weights = weights, bins = 100, range = (mean_break - 5.0, mean_break + 5.0), density = True)
        bx.hist(samples[:, 1], weights=weights, bins=300, range=(0, 130), density=True)
        bx.set_xlim(0, 130)
        bx.axvline(mean_break, color = 'yellow')
        bx.text(0.01, 0.85, '(d)', transform=bx.transAxes)
        bx.set_xlabel('Depth (km)')
        bx.set_ylabel('P(break)')


    if args.pdf is None:
        P.show()
    else:
        fig.savefig(args.pdf + '.pdf', format = 'PDF')
        if not (figb is None):
            figb.savefig(os.path.splitext(args.pdf)[0] + '_break.pdf', format = 'PDF')
    
    
