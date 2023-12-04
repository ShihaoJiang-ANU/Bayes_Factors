
import sys
import argparse

import numpy
import matplotlib.pyplot as P

import dynesty

import model
import nested_util

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input data file')

    parser.add_argument('-x', '--xmin', type = float, default = 0.0, help = 'Min x')
    parser.add_argument('-X', '--xmax', type = float, default = 130.0, help = 'Max x')

    parser.add_argument('-y', '--ymin', type = float, default = -10.0, help = 'Min y')
    parser.add_argument('-Y', '--ymax', type = float, default = 10.0, help = 'Max y')

    parser.add_argument('-b', '--xbins', type = int, default = 100, help = 'Hist bins x')
    parser.add_argument('-B', '--ybins', type = int, default = 100, help = 'Hist bins x')

    parser.add_argument('-T', '--threshold', type = int, default = 12, help = 'Exp weight threshold')
    
    parser.add_argument('-m', '--model', type = int, default = 0, help = 'Model')

    parser.add_argument('-n', '--noise-max', type = float, default = 0.0, help = 'Noise std. dev.')

    parser.add_argument('-N', '--nlive', type = int, default = 50, help = 'No. live')

    parser.add_argument('-o', '--output', type = str, default = None, help = 'Output files prefix')

    parser.add_argument('-c', '--change', type = float, default = 60.0, help = 'Change point mean')
    parser.add_argument('-C', '--change-std', type = float, default = 5.0, help = 'Change point std')

    args = parser.parse_args()

    data = numpy.loadtxt(args.input)

    noise_std = numpy.std(data[:, 1])
    print('Mean of data         :', numpy.mean(data[:, 1]))
    print('Std deviation of data:', noise_std)
    
    gdata = nested_util.group_data(data)

#    print('Groups:')
#    for x, ys in gdata.items():
#        print('  %6.3f %d' % (x, ys.size))

    plot_groups = False
    if plot_groups:
        fig, ax = P.subplots()
        ax.scatter(data[:, 0], data[:, 1])
        for x, ys in gdata.items():
            ax.scatter([x] * ys.size, ys, marker = '.')

        P.show()
        sys.exit(0)
        
    model = model.Models[args.model]

    prior_mu, prior_std = model.generate_prior(args.xmin, args.xmax, data, args.change, args.change_std)

    print('Prior:', prior_mu, prior_std)

    noise_max = args.noise_max
    if noise_max <= 0.0:
        noise_max = noise_std * 1.25
        
    pt = nested_util.mkhierarchicalpriortransform(noise_max, prior_mu, prior_std)
    #like = nested_util.mkchiloglikelihood(model, gdata, args.noise, args.xmin, args.xmax)
    like = nested_util.mkhierarchicalgrouploglikelihood(model, gdata, args.xmin, args.xmax)

    ndim = model.dimension() + 1
    print(ndim, args.xmin, args.xmax)

    #test_x = [10.0]
    #print('Like:', like(test_x))

    sampler = dynesty.NestedSampler(like, pt, ndim, nlive = args.nlive, sample = 'rwalk')
    sampler.run_nested(dlogz = 1.0)
    res = sampler.results

    #
    # The important number
    #
    evidence = res.logz[-1]
    print('Evidence:', evidence)
    print(res.samples[-1, :])

    image, xc, yc, mu, mapmodel = nested_util.build_curve_histogram(res.importance_weights(),
                                                                    res.samples, model,
                                                                    args.xmin, args.xmax, args.xbins,
                                                                    args.ymin, args.ymax, args.ybins)
    
    if (args.output is None):

        fig, ax = P.subplots()

        ax.scatter(data[:, 0], data[:, 1], marker = 'x', color = 'black')


        ax.imshow(numpy.log(image), extent = [args.xmin, args.xmax, args.ymin, args.ymax], origin = 'lower',
                  cmap = 'Blues', vmin = -200, vmax = 0)
        
        ax.plot(xc, mu, 'y-')
        ax.plot(xc, mapmodel, 'r-')
        

        P.show()

    else:

        numpy.save(args.output + '_logz.npy', res.logz)
        numpy.save(args.output + '_samples.npy', res.samples)
        numpy.save(args.output + '_weights.npy', res.importance_weights())
        numpy.save(args.output + '_image.npy', image)
        numpy.save(args.output + '_xc.npy', xc)
        numpy.save(args.output + '_yc.npy', yc)
        numpy.save(args.output + '_mu.npy', mu)
        

    

    

    
    
    
