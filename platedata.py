import os
import numpy
from pandas import read_excel
import matplotlib.pyplot as P

class Data:

    def __init__(self, ycol, global_plate=False, basin_plate=False, petrolog=False, seismic=False, seis_correction=False):

        if seismic:
            age = read_excel('file://localhost%s/%s' % (os.getcwd(), 'OIB_compiled_location_x_seis.xlsx'))
            xcol = 'SL2013sv'
        elif seis_correction:
            age = read_excel('file://localhost%s/%s' % (os.getcwd(), 'OIB_compiled_location_x_seis.xlsx'))
            xcol = 'x_mean'
        else:
            age = read_excel('file://localhost%s/%s' % (os.getcwd(), 'OIB_compiled_location_x_basin.xlsx'))
            if basin_plate:
                xcol = 'x_mean_plate'
            if global_plate:
                xcol = 'global_mean_plate'

        if petrolog:
            data = read_excel('file://localhost%s/%s' % (os.getcwd(), 'OIB_geochemistry_petrolog_final2.xlsx'))
        else:
            data = read_excel('file://localhost%s/%s' % (os.getcwd(), 'OIB_geochemistry_final2.xlsx'))

        x = []
        y = []
        xerr = []

        print(enumerate(data['Location']))
        for j, l in enumerate(data['Location']):
            i = numpy.where(age['Island'] == l)[0]
            if len(i) == 0:
                raise Exception('Failed to find %s' % l)
            
            i = i[0]
            mu = age[xcol][i]
            if numpy.isnan(mu):
                print('Location %s has nan' % l)
            else:
                if (ycol in ['P2O5', 'Yb', 'Lu', 'Th']) and (data[ycol][j] > 0) and (data[ycol][j] < 50):
                    y.append(data[ycol][j])
                elif ycol in ['P2O5', 'Yb', 'Lu', 'Th']:
                    continue
                else:
                    y.append(data[ycol][j])
                x.append(mu)

                if seismic:
                    xerr.append(age['x_std'][i])
                else:
                    xerr.append(age['global_std_plate'][i])
                

        self.ndata = len(x)
        self.xy = numpy.zeros((self.ndata, 2))
        self.xy[:, 0] = x
        self.xy[:, 1] = y

        groups = {}
        for xi, yi in zip(self.xy[:, 0], self.xy[:, 1]):

            if not (xi in groups):
                groups[xi] = []

            groups[xi].append(yi)

        self.groups = {}
        ngroup = 0
        groupstd = []
        for k, v in groups.items():
            v = numpy.array(v)
            self.groups[k] = (v.size, numpy.mean(v**2), numpy.mean(v))
            s = numpy.std(v)
            # print the lithospheric thickness (k)
            # and standard variation of each chemical diagnostic (s) on each island
            # print('%10.6f std %10.6f' % (k, s))

            ngroup += 1
            groupstd.append(s)
        print('%d %10.6f' % (ngroup, numpy.sum(groupstd)/ngroup))

        

    def likelihood(self, model, sigma_x, sigma_y):

        den = 2.0 * sigma_y**2

        ev = model.get_evaluator()

        like = 0.0
        
        for i in range(self.N):

            py = ev(self.x[i])

            like += ((self.y[i] - py)**2)/den

        norm = 0.5*self.ndata * numpy.log(den)
        return like, norm
    
    def likelihood_groups(self, model, sigma_x, sigma_y):

        den = 2.0 * sigma_y**2

        ev = model.get_evaluator()

        like = 0.0

        for x, (N, v2, v) in self.groups.items():

            py = ev(x)

            like += (N * (py**2 + v2 - 2.0*py*v))/den

        norm = 0.5*self.ndata * numpy.log(den)
        return like, norm

    def get_points(self):

        return self.xy[:, 0], self.xy[:, 1]
    
if __name__ == '__main__':

    #
    # Figure 10/11:
    #variable = 'SiO2'
    variable = 'P2O5'
    # variable = 'TiO2'
    #
    #variable = 'TiO2'
    #variable = 'lambda1'
    #variable = 'λ1'
    data = Data(variable, xcol = 'x_mean_plate')
    
    fig, ax = P.subplots()

    x, y = data.get_points()

    ax.scatter(x, y)

    P.show()
    
    

        
        

            
        
    
