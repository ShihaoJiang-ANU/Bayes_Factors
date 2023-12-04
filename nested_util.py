
import numpy
import scipy.stats
import scipy.special

def mkloglikelihood(model, data, sigma, xmin, xmax):

    ndata, _ = data.shape
    lnorm = -0.5 * ndata * numpy.log(2.0*numpy.pi) - numpy.log( sigma)*ndata
    
    def loglikelihood(x):

        eval = model(x, xmin, xmax).get_evaluator()
        r = eval(data[:, 0]) - data[:, 1]
        
        return -0.5 * r.dot(r)/sigma**2 + lnorm

    return loglikelihood

def group_data(data):

    r, _ = data.shape

    gddata = {}
    for i in range(r):
        x, y = data[i, :]
        if not (x in gddata):
            gddata[x] = []
        gddata[x].append(y)

    # Pre-convert to numpy array to prevent doing it in likelihood loop
    gdata = {}
    for k, v in gddata.items():
        gdata[k] = numpy.array(v)

    return gdata

def mkgrouploglikelihood(model, gdata, sigma, xmin, xmax):

    if not (type(gdata) is dict):
        raise Exception('Expected dictionary for data parameter')
    
    def loglikelihood(x):
        eval = model(x, xmin, xmax).get_evaluator()

        like = 0.0
        
        for xi, ys in gdata.items():

            mui = eval(xi)


            Ki = float(ys.size)
            Qi = numpy.sum(((mui - ys)/sigma)**2)

            #like += (-Ki/2.0 * numpy.log(2.0)) + (-(scipy.special.gammaln(Ki/2.0))) + ((Ki/2.0)*numpy.log(Qi)) - Qi/2.0
            like += -Ki/2.0 * numpy.log(2.0*numpy.pi) - Ki*numpy.log(sigma) - 0.5*Qi
            
        return like
    
    return loglikelihood

def mkhierarchicalgrouploglikelihood(model, gdata, xmin, xmax):

    if not (type(gdata) is dict):
        raise Exception('Expected dictionary for data parameter')
    
    def loglikelihood(x):
        sigma = x[0]
        eval = model(x[1:], xmin, xmax).get_evaluator()

        like = 0.0
        
        for xi, ys in gdata.items():

            mui = eval(xi)


            Ki = float(ys.size)
            Qi = numpy.sum(((mui - ys)/sigma)**2)

            #like += (-Ki/2.0 * numpy.log(2.0)) + (-(scipy.special.gammaln(Ki/2.0))) + ((Ki/2.0)*numpy.log(Qi)) - Qi/2.0
            like += -Ki/2.0 * numpy.log(2.0*numpy.pi) - Ki*numpy.log(sigma) - 0.5*Qi
            
        return like
    
    return loglikelihood
    
    
def mkchiloglikelihood(model, gdata, sigma, xmin, xmax):

    if not (type(gdata) is dict):
        raise Exception('Expected dictionary for data parameter')
    
    def loglikelihood(x):
        eval = model(x, xmin, xmax).get_evaluator()

        like = 0.0
        
        for xi, ys in gdata.items():

            mui = eval(xi)

            
            Ki = ys.size
            Qi = numpy.sum(((mui - ys)/sigma)**2)

            #like += (-Ki/2.0 * numpy.log(2.0)) + (-(scipy.special.gammaln(Ki/2.0))) + ((Ki/2.0)*numpy.log(Qi)) - Qi/2.0
            like += scipy.stats.chi2.logpdf(Qi, Ki)
            
        return like
    
    return loglikelihood

def mkpriortransform(mu, std):

    def priortransform(u):
        return scipy.stats.norm.ppf(u, loc = mu, scale = std)

    return priortransform


def mkhierarchicalpriortransform(noisemax, mu, std):

    def priortransform(u):

        x = numpy.ones((len(mu) + 1,)) * noisemax * u[0]
        x[1:] = scipy.stats.norm.ppf(u[1:], loc = mu, scale = std)

        return x

    return priortransform

#def gradient(x):
#    """Gradient of multivariate normal log-likelihood."""
#    return -np.dot(Cinv, x) / stats.norm.pdf(x)

def normalize_and_confidence(density, cutoff = 0.9):

    ndensity = density/numpy.sum(density)

    i = 0
    j = ndensity.size - 1

    while i < j and numpy.sum(ndensity[i:j + 1]) > cutoff:
        if ndensity[i] < ndensity[j]:
            i += 1
        else:
            j -= 1

    return ndensity, i, j

def build_curve_histogram(weights, samples, model,
                          xmin, xmax, xbins,
                          ymin, ymax, ybins):

    indices = numpy.where(weights > 0.0)[0]

    mapi = numpy.argmax(weights)
    
    image = numpy.zeros((ybins, xbins))

    xcentres = numpy.linspace(xmin, xmax, xbins + 1)[:-1]
    xcentres += (xcentres[1] - xcentres[0])/2.0

    ycentres = numpy.linspace(ymin, ymax, ybins + 1)[:-1]
    ycentres += (ycentres[1] - ycentres[0])/2.0
    
    xmu = numpy.zeros((xbins,))
    xacc = numpy.zeros((xbins,))
    mapmodel = numpy.zeros((xbins,))
    
    for i in indices:
        e = model(samples[i, 1:], xmin, xmax).get_evaluator()

        w = weights[i]

        for j in range(xbins):
            y = e(xcentres[j])

            yi = numpy.clip(int((y - ymin)/(ymax - ymin) * ybins), 0, ybins - 1)

            image[yi, j] += w

            xacc[j] += w
            xmu[j] += w*y

        if i == mapi:
            for j in range(xbins):
                mapmodel[j] = e(xcentres[j])

            

    conf_min = numpy.zeros((xbins,))
    conf_max = numpy.zeros((xbins,))
    for j in range(xbins):
        image[:, j], low, high = normalize_and_confidence(image[:, j])

        conf_min[j] = ycentres[low]
        conf_max[j] = ycentres[high]
    
    return image, xcentres, ycentres, numpy.divide(xmu, xacc), mapmodel
        
