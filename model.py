import numpy
import matplotlib.pyplot as P

class EvaluatorConstant:

    def __init__(self, x, xmin, xmax):

        self.y = x

        self.xmin = xmin
        self.xmax = xmax

    def get_evaluator(self):

        def eval(x):

            return 0.0*x + self.y

        return eval

    def dimension():
        return 1

    def generate_prior(xmin, xmax, data, change_mean, change_std):

        return [numpy.mean(data[:, 1])], [numpy.std(data[:, 1])]
    
class EvaluatorLinear:

    def __init__(self, x, xmin, xmax):

        self.y1, self.y2 = x

        self.xmin = xmin
        self.xmax = xmax

    def get_evaluator(self):

        def eval(x):

            return self.y1 + (self.y2 - self.y1)*(x - self.xmin)/(self.xmax - self.xmin)

        return eval

    def dimension():
        return 2

    def generate_prior(xmin, xmax, data, change_mean, change_std):

        return [numpy.mean(data[:, 1])]*2, [numpy.std(data[:, 1])]*2
    

class EvaluatorTwoLine:

    def __init__(self, x, xmin, xmax):

        self.x1, self.y1, self.y2, self.y3 = x

        self.xmin = xmin
        self.xmax = xmax

    def get_evaluator(self):

        def eval(x):

            if numpy.isscalar(x):

                if x < self.x1:

                    y = self.y1 + (self.y2 - self.y1)*(x - self.xmin)/(self.x1 - self.xmin)

                else:

                    y = self.y2 + (self.y3 - self.y2)*(x - self.x1)/(self.xmax - self.x1)
                    
            else:
                y = self.y1 + (x - self.xmin)*(self.y2 - self.y1)/(self.x1 - self.xmin)

                indices = numpy.where(x > self.x1)[0]
                y[indices] = self.y2 + (self.y3 - self.y2)*(x[indices] - self.x1)/(self.xmax - self.x1)

            return y

        return eval

    def dimension():
        return 4

    def generate_prior(xmin, xmax, data, change_mean, change_std):

        return [change_mean] + [numpy.mean(data[:, 1])]*3, [change_std] + [numpy.std(data[:, 1])]*3
    

Models = [EvaluatorConstant,
          EvaluatorLinear,
          EvaluatorTwoLine]


        
