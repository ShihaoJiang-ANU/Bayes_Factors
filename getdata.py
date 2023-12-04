import argparse

import numpy
import matplotlib.pyplot as P

import platedata

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--column', type=str, required=True, help='Column Name')

    parser.add_argument('-G', '--global-plate', action='store_true', default=False, help='Use global plate model')
    parser.add_argument('-B', '--basin-plate', action='store_true', default=False, help='Use basin-based plate model')

    parser.add_argument('-P', '--petro', action='store_true', default=False, help='Use petrolog correction')

    parser.add_argument('-L', '--seismic', action='store_true', default=False,
                        help='Use seismic lithospheric thickness')
    parser.add_argument('-C', '--seis-correction', action='store_true', default=False,
                        help='Use seismic lithospheric thickness with correction')

    parser.add_argument('-S', '--show', action = 'store_true', default = False, help = 'Plot data')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')
    
    args = parser.parse_args()

    data = platedata.Data(args.column,
                          basin_plate = args.basin_plate,
                          global_plate = args.global_plate,
                          petrolog = args.petro,
                          seismic = args.seismic,
                          seis_correction = args.seis_correction)

    numpy.savetxt(args.output, data.xy)

    if args.show:
        fig, ax = P.subplots()

        ax.scatter(data.xy[:, 0], data.xy[:, 1])

        P.show()

    
    

    
