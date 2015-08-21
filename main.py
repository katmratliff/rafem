#! /usr/bin/env python
import numpy as np


def plot_elevation(avulsion):
    import matplotlib.pyplot as plt

    z = avulsion.get_value('land_surface__elevation')

    plt.imshow(z, origin='lower', cmap='terrain')
    plt.colorbar().ax.set_label('Elevation (m)')
    plt.show()


def main():
    import argparse
    from avulsion_bmi import BmiRiverModule

    parser = argparse.ArgumentParser('Run the avulsion model')
    parser.add_argument('file', help='YAML-formatted parameters file')
    parser.add_argument('--days', type=int, default=0,
                        help='Run model for DAYS')
    parser.add_argument('--years', type=int, default=0,
                        help='Run model for YEARS')
    parser.add_argument('--plot', action='store_true',
                        help='Plot final elevations')

    args = parser.parse_args()

    np.random.seed(1945)

    avulsion = BmiRiverModule()
    avulsion.initialize(args.file)

    n_steps = int((args.days + args.years * 365.) / avulsion.get_time_step())
    for _ in xrange(n_steps):
        avulsion.update()

    if args.plot:
        plot_elevation(avulsion)

    avulsion.finalize()


if __name__ == '__main__':
    main()
