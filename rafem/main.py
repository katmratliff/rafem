#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import shutil

import click
import numpy as np
from six.moves import range


def empty_bmi_var_array(bmi, name):
    return np.empty(bmi.get_var_nbytes(name), dtype=np.uint8).view(
        dtype=bmi.get_var_type(name)
    )


def plot_elevation(avulsion):
    import matplotlib.pyplot as plt

    z = empty_bmi_var_array(avulsion, "land_surface__elevation")
    avulsion.value("land_surface__elevation", z)

    plt.imshow(z, origin="lower", cmap="terrain")
    plt.colorbar().ax.set_label("Elevation (m)")
    plt.show()


def plot_profile(avulsion):
    import matplotlib.pyplot as plt

    prof = empty_bmi_var_array(avulsion, "channel_centerline__elevation")
    avulsion.value("channel_centerline__elevation", prof)

    plt.plot(prof)
    plt.show()


class RafemOutputWriter:
    def __init__(self, bmi, run_id=None, output_interval=1):
        self._bmi = bmi
        self._run_id = run_id
        self._output_interval = output_interval
        self._steps = 0

        self._outputs = {
            "elev_grid": os.path.join(self.prefix, "elev_grid{0}".format(self._run_id)),
            "riv_course": os.path.join(
                self.prefix, "riv_course{0}".format(self._run_id)
            ),
            "riv_profile": os.path.join(
                self.prefix, "riv_profile{0}".format(self._run_id)
            ),
        }

        self._make_dirs()
        self._make_buffers()

        with open(os.path.join(self.prefix, "input.yaml"), "w") as fp:
            print(bmi._model.to_yaml(), file=fp)

    def _make_dirs(self):
        os.mkdir(self.prefix)
        for dir_ in self._outputs.values():
            os.mkdir(dir_)

    def _make_buffers(self):
        self._z = empty_bmi_var_array(self._bmi, "land_surface__elevation")
        self._x = empty_bmi_var_array(self._bmi, "channel_centerline__x_coordinate")
        self._y = empty_bmi_var_array(self._bmi, "channel_centerline__y_coordinate")
        self._prof = empty_bmi_var_array(self._bmi, "channel_centerline__elevation")

    def _update_buffers(self):
        self._bmi.get_value("land_surface__elevation", self._z)
        self._bmi.get_value("channel_centerline__x_coordinate", self._x)
        self._bmi.get_value("channel_centerline__y_coordinate", self._y)
        self._bmi.get_value("channel_centerline__elevation", self._prof)

    def time_stamp(self):
        return str(self._steps)

    def _save_output(self):
        self._update_buffers()

        np.savetxt(
            os.path.join(
                self._outputs["elev_grid"], "elev_{0}.out".format(self.time_stamp())
            ),
            self._z,
            fmt="%.5f",
        )
        np.savetxt(
            os.path.join(
                self._outputs["riv_course"], "riv_{0}.out".format(self.time_stamp())
            ),
            list(zip(self._x, self._y)),
            fmt="%.5f",
        )
        np.savetxt(
            os.path.join(
                self._outputs["riv_profile"], "prof_{0}.out".format(self.time_stamp())
            ),
            self._prof,
            fmt="%.5f",
        )

    @property
    def prefix(self):
        return "run{0}".format(self._run_id)

    def update(self, n_steps):
        self._steps += n_steps
        if self._steps % self._output_interval == 0:
            self._save_output()


@click.command()
@click.version_option()
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="don't emit messages to stderr unless they are errors",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="print input parameters that would be used but don't actually do anything.",
)
@click.option("--days", type=int, default=0, help="number of days to run model")
@click.option("--years", type=int, default=0, help="number of years to run model")
@click.option("--plot-elev/--no-plot-elev", default=False, help="plot final elevations")
@click.option("--plot-prof/--no-plot-prof", default=False, help="plot final profile")
@click.option("--save/--no-save", default=False, help="save output files")
@click.option(
    "--spacing", type=int, default=1, help="spacing for saved files (timesteps)"
)
@click.option("--run-id", type=int, default=1, help="experiment id number")
@click.argument(
    "file", type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True)
)
def main(
    file, quiet, dry_run, days, years, plot_elev, plot_prof, save, spacing, run_id
):
    from .riverbmi import BmiRiverModule

    avulsion = BmiRiverModule()
    avulsion.initialize(file)

    n_days = days + years * 365
    n_steps = int(n_days / avulsion.get_time_step())

    if dry_run:
        click.secho(avulsion._model.to_yaml(), err=False)
    elif n_steps == 0:
        click.secho("Nothing to do (years == 0). 😴", err=True, fg="green")
    else:
        if save:
            output = RafemOutputWriter(avulsion, run_id=run_id, output_interval=spacing)

        with click.progressbar(
            range(n_steps),
            label=" ".join(["🚀", os.path.basename(file)]),
            item_show_func=lambda step: "day {0} of {1}".format(
                int(0 if step is None else step * avulsion.get_time_step()), n_days
            ),
        ) as bar:
            for step in bar:
                avulsion.update()
                save and output.update(1)

        not quiet and click.secho("💥 Finished! 💥", err=True, fg="green")
        if save and not quiet:
            click.secho(
                "Output written to {0}".format(output.prefix), err=True, fg="green"
            )

        plot_elev and plot_elevation(avulsion)
        plot_prof and plot_profile(avulsion)

    avulsion.finalize()
