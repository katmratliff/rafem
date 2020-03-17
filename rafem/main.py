#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pathlib
import re
import sys
from collections import OrderedDict
from functools import partial

import click
import matplotlib.pyplot as plt
import numpy as np
import yaml

from ._version import get_versions
from .riverbmi import BmiRiverModule

__version__ = get_versions()["version"]

out = partial(click.secho, bold=True, err=True)
err = partial(click.secho, fg="red", err=True)


def setup_yaml_with_canonical_dict():
    """ https://stackoverflow.com/a/8661021 """
    yaml.add_representer(
        OrderedDict,
        lambda self, data: self.represent_mapping(
            "tag:yaml.org,2002:map", data.items()
        ),
        Dumper=yaml.SafeDumper,
    )

    def repr_ordered_dict(self, data):
        return self.represent_mapping("tag:yaml.org,2002:map", data.items())

    yaml.add_representer(dict, repr_ordered_dict, Dumper=yaml.SafeDumper)

    def repr_dict(self, data):
        return self.represent_mapping(
            "tag:yaml.org,2002:map", sorted(data.items(), key=lambda t: t[0])
        )

    yaml.add_representer(dict, repr_dict, Dumper=yaml.SafeDumper)

    # https://stackoverflow.com/a/45004464
    def repr_str(dumper, data):
        if "\n" in data:
            return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
        return dumper.represent_str(data)

    yaml.add_representer(str, repr_str, Dumper=yaml.SafeDumper)

    def repr_tuple(dumper, data):
        return dumper.represent_sequence("tag:yaml.org,2002:seq", list(data))

    yaml.add_representer(tuple, repr_tuple, Dumper=yaml.SafeDumper)

    yaml.add_implicit_resolver(
        "tag:yaml.org,2002:float",
        re.compile(
            r"""^(?:
         [-+]?(?:[0-9][0-9_]*)\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |[-+]?\.[0-9_]+(?:[eE][-+]?[0-9]+)
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\.[0-9_]*
        |[-+]?\.(?:inf|Inf|INF)
        |\.(?:nan|NaN|NAN))$""",
            re.X,
        ),
        list("-+0123456789."),
        Loader=yaml.SafeLoader,
    )


setup_yaml_with_canonical_dict()


def empty_bmi_var_array(bmi, name):
    return np.empty(bmi.get_var_nbytes(name), dtype=np.uint8).view(
        dtype=bmi.get_var_type(name)
    )


class RafemOutputWriter:
    def __init__(self, bmi, output_interval=1):
        self._bmi = bmi
        # self._run_id = run_id
        self._output_interval = output_interval
        self._steps = 0
        self._prefix = pathlib.Path("output")

        self._outputs = {
            "elevation": self._prefix / "elevation",
            "river": self._prefix / "river",
            "profile": self._prefix / "profile",
        }

        self._make_dirs()
        self._make_buffers()

    def _make_dirs(self):
        os.mkdir(self._prefix)
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
            self._outputs["elevation"] / "elev_{0}.out".format(self.time_stamp()),
            self._z,
            fmt="%.5f",
            header=os.linesep.join(
                ["version: {0}".format(__version__), "Elevation [m]"]
            ),
        )
        np.savetxt(
            self._outputs["river"] / "riv_{0}.out".format(self.time_stamp()),
            list(zip(self._x, self._y)),
            fmt="%.5f",
            header=os.linesep.join(
                ["version: {0}".format(__version__), "X [m], Y [m]"]
            ),
        )
        np.savetxt(
            self._outputs["profile"] / "prof_{0}.out".format(self.time_stamp()),
            self._prof,
            fmt="%.5f",
            header=os.linesep.join(
                ["version: {0}".format(__version__), "Elevation [m]"]
            ),
        )

    @property
    def prefix(self):
        return str(self._prefix)

    def update(self, n_steps):
        self._steps += n_steps
        if self._steps % self._output_interval == 0:
            self._save_output()


@click.group(chain=True)
@click.version_option()
@click.option(
    "--cd",
    default=".",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="chage to directory, then execute",
)
def rafem(cd) -> None:
    """The River Avulsion and Floodplain Evolution Model (RAFEM)

    \b
    Examples:

      Create a folder with example input files,

        $ mkdir example && rafem setup

      Run a simulation using the examples input files,

        $ cd example && rafem run

      Commands can also be chained together. The folling will setup
      a simulation, run it, and then plot elevations.

        $ mkdir example && rafem setup run plot elevation
    """
    os.chdir(cd)


@rafem.command()
@click.option("-v", "--verbose", is_flag=True, help="Emit status messages to stderr.")
@click.option("--dry-run", is_flag=True, help="Do not actually run the model")
def run(dry_run: bool, verbose: bool) -> None:
    """Run a simulation."""
    run_dir = pathlib.Path.cwd()
    config_file = run_dir / "rafem.yaml"

    message = []
    if not config_file.is_file():
        message.append("missing RAFEM configuration file: {0}".format(config_file))
    if (run_dir / "output").exists():
        message.append(
            "RAFEM output directory already exists: {0}".format(run_dir / "output")
        )
    if message:
        err(os.linesep.join(message))
        raise click.Abort(os.linesep.join(message))

    if dry_run:
        out("Nothing to do. 😴")
    else:
        avulsion = BmiRiverModule()
        avulsion.initialize("rafem.yaml")

        with open("rafem.yaml", "r") as fp:
            params = yaml.safe_load(fp)

        spacing = 10
        n_days = params["days"]  # + years * 365
        n_steps = int(n_days / avulsion.get_time_step())

        if dry_run:
            out(avulsion._model.to_yaml())
        elif n_steps == 0:
            out("Nothing to do (years == 0). 😴")
        else:
            output = RafemOutputWriter(avulsion, output_interval=spacing)

            with click.progressbar(
                range(n_steps),
                label=" ".join(["🚀", str(run_dir)]),
                item_show_func=lambda step: "day {0} of {1}".format(
                    int(0 if step is None else step * avulsion.get_time_step()), n_days
                ),
            ) as bar:
                for step in bar:
                    avulsion.update()
                    output.update(1)

            out("💥 Finished! 💥")
            out("Output written to {0}".format(output.prefix))

        avulsion.finalize()


@rafem.command()
@click.argument("infile", type=click.Choice(["rafem"]))
def show(infile: str) -> None:
    """Show example input files."""
    print(_contents_of_input_file(infile))


@rafem.command()
def setup() -> None:
    """Setup a folder of input files for a simulation."""
    files = [pathlib.Path("rafem.yaml")]

    existing_files = [name for name in files if name.exists()]
    if existing_files:
        for name in existing_files:
            err(
                f"{name}: File exists. Either remove and then rerun or setup in a different folder",
            )
    else:
        for fname in files:
            with open(fname, "w") as fp:
                print(_contents_of_input_file(fname.stem), file=fp)

    if existing_files:
        raise click.Abort()


def _contents_of_input_file(infile: str) -> str:

    # Avulsion module parameters
    conf = OrderedDict(
        (
            ("_version", __version__),
            # Space
            ("shape", [120, 200]),  # Length x Width (km)
            ("spacing", [0.1, 0.1]),  # dy, dx (km)
            ("n0", 5.0),  # upstream elevation
            ("nslope", 0.001),  # initial landscape slope
            (
                "max_rand",
                0.1,
            ),  # multiply by slope for max height of a random perturbation
            # Time
            ("days", 7),  # number of days to run for
            ("dt_day", 0.01),  # timestep (days)
            # Random seed
            ("rand_seed", 623),  # seed for random number generator
            # Sea level and subsidence parameters
            ("Initial_SL", 0.0),  # initial sea level
            ("SLRR_m", 0.0),  # sea level rise rate (m/yr)
            ("SubRate_m", 0.0),  # subsidence rate (m/yr)
            ("Sub_Start", 0),  # row where subsidence starts
            # River characteristics
            ("ch_width", 10.0),  # characteristic channel width (m)
            ("ch_depth", 1.0),  # characteristic channel depth (m)
            ("ch_discharge", 10.0),  # long-term averaged discharge (m^3/s)
            ("A", 1.0),  # river-dependent const. (1 for meandering; 1.4 for braided)
            ("c_f", 0.01),  # drag coefficient
            ("C_0", 1.0),  # sediment concentration on bed
            ("sed_sg", 2.65),  # sediment specific gravity
            ("init_cut_frac", 1),  # initial cut of the channel into land surface
            # Avulsion parameters
            ("super_ratio", 1.0),  # normalized SE ratio to trigger avulsion
            ("short_path", 1),  # flag for using shortest path to complete avulsion
            # Floodplain and Wetland parameters
            ("WL_Z", 0.0),  # elevation that wetlands maintain above SL
            ("WL_dist", 0),  # cell distance beyond channel that wetlands exist
            ("blanket_rate_m", 0.0),  # "blanket" deposition rate (m)
            (
                "fine_dep_frac",
                0.0,
            ),  # fraction of channel deposit for adjacent fine deposition
            ("splay_type", 2),  # size/type of splay
            # Splay types:
            # splay_type = 0: no splay deposition
            # splay_type = 1: just the first failed avulsion river cell
            # splay_type = 2: first failed cell and adjacent cells
            # Saving information
            ("saveavulsions", False),  # flag for saving avulsion info
            ("savecourseupdates", False),  # flag for saving course updates
        )
    )

    contents = {
        "rafem": yaml.safe_dump(conf, default_flow_style=False),
    }

    return contents[infile]


@rafem.command()
@click.option("--time", default=-1, type=int, help="time step to plot")
@click.argument("value", type=click.Choice(["elevation", "profile"]))
def plot(value: str, time: int) -> None:
    """Plot output from a simulation."""
    basepath = pathlib.Path("output") / value

    if time == -1:

        def _alphanum_key(key):
            def _convert(text):
                return int(text) if text.isdigit() else text

            return [_convert(c) for c in re.split("([0-9]+)", key)]

        files = sorted(
            [name.name for name in basepath.glob("*.out")], key=_alphanum_key
        )
        filepath = basepath / files[-1]
    else:
        filepath = list(basepath.glob("*_{0}.out".format(time)))[0]
    out(str(filepath))

    with open("rafem.yaml", "r") as fp:
        params = yaml.safe_load(fp)
    shape = params["shape"]

    z = np.loadtxt(filepath, delimiter=",", comments="#")
    if value == "elevation":
        plt.imshow(z.reshape(shape), origin="lower", cmap="terrain")
        plt.colorbar().ax.set_label("Elevation (m)")
    elif value == "profile":
        plt.plot(z)

    plt.show()
