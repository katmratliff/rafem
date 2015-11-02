#! /usr/local/bin/python
""" Basic Model Interface implementation for River Module"""

import numpy as np

from bmi import Bmi
from rivermodule import RiverModule


class BmiRiverModule(Bmi):

    _name = 'Avulsion Module'
    _input_var_names = ()
    _output_var_names = ('channel_centerline__x_coordinate',
                         'channel_centerline__y_coordinate',
                         'channel_water_sediment~bedload__mass_flow_rate',
                         'channel_exit__x_coordinate',
                         'channel_exit__y_coordinate',
                         'land_surface__elevation',
                         'channel_profile'
                         'avulsion_record'
                        )
    def __init__(self):
        """Create a BmiRiver module that is ready for initialization."""
        self._model = None
        self._values = {}
        self._var_units = {}

    def initialize(self, filename):
        """Initialize the River module"""
        self._model = RiverModule.from_path(filename)

        self._values = {
            'channel_centerline__x_coordinate': lambda: self._model.river_x_coordinates,
            'channel_centerline__y_coordinate': lambda: self._model.river_y_coordinates,
            'channel_water_sediment~bedload__volume_flow_rate': lambda: self._model.sediment_flux,
            'channel_exit__x_coordinate': lambda: self._model.river_x_coordinates[-1],
            'channel_exit__y_coordinate': lambda: self._model.river_y_coordinates[-1],
            'land_surface__elevation': lambda: self._model.elevation,
            'channel_profile': lambda: self._model.profile,
            'avulsion_record': lambda: self._model.avulsions,
        }

        self._var_units = {
            'channel_centerline__x_coordinate': 'm',
            'channel_centerline__y_coordinate': 'm',
            'channel_water_sediment~bedload__volume_flow_rate': "kg s^-1",
            'channel_exit__x_coordinate': 'm',
            'channel_exit__y_coordinate': 'm',
            'land_surface__elevation': 'm',
            'channel_profile': 'm',
            'avulsion_record': 'none',
        }

    def update(self):
        """Advance model by one time step."""
        self._model.advance_in_time()
        self.river_mouth_location = (self._model.river_x_coordinates[-1],
                                     self._model.river_y_coordinates[-1])

    def update_frac(self, time_frac):
        """Update model by a fraction of a time step."""
        time_step = self.get_time_step()
        self._model.time_step = time_frac * time_step
        self.update()
        self._model.time_step = time_step

    def update_until(self, then):
        """Update model until a particular time."""
        n_steps = (then - self.get_current_time()) / self.get_time_step()

        for _ in xrange(int(n_steps)):
            self.update()
        self.update_frac(n_steps - int(n_steps))

    def finalize(self):
        """Clean up & save avulsion file"""
        pass

    def get_var_type(self, var_name):
        """Data type of variable."""
        return str(self.get_value(var_name).dtype)

    def get_var_units(self, var_name):
        """Get units of variable."""
        return self._var_units[var_name]

    def get_var_nbytes(self, var_name):
        """Get units of variable."""
        return self.get_value(var_name).nbytes

    def get_value(self, var_name):
        """Copy of values."""
        return self._values[var_name]()

    def set_value(self, var_name, new_vals):
        """Set model values."""
        if var_name == 'land_surface__elevation':
            self._model.elevation[:] = new_vals[:]

    def get_component_name(self):
        """Name of the component."""
        return self._name

    def get_input_var_names(self):
        """Get names of input variables."""
        return self._input_var_names

    def get_output_var_names(self):
        """Get names of output variables."""
        return self._output_var_names

    def get_start_time(self):
        """Start time of model."""
        return 0.

    def get_end_time(self):
        """End time of model."""
        return np.finfo('d').max

    def get_current_time(self):
        """Current time of model."""
        return self._model.time

    def get_time_step(self):
        """Time step of model."""
        return self._model.time_step

    def get_time_units(self):
        """Units of time"""
        return 'd'
