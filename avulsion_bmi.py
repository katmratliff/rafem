#! /usr/local/bin/python
""" Basic Model Interface implementation for River Module"""

import numpy as np

from bmi import Bmi
from avulsion_main import RiverModule


class BmiRiverModule(Bmi):

    _name = 'Avulsion Module'
    _input_var_names = ('sea_shoreline')
    # not sure what's the most appropriate CSDMS river mouth stuff?
    _output_var_names = ('river_x_coordinates',
                         'river_y_coordinates',
                         'river_mouth_location',
                         'channel_water_sediment~bedload__mass_flow_rate')

    def __init__(self):
    	"""Create a BmiRiver module that is ready for initialization."""
    	self._model = None
        self._values = {}
        self._var_units = {}

    def initialize(self, filename):
    	"""Initialize the River module"""
    	self._model = RiverModule.params_from_file(filename)

    	self._values = {
            'river_x_coordinates': self._model.river_x_coordinates,
            'river_y_coordinates': self._model.river_y_coordinates,
            'river_mouth_location': self.river_mouth_location,
            'channel_water_sediment~bedload__volume_flow_rate': self._model.sed_flux,
            'sea_shoreline': self._model._shoreline
        }

        self._var_units = {
            'channel_water_sediment~bedload__volume_flow_rate': "kg s^-1"
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
        
        if self._model.savefiles == 1:
            np.savetxt('avulsions', self._model.avulsions,
                       fmt='%i %i %i %.3f %.3f %.3f')
        pass

    def get_var_type(self, var_name):
        """Data type of variable."""
        return str(self.get_value_ref(var_name).dtype)

    def get_var_units(self, var_name):
        """Get units of variable."""
        return self._var_units[var_name]

    def get_var_nbytes(self, var_name):
        """Get units of variable."""
        return self.get_value_ref(var_name).nbytes

    def get_value_ref(self, var_name):
        """Reference to values."""
        return self._values[var_name]

    def get_value(self, var_name):
        """Copy of values."""
        return self.get_value_ref(var_name).copy()

    def get_value_at_indices(self, var_name, indices):
        """Get values at particular indices."""
        return self.get_value_ref(var_name).take(indices)

    def set_value(self, var_name, src):
        """Set model values."""
        val = self.get_value_ref(var_name)
        val[:] = src

    def set_value_at_indices(self, var_name, src, indices):
        """Set model values at particular indices."""
        val = self.get_value_ref(var_name)
        val.flat[indices] = src

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
        """Current time of model (days)."""
        return self._model.time

    def get_time_step(self):
        """Time step of model (days)."""
        return self._model.time_step

