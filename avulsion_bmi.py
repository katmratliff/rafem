#! /usr/local/bin/python
""" Basic Model Interface implementation for River Module"""

from bmi import bmi
from avulsion_main import RiverModule


class BmiRiverModule(Bmi):

    _name = 'Avulsion Module'
    _input_var_names = ('sea_shoreline')
    # not sure what's the most appropriate CSDMS river mouth location name?
    _output_var_names = ('river_mouth_location',
                         'channel_water_sediment~bedload__mass_flow_rate')

    def __init__(self):
    	"""Creat a BmiRiver module that is ready for initialization."""
    	self._model = None
        self._values = {}
        self._var_units = {}

    def initialize(self, filename):
    	"""Initialize the River module"""
    	self._model = RiverModule.params_from_file(filename)

    	self._values = {
            'river_x_coordinates': self._model.river_x_coordinates,
            'river_y_coordinates': self._model.river_y_coordinates,
            'channel_water_sediment~bedload__volume_flow_rate': self._model.sed_flux,
            'sea_shoreline': self._model._shoreline
        }

        self._var_units = {
            'channel_water_sediment~bedload__volume_flow_rate': "kg s^-1"
        }

    def update(self):
    	"""Advance model by one time step."""
    	self._model.advance_in_time()
    	river_mouth_location = (self._model.river_x_coordinates,
    							self._model.river_y_coordinates)

    def finalize(self):
        """Clean up & save avulsion file"""
        
        if self._model.savefiles == 1:
            np.savetxt('avulsions', self._model.avulsions, fmt='%i %i %i %.3f %.3f %.3f')
        pass