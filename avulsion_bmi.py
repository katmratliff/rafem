#! /usr/local/bin/python
""" Basic Model Interface implementation for River Module"""

from bmi import bmi
from avulsion_main import RiverModule


class BmiRiverModule(Bmi):

    _name = 'Avulsion Module'
    _input_var_names = ('sea_shoreline')
    # not sure what's the most appropriate CSDMS river mouth location name?
    _output_var_names = ('river_mouth_location',
                         'channel_water_sediment~bedload__volume_flow_rate')

    def __init__(self):
    	"""Creat a BmiRiver module that is ready for initialization."""
    	self._model = None
        self._values = {}
        self._var_units = {}

    def initialize(self, filename):
    	"""Initialize the River module"""
    	self._model = RiverModule.params_from_file(filename)

    	self._values = {
            'river_mouth_location': self._modelriv_mouth,
            'channel_water_sediment~bedload__volume_flow_rate': self._model.sed_flux
        }

        self._var_units = {
            'channel_water_sediment~bedload__volume_flow_rate': "m^3 s^-1"
        }

    def update(self):
    	"""Advance model by one time step."""
    	self._model.advance_in_time()
    	