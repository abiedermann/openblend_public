"""
:Authors: Andrew Biedermann, Isabella Gengaro
:Year: 2020
:Copyright: MIT License
"""

from opentrons import protocol_api
import pandas as pd
import numpy as np

class stock_solution_plate(object):
    """
    Common wrapper for all stock solution plate setups
    """

    def __init__(self,p,sheet,location):
        self.plate_type = sheet[0][0]
        self.loc = location
        self.plate = p.load_labware(self.plate_type,self.loc)
        self.df = sheet[3:].fillna(0)
        self.df.columns = sheet.iloc[2]
        self.df.set_index("Location",inplace=True)
        self.df_volume = None
        print(self.df)

class wellplate24(object):
    """
    Common wrapper for all 24-well-plate objects
    """
    def __init__(self,p,sheet,location):
        self.plate_type = sheet[0][0]
        self.loc = location
        self.plate = p.load_labware(self.plate_type,self.loc)
        self.df = sheet[3:]
        self.df.columns = sheet.iloc[2]
        self.df.fillna(0,inplace=True)
        self.df.set_index("Location",inplace=True)
        self.df_volume = None
        print(self.df)

class water_plate(object):
    """
    Common wrapper for all water plate setups
    """

    def __init__(self,p,location):
        self.plate_type = "agilent_1_reservoir_290ml" 
        self.loc = location
        self.plate = p.load_labware(self.plate_type,self.loc)

