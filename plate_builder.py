"""
:Authors: Andrew Biedermann, Isabella Gengaro
:Year: 2020
:Copyright: MIT License
"""

from opentrons import protocol_api, robot
from openblend_experiment_wrapper import experiment
import pandas as pd
import numpy as np

class plate_builder(object):
    """
    Executes an experimental plan on the Opentrons.
    """
    def __init__(self,p,xlsx_name):
        """
        Aquires experimental plan, initializes instrument, and derives
        component list
        """
        self.p = p
        robot.reset()
        robot.home()
        self.plan = experiment(p,xlsx_name)
        self.plates = self.plan.get_plates()
        self.stocks = self.plan.get_stocks()
        self.water = self.plan.get_water()
        self.total_volume = 0.003 # L
        self.working_volumes = []
        for s in self.stocks:
            self.working_volumes.append(s.df[["Component","Working Volume (mL)"]])
        self.adjust_working_volumes(self.stocks[0],"A2",1000)

        # initializing instruments
        try:
            self.p10 = self.p.load_instrument('p10_single',mount="left",
                        tip_racks=[self.plan.t10])
        except:
            self.p10 = self.p.load_instrument("p20_single",mount="left",
                        tip_racks=[self.plan.t10])
        self.p300 = self.p.load_instrument('p300_single',mount="right",
                        tip_racks=[self.plan.t300])
        
        self.p10.home()
        self.p300.home()

        # deriving a list of components
        comps = []
        for s in self.stocks:
            comps.append(s.df["Component"].values)
        self.comps = np.unique(np.concatenate(comps).tolist())       
        self.calculate_volumes()

    def calculate_volumes(self):
        """
        Note: an error is thrown if any water additions are negative, as this
        suggests that total volume of the well will exceed the max volume limit
        """
        for this_plate in self.plates:
            this_comps = this_plate.df.columns.values
            this_cdomps = this_comps[this_comps!="Location"]
            this_plate.df_volume = this_plate.df
            for c in this_comps:
                this_plate.df_volume[c] = (this_plate.df_volume[c]*self.total_volume/
                                           self.get_stock_concentration(c))*1e6
            this_plate.df_volume['water'] = self.total_volume*1e6 - this_plate.df_volume.sum(axis=1)
        self.check_for_volume_issues()
        return            

    def get_stock_concentration(self,c):
        for s in self.stocks:
            this_comp = s.df["Component"].values
            if c in this_comp:
                this_conc = s.df["Concentration (g/L)"].values
                this_conc = this_conc[np.where(c==this_comp)[0]]
                break
        return this_conc

    def calculate_total_reagent_volumes(self):
        """
        Prints total volumes needed in microliters
        """
        vols = dict()
        keys = []
        for this_plate in self.plates:
            keys.append(this_plate.df_volume.columns.values[:])
        keys = np.unique(keys)
        for k in keys:
            this_sum = 0
            for this_plate in self.plates:
                this_sum += this_plate.df_volume[k].values.sum()
            vols[k] = this_sum
        print(pd.Series(vols))
        return

    def check_for_volume_issues(self):
        for this_plate in self.plates:
            try:
                w = this_plate.df_volume['water'].values
                assert(np.all(w>0)) # No negative water additions
            except:
                print("Error: Volume exceeded on plate "+str(this_plate.loc))
                print(this_plate.df_volume)
                quit()
            
            try:
                assert(not np.any(np.logical_and(this_plate.df_volume.values>0,
                      this_plate.df_volume.values<1)))
            except:
                print("Error: Sub-microliter volumes on plate "+str(this_plate.loc))
                print(this_plate.df_volume)
                quit()
        return


    def make_plates(self):
        """
        Create media, adding water first
        """
        assert('water' not in self.comps)
        self.p10.pick_up_tip()
        self.p300.pick_up_tip()
        for this_plate in self.plates:
            self.dispense_comp_to_plate(this_plate,self.water,'water','A1')
        self.p10.drop_tip()
        self.p300.drop_tip()
        
        for c in self.comps:
            this_stock=None
            for s in self.stocks:
                this_comp = s.df["Component"].values
                if c in this_comp:
                    this_loc = s.df.index.values
                    this_loc = this_loc[np.where(c==this_comp)[0][0]]
                    this_stock = s
                    break
            self.p10.pick_up_tip()
            self.p300.pick_up_tip()
            for p in self.plates:
                if c not in p.df.columns:
                    continue
                self.dispense_comp_to_plate(p,this_stock,c,this_loc)
            self.p10.drop_tip()
            self.p300.drop_tip()
        robot.home()
        self.p10.home()
        self.p300.home()
        return

    def dispense_comp_to_plate(self,this_plate,stock_plate,
                               comp_name,comp_loc):
        plate_loc = this_plate.df_volume.index.values
        v_add = this_plate.df_volume[comp_name].values
        for i in range(len(v_add)):
            if v_add[i]==0:
                continue
            self.dispense_well(this_plate,stock_plate,plate_loc[i],
                               comp_loc,v_add[i])
        return

    def dispense_well(self,this_plate,stock_plate,ploc,sloc,v_add):
        sloc = str(sloc)
        h = self.get_tip_height(stock_plate,sloc)
        if v_add<=30:
            self.p10.transfer([v_add],source=stock_plate.plate[sloc].bottom(h),
                    dest=this_plate.plate[ploc],blow_out=True,
                    touch_tip=False,carryover=True,new_tip='never')
        else:
            self.p300.transfer([v_add],source=stock_plate.plate[sloc].bottom(h),
                    dest=this_plate.plate[ploc],blow_out=True,
                    touch_tip=False,carryover=True,new_tip='never')
        self.adjust_working_volumes(stock_plate,sloc,v_add)
        return

    def get_working_volume(self,stock_plate,sloc):
        idx = None
        for i,wv in enumerate(self.working_volumes):
            if(np.all(stock_plate.df["Component"].values
            ==wv["Component"].values)):
                idx=i
                return self.working_volumes[idx].loc[sloc,"Working Volume (mL)"]
        return

    def get_tip_height(self,stock_plate,sloc):
        if stock_plate.plate_type == "agilent_1_reservoir_290ml":
            return 1
        wv = self.get_working_volume(stock_plate,sloc)
        if stock_plate.plate[sloc].diameter > 20:
            return self.get_50_height(wv)
        else:
            return self.get_15_height(wv)

    def get_50_height(self,wv):
        dx = (110)/50
        height = (wv-5)*dx
        if height < 2:
            height = 2 # mm
        return height

    def get_15_height(self,wv):
        dx = (110)/15
        height = (wv-1.5)*dx
        if height < 2:
            height = 2 # mm
        return height

    def adjust_working_volumes(self,stock_plate,sloc,v_add):
        if stock_plate.plate_type == "agilent_1_reservoir_290ml":
            return
        idx = None
        for i,wv in enumerate(self.working_volumes):
            if(np.all(stock_plate.df["Component"].values
            ==wv["Component"].values)):
                idx=i
        self.working_volumes[idx].loc[sloc,"Working Volume (mL)"] -= v_add*1e-3
        return




