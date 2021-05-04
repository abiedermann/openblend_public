"""
:Author: Andrew Biedermann
:Year: 2020
:Copyright: MIT License
"""

import openturns as ot
import pandas as pd
import numpy as np
import xlsxwriter

class experiment_designer(object):
    """
    Creates an experimental plan, given a dictionary of factors mapped to bounds,
    and the number of plates to use.
    """

    def __init__(self,factor_map=None,constant_comps=None,stock_conc=None,
                 n_plates=None,output_filename=None,input_filename=None,well_volume=0.003):
        """
        Parameters
        ----------
        factor_map: dict, keys are names of components, value is a tuple of
                    minimum and maximum values to explore (inclusive).
        constant_comps: dict, keys are names of components, value is the
                    target concentration
        stock_conc: the concentration of all stock solutions used in this
                    experiment
        n_plates: number of plates available for this experiment
        output_filename: name of the Excel file the plan is written to. Note:
                    This string should have a .xlsx suffix
        """
        if stock_conc==None:
            r=pd.read_excel(input_filename,sheet_name=None,header=None)
            s = r["Stocks"]
            this_stocks = s[3:]
            self.stock_conc = dict()
            comp = this_stocks.iloc[:,0].values
            conc = this_stocks.iloc[:,1].values
            print(conc)
            print(comp)
            for i in range(len(comp)):
                self.stock_conc[comp[i]] = conc[i]
        else:
            self.stock_conc = stock_conc
        self.well_volume = well_volume # L
        self.output_filename = output_filename
        # running sensitivity analysis
        if factor_map:
            self.factor_map = factor_map
            self.constant_comps = constant_comps
            self.n_factors = len(self.factor_map.keys())
            self.n_plates = n_plates
            self.n_samples = 23*self.n_plates # Note: we need one well as a blank
            self.volumes = None
            self.generate_design()
        # running a manually generated experiment
        else:
            self.design=pd.DataFrame()
            self.n_plates=0
            r = pd.read_excel(input_filename,sheet_name=None,header=None)
            for i in range(1,12):
                if str(i) not in r.keys():
                    continue
                s = r[str(i)]
                sheet_type = s[0][0]
                if sheet_type=="corning_24_wellplate_3.4ml_flat":
                    if self.design.empty:
                        self.design = s[3:]
                        self.design.columns = s.iloc[2]
                        self.design.fillna(0,inplace=True)
                        self.design.drop("Location",axis=1,inplace=True)
                        self.n_plates+=1
                    else:
                        this_df = s[3:]
                        this_df.columns = s.iloc[2]
                        this_df.fillna(0,inplace=True)
                        this_df.drop("Location",axis=1,inplace=True)
                        self.design = self.design.append(this_df,ignore_index=True)
                        self.n_plates+=1
            self.n_samples=self.n_plates*23
            # removing empty lines from design
            self.design = self.design.loc[~(self.design==0).all(axis=1)]
        self.calculate_total_reagent_volumes()
        self.configure_stocks_setup()
        return

    def generate_design(self):
        self.design = pd.DataFrame(index=np.arange(self.n_samples))
        if self.n_factors > 0:
            seq = ot.SobolSequence(self.n_factors)
            self.design = pd.DataFrame(np.array(seq.generate(self.n_samples)))
            self.design.columns = self.factor_map.keys()

            # transforming design into the ranges specified by factor_map
            for k in self.design.columns.values:
                self.design[k] = (self.design[k]
                        *(self.factor_map[k][1]-self.factor_map[k][0])
                        +self.factor_map[k][0])
        
        # adding constant components
        for k in self.constant_comps.keys():
            self.design[k] = self.constant_comps[k]*np.ones(len(self.design.index))
        return

    def print_design(self):
        plate_designs = []
        new_index = ["A1","A2","A3","A4","A5","A6",
                     "B1","B2","B3","B4","B5","B6",
                     "C1","C2","C3","C4","C5","C6",
                     "D1","D2","D3","D4","D5"]
                     
        for i in range(self.n_plates):
            this_design = self.design[i*23:(i+1)*23]
            this_design.index = new_index
            plate_designs.append(this_design)

        writer = pd.ExcelWriter(self.output_filename,engine="xlsxwriter")
        for i,df in enumerate(plate_designs):
            this_sheet_name = str(i+1)
            if int(this_sheet_name) > 8:
                print("Error: design calls for too many plates")
                quit()
            df.to_excel(writer,sheet_name=str(this_sheet_name),startrow=2,startcol=0)
            this_sheet = writer.sheets[str(this_sheet_name)]
            this_sheet.write("A1","corning_24_wellplate_3.4ml_flat")
            this_sheet.write("A3","Location")

        for i in range(len(self.stocks)):
            this_sheet_name = str(i+len(plate_designs)+1)
            if int(this_sheet_name) > 8:
                print("Error: design calls for too many plates")
                quit()
            self.stocks[i][1].to_excel(writer,sheet_name=str(this_sheet_name),
                                       startrow=2,startcol=0,index=False)
            this_sheet = writer.sheets[this_sheet_name]
            this_sheet.write("A1",self.stocks[i][0])

        for i in range(9,12):
            this_sheet = writer.book.add_worksheet(str(i))
            if i==9:
                this_sheet.write("A1","water")
            if i==10:
                this_sheet.write("A1","opentrons-tiprack-10ul")
            if i==11:
                this_sheet.write("A1","opentrons-tiprack-300ul")

        writer.save()
        print("Successfully completed experiment plan construction!")
        return

    def calculate_total_reagent_volumes(self):
        """
        Prints total volumes needed in microliters
        Assumes well volume of 3 mL
        """
        # Very important to pass by value and not reference here
        self.volumes = self.design.copy()
        for k in self.stock_conc.keys():
            self.volumes[k] = self.design[k]*self.well_volume/self.stock_conc[k]*1e3 # mL

        try:
            total_vols = self.volumes.sum(axis=1)
            assert(np.all(np.logical_and(total_vols>=0,
                total_vols<1e3*self.well_volume)))
        except:
            print(self.design)
            print(self.volumes.sum())
            print("Error: Max well volume exceeded: "+str(max(self.volumes.sum(axis=1))))
            quit()
        try:
            assert(not np.any(np.logical_and(self.volumes.values>0,
                      self.volumes.values<0.001)))
        except:
            print(self.volumes)
            print(self.volumes.sum())
            print("Error: Sub uL volume required")
            quit()

        r_vols = self.volumes.sum()

        self.working_volumes = 5*np.ceil((r_vols+2)/5)
        return


    def handle_large_stock_volumes(self):
        large_volumes = self.volumes.columns.values[self.volumes.sum()>45]
        idx = np.arange(len(self.volumes.index))
        np.random.shuffle(idx)
        for lv in large_volumes:
            n_stocks = np.ceil(self.volumes[lv].sum()/45)
            new_rows = np.array_split(idx,n_stocks)
            for i,this_idx in enumerate(new_rows):
                name = lv+"_"+str(i)
                self.volumes[name] = np.zeros(len(self.volumes.index))
                self.design[name] = np.zeros(len(self.design.index))
                self.volumes.loc[this_idx,name] = self.volumes.loc[this_idx,lv]
                self.design.loc[this_idx,name] = self.design.loc[this_idx,lv]
                self.stock_conc[name] = self.stock_conc[lv]
            self.volumes.drop(lv,axis=1,inplace=True)
            self.design.drop(lv,axis=1,inplace=True)
            del self.stock_conc[lv]
        self.calculate_total_reagent_volumes()
        return

    def configure_stocks_setup(self):
        # Calculate the most space efficient setup
        c6x50 = 0
        c4x50_6x15 = 0
        c15x15 = 0

        if not np.all(self.working_volumes <= 45):
            self.handle_large_stock_volumes()

        # check that we've fixed the issue
        assert(np.all(self.working_volumes <= 45))
        n_15 = sum(self.working_volumes <= 10)
        n_50 = len(self.working_volumes)-n_15
        while n_50 > 4:
            c6x50 += 1
            n_50 -= 6
        if n_50 > 0:
            c4x50_6x15 += 1
            n_50 -= 4
            n_15 -= 6
        n_15 += n_50 # accounting for # of 50's that can be filled
        while n_15 > 0:
            c15x15 += 1
            n_15 -= 15

        # Assign 50s to 50s and 15s to everything else
        comps = self.working_volumes.index.values
        wv = self.working_volumes.values
        comps_50 = comps[self.working_volumes>10]
        comps_15 = comps[self.working_volumes<=10]
        wv_50 = wv[self.working_volumes>10]
        wv_15 = wv[self.working_volumes<=10]
        stocks = []

        while c6x50 > 0:
            this_wells = ["A1","A2","A3",
                     "B1","B2","B3"]
            if len(comps_50)<6:
                this_comps = comps_50
                this_wv = wv_50
            else:
                this_comps = comps_50[:6]
                comps_50 = comps_50[6:]
                this_wv = wv_50[:6]
                wv_50 = wv_50[6:]
            this_conc = [self.stock_conc[s] for s in this_comps]
            this_df = pd.DataFrame.from_dict({"Location":this_wells,
                                    "Component":this_comps,
                                    "Concentration (g/L)":this_conc,
                                    "Working Volume (mL)":this_wv},
                                    orient='index').T
            stocks.append(("opentrons_6_tuberack_falcon_50ml_conical",
                           this_df))
            c6x50 -= 1

        if c4x50_6x15 > 0:
            this_wells = ["A1","A2","A3","A4",
                          "B1","B2","B3","B4",
                          "C1","C2"]
            this_comps = [None for i in range(len(this_wells))]
            this_wv = [None for i in range(len(this_wells))]
            well_50 = [2,3,6,7]
            counter_50 = 0
            well_15 = [0,1,4,5,8,9]
            counter_15 = 0
            while len(comps_50)>0:
                this_comps[well_50[counter_50]] = comps_50[0]
                this_wv[well_50[counter_50]] = wv_50[0]
                if len(comps_50)>1:
                    comps_50 = comps_50[1:]
                    wv_50 = wv_50[1:]
                else:
                    comps_50 = []
                counter_50 += 1
            while((len(comps_15)>0)&(counter_15<6)):
                this_comps[well_15[counter_15]] = comps_15[0]
                this_wv[well_15[counter_15]] = wv_15[0]
                if len(comps_15)>1:
                    comps_15 = comps_15[1:]
                    wv_15 = wv_15[1:]
                else:
                    comps_15 = []
                counter_15 += 1
            while((len(comps_15)>0)&(counter_50<4)):
                this_comps[well_50[counter_50]] = comps_15[0]
                this_wv[well_50[counter_50]] = wv_15[0]
                if len(comps_15)>1:
                    comps_15 = comps_15[1:]
                    wv_15 = wv_15[1:]
                else:
                    comps_15 = []
                counter_50 += 1
            this_conc = [] 
            for s in this_comps:
                if s:
                    this_conc.append(self.stock_conc[s])
                else:
                    this_conc.append(None)
            this_df = pd.DataFrame.from_dict({"Location":this_wells,
                                    "Component":this_comps,
                                    "Concentration (g/L)":this_conc,
                                    "Working Volume (mL)":this_wv},
                                    orient='index').T
            stocks.append(("opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical",
                           this_df))
            c6x50 -= 1
                
        while c15x15 > 0:
            this_wells = ["A1","A2","A3","A4","A5",
                          "B1","B2","B3","B4","B5",
                          "C1","C2","C3","C4","C5"]
            if len(comps_15)<15:
                this_comps = comps_15
                this_wv = wv_15
            else:
                this_comps = comps_15[:15]
                comps_15 = comps_15[15:]
                this_wv = wv_15[:15]
                wv_15 = wv_15[15:]
            this_conc = [self.stock_conc[s] for s in this_comps]
            this_df = pd.DataFrame.from_dict({"Location":this_wells,
                                    "Component":this_comps,
                                    "Concentration (g/L)":this_conc,
                                    "Working Volume (mL)":this_wv},
                                    orient='index').T
            stocks.append(("opentrons_15_tuberack_falcon_15ml_conical",
                           this_df))
            c15x15 -= 1
            break

        self.stocks = stocks
        return

