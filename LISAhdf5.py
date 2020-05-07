"""
    Copyright (C) 2017 Stas Babak, Antoine Petiteau for the LDC team

    This file is part of LISA Data Challenge.

    LISA Data Challenge is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

##################################################
#                                                #
#     Structure for the LISA hdf5 file/object    #
#                 version 1.1                    #
#                                                #
#         S. Babak, A. Petiteau, ...             #
#      for the LISA Data Challenge team          #
#                                                #
##################################################





import numpy as np
import h5py as h5
import os, sys, re

import LISAConstants as LC
from LISACommonFunctions import *


author_default = "Stas Babak, Antoine Petiteau (2017)"

def Str(a):
    if type(a)==bytes:
        return a.decode('utf8')
    if type(a)!=str:
        return str(a)
    return a

class ParsUnits():
    """
    This class defines the small object to manage parameters and units.
    """

    def __init__(self, pars_i=None, units_i=None, name='', value=0., unit=''):
        """
        Initalize `ParsUnits`
        @param pars_i is an optional dictionary of parameters
        @param units_i is an optional dictionary of units (same size as pars_i)
        @param name is an optional name
        @param value is an optional value
        @param unit is an optional unit
        """
        self.pars = {}
        self.units = {}
        if pars_i is not None:
            self.addDict(pars_i,units_i)
        elif name!='':
            self.addPar(name,value,unit)

    def __del__(self):
        # so far is empty, see what we need to add here
        pass

    def display(self,ReturnStr=False):
        """
        Display all parameters
        @param ReturnStr is true to return string instead of display
        """
        r = ""
        for i,k in enumerate(self.pars):
            r = r + "\t"+str(k)+" "
            if type(self.pars[k])==str:
                r = r + self.pars[k]
            else:
                r = r + str(self.pars[k])
            r = r + " ["+Str(self.units[k])+"]\n"
        if ReturnStr:
            return r
        else:
            print(r)

    def addPar(self,name,value,unit):
        """
        Add parameter, its value and unit
        @param name is the name of the parameter
        @param value is the value of the parameter
        @param unit is the unit of the parameter
        """
        self.pars.update({name : value})
        self.units.update({name : unit})

    def addDict(self,pars_i,units_i):
        """
        Add dictionnary and its unit
        @param pars_i is a dictionary of parameters
        @param units_i is a dictionary of units(same size as pars_i)
        """
        if len(pars_i)==len(units_i):
            self.pars.update(pars_i)
            self.units.update(units_i)
        else:
            raise Exception('addDict : parameters and units should have the same number of elements.')

    def get(self,parName):
        """
        Get parameter value
        @param parName parameter name
        @return value
        """
        ### TODO : Use proper error system
        if parName not in self.pars :
            print("WARNING: ParsUnits.get :",parName,"is not a parameter name (",self.pars,")")
            return None
        else:
            return self.pars[parName]

    def getConvert(self,parName,conversion,requiredUnit):
        """
        Get parameter value for distance parameters
        @param parName is the parameter name
        @param conversion is the conversion dictionary:
            + LC.convT for mass, time and distance (everything in time)
            + LC.convMass for mass only
            + LC.convTime for time only
            + LC.convDistance for distance only
        @param requiredUnit is the required unit [default:s]
        @return value
        """
        ### TODO : Use proper error system
        v = self.get(parName)
        if type(v)!=type(np.zeros(10)) and v == None:
            return None
        else:
            uV = self.units[parName]
            uV = uV.lower()
            requiredUnit = requiredUnit.lower()
            if requiredUnit not in conversion:
                print("WARNING: ParsUnits.getConvert : parameter unit",requiredUnit,"is not in the conversion list (",conversion,")")
                return None
            if uV not in conversion:
                print("WARNING: ParsUnits.getConvert : required unit",uV,"is not in the conversion list (",conversion,")")
                return None
            return v * ( conversion[uV] / conversion[requiredUnit] )

class LISAhdf5():
    """
    This class defines the LISA hdf5 file/object which stores information about
    LISA: GW sources, instrument, ...
    Here we attempt to make it rather flexible.
    Basically it is IO between user and hdf5.
    """

    h5file = None

    def __init__(self, filename, author=""):
        """
        Create hdf5 file and create the main group
        @param filename: name of hdf5 file
        """
        self.filename = filename
        self.h5file = h5.File(self.filename, mode='a')
        self.MainGroupName = 'H5LISA'
        self.MainGroup = self.h5file.require_group(self.MainGroupName)
        if not 'Author' in self.MainGroup:
            if author!="":
                self.MainGroup.create_dataset('Author', data=author)
            else:
                self.MainGroup.create_dataset('Author', data=author_default)
        self.h5file.close()


    def __del__(self):
        # so far is empty, see what we need to add here
        pass


    def display(self, ReturnStr=False, Print=True, DispPath=False):
        """
        Display the content of LISA hdf5
        @param Print is True to display [default: True]
        @param ReturnStr is True to ReturnStr [default: False]
        @param DispPath is True to display the path of each element[default: False]
        """
        list_of_names = []
        self.h5file = h5.File(self.filename,'r')
        self.h5file.visit(list_of_names.append)
        g0 = self.h5file.get('H5LISA')
        d = [] # list of list of line to display and path
        mc = 0 # max length line
        for x in list_of_names:
            r = ''
            ntab = len(x.split('/'))
            r = r + ntab*'  ' + '- '+ os.path.basename(x)
            p = self.h5file.get(x)
            units = ""
            if 'Units' in list(p.attrs.keys()):
                u = Str(p.attrs['Units'])
                units = " [" + u +"]"
            if type(g0)!=type(p):
                ## If not a group
                r = r + ' : '
                if len(p.shape)==0:
                    ## If standard parameters
                    v = p.value
                    if type(v)==bytes:
                        v = Str(v)
                    if type(v)==str:
                        lines = v.split('\n')
                        if len(lines) == 1:
                            r = r + v + " " + units
                        else:
                            r = r + '[multiple lines]' + units
                            #for line in lines:
                            #    r = r + (ntab+1)*'  ' + line + '\n'
                    else:
                        r = r + str(v) + units
                else:
                    ## If array
                    if len(p.shape)==1 and p.shape[0]<=20 :
                        r = r + '['
                        for i in range(p.shape[0]):
                            r = r + ' ' + str(p[i])
                        r = r + ' ]' + units
                    else:
                        r = r + '[array '
                        for i in range(len(p.shape)):
                            if i!=0:
                                r = r + ' x '
                            r = r + str(p.shape[i])
                        r = r + ' ]' + units
            d.append([r,x])
            mc = max(mc,len(r))
        rs = ''
        for x in d:
            if DispPath:
                rs = rs + x[0]
                for i in range(len(x[0]),mc+3):
                    rs = rs + " "
                rs = rs + x[1] + os.linesep
            else:
                rs = rs + x[0] + os.linesep

        self.h5file.close()

        if Print :
            print(rs)
        if ReturnStr:
            return rs


    def compare(self, ref, verbose=False):
        """
        Compare with another LISAhdf5
        @param ref pther LISAhdf5 object
        """
        ### Check
        r_loc = self.display(ReturnStr=True,Print=False,DispPath=True)
        r_ref = ref.display(ReturnStr=True,Print=False,DispPath=True)
        l_loc = r_loc.split('\n')
        l_ref = r_ref.split('\n')
        ## Compare number of lines
        if len(l_loc) != len(l_ref):
            if verbose:
                print("LISAhdf5.compare : different number object: ",len(l_loc),"vs",len(l_ref))
            return False
        ## Compare lines
        DiffLine = False
        for i in range(len(l_loc)):
            ll = l_loc[i].replace(" ", "")
            lr = l_ref[i].replace(" ", "")
            print(">>>"+ll+"<<<")
            print(">>>"+lr+"<<<")
            if ll != lr and (re.search('H5LISA/History',ll)==None):
                if verbose:
                    print("LISAhdf5.compare : difference in line ",i,": ",len(ll),"vs",len(lr))
                    print(">>>"+ll+"<<<")
                    print(">>>"+lr+"<<<")
                DiffLine = True
        if DiffLine:
            return False
        ## Compare data

        loc_list = []
        self.h5file = h5.File(self.filename,'r')
        self.h5file.visit(loc_list.append)
        g0 = self.h5file.get('H5LISA')
        for x in loc_list:
            if x[:14]!='H5LISA/History':
                p = self.h5file.get(x)
                if type(g0)!=type(p) and len(p.shape)!=0:
                    d_l = np.copy(p.value)
                    d_r = ref.get(x)
                    ib = np.where(np.isclose(d_l,d_r)==False)[0]
                    if len(ib)!=0:
                        if verbose:
                            print("LISAhdf5.compare : difference for array ",x,"at index",ib)
                        self.h5file.close()
                        return False
        self.h5file.close()
        ## All comparison are ok
        return True


    def get(self, path, returnUnit=False):
        """
        Return the parameter corresponding to the path
        @param path is the path to the parameter within the hdf5
        @param ReturnUnit is True to return a ParsUnits object including the Unit
        """
        self.h5file = h5.File(self.filename, 'r')
        if path in self.h5file:
            v = self.h5file[path].value
            if type(v)==bytes or type(v)==str :
                v = Str(v)
            if returnUnit:
                name = os.path.basename(path)
                u = Str(self.h5file[path].attrs['Units'])
                self.h5file.close()
                return ParsUnits({name : v}, {name : u})
            else:
                self.h5file.close()
                return v
        else:
            self.h5file.close()
            self.display()
            print("ERROR: LISAhdf5:get : No parameter in the hdf5 file ("+self.filename+") corresponding to the path ("+path+")")
            return None


    def set(self, path, newPar, newUnit='NotGiven'):
        """
        Set the parameter corresponding to the path
        @param path is the pathh to the parameter within the hdf5
        @param newPar is the new parameter
        @param newUnit if specified, the new unit
        """
        self.h5file = h5.File(self.filename, 'a')

        if path in self.h5file:
            data = self.h5file[path]
            data[...] =  newPar
            if newUnit!='NotGiven':
                data.attrs["Units"] = newUnit
        else:
            self.display(DispPath=True)
            print("ERROR: LISAhdf5:set : No parameter in the hdf5 file ("+self.filename+") corresponding to the path ("+path+")")
        self.h5file.close()


    def delete(self, path):
        """
        Delete the parameter corresponding to the path
        @param path is the pathh to the parameter within the hdf5
        """
        self.h5file = h5.File(self.filename, 'a')

        if path in self.h5file:
            self.h5file.__delitem__(path)
        else:
            self.display(DispPath=True)
            print("ERROR: LISAhdf5:set : No parameter in the hdf5 file ("+self.filename+") corresponding to the path ("+path+")")
        self.h5file.close()


    def addSource(self, sourceName, parameters, overwrite=False, hphcData=None, XYZdata=None):
        """
        Add one single GW source
        @param sourceName name (string) to be given to this source
        @param parameters is a set of parameters with their units
        @param hpcFile file containing time, h_plus, h_cross, or freq, h_p, h_c
        """

        ### Open hdf5 file
        self.h5file = h5.File(self.filename, 'a')

        ### Check if pathGroup exist, if not create it
        pathGroupSrc = "/"+self.MainGroupName+"/GWSources/"+sourceName

        try:
            self.SrcGroup = self.h5file.create_group(pathGroupSrc)
        except:
            print("WARNING in LISAhdf5:addSource :")
            print("\tLooks like this source (",sourceName,") exist or name is used twice")
            if (overwrite):
                print("\t=> proceed with overwriting it")
                self.SrcGroup = self.h5file.require_group(pathGroupSrc)
            else:
                print("ERROR: LISAhdf5:addSource : Cannot overwrite an existing source ",sourceName,"!")
                sys.exit(1)

        # TODO Do we want a sanity check here that parameters and their values make sense?
        # TODO need to check that the size of units and par dictionary is the same
        ks =  list(self.SrcGroup.keys())
        for i, key in enumerate(parameters.pars.keys()):

            if (key in ks):
                if (overwrite):
                    #print self.SrcGroup[key]
                    dataset = self.SrcGroup[key]
                    dataset[...] = parameters.pars[key]
            else:
                dataset = self.SrcGroup.create_dataset(key, data=parameters.pars[key])
            dataset.attrs["Units"] = parameters.units[key]


        if (np.any(hphcData)):
            if ('hphcData' not in ks):
                HpHcData = self.SrcGroup.create_dataset('hphcData', np.shape(hphcData), dtype=h5.h5t.IEEE_F64LE)
                HpHcData[...] = hphcData
                HpHcData.attrs['h+, hx struct'] = 't, h+, hx'
            elif (overwrite):
                    HpHcData = self.SrcGroup['hphcData']
                    HpHcData[...] = hphcData
                    HpHcData.attrs['h+, hx struct'] = 't, h+, hx'
            else:
                print("Nothing to do for hp hc data: data already exist and not overwrite.")

        #if (np.any(XYZdata['X'])):
        if (np.any(XYZdata)):
            if ('XYZdata' not in ks):
                tdiXYZdata = self.SrcGroup.create_dataset('XYZdata', np.shape(XYZdata), dtype=h5.h5t.IEEE_F64LE)
                tdiXYZdata[...] = XYZdata
                tdiXYZdata.attrs['X, Y, Z struct'] = 't, X, Y, Z'
            elif (overwrite):
                    tdiXYZdata = self.SrcGroup['XYZdata']
                    tdiXYZdata[...] = XYZdata
                    tdiXYZdata.attrs['X, Y, Z struct'] = 't, X, Y, Z'
            else:
                print("Nothing to do for XYZ data: data already exist and not overwrite.")


        """
        if (TDIfile != None):
            datasetTDI = np.genfromtxt(TDIfile, skip_header=1)
            if ('TDIData' not in ks):
                TDIdata = self.SrcGroup.create_dataset('TDIdata', np.shape(datasetTDI), dtype=h5.h5t.IEEE_F64LE)
                TDIdata[...] = datasetTDI
                TDIdata.attrs['TDI struct'] = 't,X,Y,Z'
                TDIdata.attrs['TDY type'] = 'frac. freq.'
            elif (overwrite):
                TDIdata = self.SrcGroup['TDIdata']
                TDIdata[...] = datasetTDI
                TDIdata.attrs['TDI struct'] = 't,X,Y,Z'
                TDIdata.attrs['TDY type'] = 'frac. freq.'
            else:
                print "Nothing to do for hp hc data: data already exist and not overwrite."
        """
        self.h5file.close()


    def getSourcesNum(self):
        """
        Return the number of GW sources
        """
        self.h5file = h5.File(self.filename, 'r')
        gpgw = self.MainGroupName+"/GWSources"
        try:
            gp = self.h5file.get(gpgw)
            Nsrc = len(gp)
            self.h5file.close()
            return (Nsrc)
        except:
            print("WARNING in LISAhdf5:getSourcesName : No GWs : ",gpgw)
            self.h5file.close()
            return 0


    def getSourcesName(self):
        """
        Return the list of names of GW sources
        """
        self.h5file = h5.File(self.filename, 'r')
        gpgw = self.MainGroupName+"/GWSources"
        try:
            gp = self.h5file.get(gpgw)
            keys = list(gp.keys())
            self.h5file.close()
            return keys
        except:
            print("WARNING in LISAhdf5:getSourcesName : No GWs : ",gpgw)
            self.h5file.close()
            return []


    def getSourceParameters(self, SrcName):
        """
        Return the parameters corresponding to the source #indexSource
        @param SrcName od is the source index starting at 0
        @return ParsUnits object containing the parameters
        """
        SrcNames = self.getSourcesName()
        self.h5file = h5.File(self.filename, 'r')
        gp = self.h5file.get(self.MainGroupName+"/GWSources")
        params = ParsUnits()
        if SrcName not in SrcNames :
            print("WARNING: LISAhdf5.getSource: ",SrcName, "is not in source names", SrcNames, "!")
        else:
            src = gp.get(SrcName)
            for ky in list(src.keys()):
                if ky!='hphcData' and ky!='XYZdata' :
                    p = src.get(ky)
                    params.addPar(ky,p.value,Str(p.attrs['Units']))
        self.h5file.close()
        return params


    def getSourceHpHc(self, SrcName):
        """
        Return the t,h+,hx numpy array corresponding to the source #indexSource
        @param SrcName od is the source index starting at 0
        @return numpy array t=[:,0], h+=[:,1], hx=[:,2]
        """
        SrcNames = self.getSourcesName()
        self.h5file = h5.File(self.filename, 'r')
        gp = self.h5file.get(self.MainGroupName+"/GWSources")
        params = ParsUnits()
        if SrcName not in SrcNames :
            print("WARNING: LISAhdf5.getSource: ",SrcName, "is not in source names", SrcNames, "!")
        else:
            src = gp.get(SrcName)
            if 'hphcData' in list(src.keys()):
                d = np.copy(src.get('hphcData'))
                self.h5file.close()
                return d
            else:
                print("ERROR in LISAhdf5:getSourceHpHc: there is no h+, hc data for the source",SrcName)
                sys.exit(1)
        self.h5file.close()


    def addLISADataSource(self, Name, Model, parameters, overwrite=False, data=None):
        """
        Add LISA data source. A data source is an object producing data (timeseries)
        by itself (without any input). Examples: noise, orbits, ...
        @param Name name (string) to be given to this data source
        @param Model model of data source (ex: WhiteFrequencyNoise, Eccentric, ...)
        @param parameters is a set of parameters with their units
        @param overwrite is true if we want to replace existing DataSource withthe same name
        @param data are potential data associated to the DataSource
        """
        ### Open hdf5 file
        self.h5file = h5.File(self.filename, 'a')

        ### Check if pathGroup DataSources exist, if not create it
        pathGroupDataSource = "/"+self.MainGroupName+"/Observatory/DataSources/"+Name

        try:
            groupDataSource = self.h5file.create_group(pathGroupDataSource)
        except:
            print("WARNING in LISAhdf5:addLISADataSource :")
            print("\tLooks like this data source",Name,"exist or name is used twice")
            if (overwrite):
                print("\t=> proceed with overwriting it")
                groupDataSource = self.h5file.require_group(pathGroupDataSource)
            else:
                print("ERROR: LISAhdf5:addLISADataSource : Cannot overwrite an existing data source ",Name,"!")
                sys.exit(1)

        ### Add the model
        if not 'Model' in groupDataSource:
            groupDataSource.create_dataset('Model', data=Model)
        elif overwrite:
            groupDataSource['Model'][...] = Model


        ### Add parameters
        ks =  list(groupDataSource.keys())
        for i, key in enumerate(parameters.pars.keys()):

            if (key in ks):
                if (overwrite):
                    dataset = groupDataSource[key]
                    dataset[...] = parameters.pars[key]
            else:
                dataset = groupDataSource.create_dataset(key, data=parameters.pars[key])
            dataset.attrs["Units"] = parameters.units[key]

        if ( np.any(data) ):
            if ('data' not in ks):
                noiseData = groupDataSource.create_dataset('data', np.shape(data), dtype=h5.h5t.IEEE_F64LE)
                noiseData[...] = data
                noiseData.attrs['struct'] = 't, value'
            elif (overwrite):
                    noiseData = groupDataSource['data']
                    noiseData[...] = data
                    noiseData.attrs['struct'] = 't, value'
            else:
                print("Nothing to do for data: data already exist and not overwrite.")

        self.h5file.close()


    def addLISADataModel(self, Name, Model, parameters, InputsData={}, OutputsData={}, overwrite=False, data=None):
        """
        Add LISA data model. A data model is an object producing data (timeseries)
        by itself (without any input). Examples: noise, orbits, ...
        @param Name name (string) to be given to this data source
        @param Model model of data source (ex: WhiteFrequencyNoise, Eccentric, ...)
        @param parDict dictionary of parameter name and its value
        @param units supplies units of parameters which will be added as attributes
        @param InputsData is dictionnary of input data
        @param OutputsData is dictionnary of output data
        @param overwrite is true if we want to replace existing DataSource withthe same name
        @param data are potential data associated to the DataSource
        """
        ### Open hdf5 file
        self.h5file = h5.File(self.filename, 'a')

        ### Check if pathGroup DataSources exist, if not create it
        pathGroupDataModel = "/"+self.MainGroupName+"/Observatory/DataModels/"+Name

        try:
            groupDataModel = self.h5file.create_group(pathGroupDataModel)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print("\tLooks like this data model",Name,"exist or name is used twice")
            if (overwrite):
                print("\t=> proceed with overwriting it")
                groupDataModel = self.h5file.require_group(pathGroupDataModel)
            else:
                print("ERROR: LISAhdf5:addLISADataModel : Cannot overwrite an existing data model ",Name,"!")
                sys.exit(1)

        ### Add the model
        if not 'Model' in groupDataModel:
            groupDataModel.create_dataset('Model', data=Model)
        elif overwrite:
            groupDataModel['Model'][...] = Model

        ### Add parameters
        ks =  list(groupDataModel.keys())
        for i, key in enumerate(parameters.pars.keys()):
            if (key in ks):
                if (overwrite):
                    dataset = groupDataModel[key]
                    dataset[...] = parameters.pars[key]
            else:
                dataset = groupDataModel.create_dataset(key, data=parameters.pars[key])
            dataset.attrs["Units"] = parameters.units[key]

        ### Add group for input data (=node)
        pathGroupInputs = pathGroupDataModel+'/Inputs'
        try:
            groupInputs = self.h5file.create_group(pathGroupInputs)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print("\tLooks like this group Inputs",Name,"exist or name is used twice")
            if (overwrite):
                print("\t=> proceed with overwriting it")
                groupInputs = self.h5file.require_group(pathGroupInputs)
            else:
                print("ERROR: LISAhdf5:addLISADataModel : Cannot overwrite an existing input ",pathGroupInputs,"!")
                sys.exit(1)
        for i, key in enumerate(InputsData.keys()):
            if (key in ks):
                if (overwrite):
                    dataset = groupInputs[key]
                    dataset[...] = InputsData[key]
            else:
                dataset = groupInputs.create_dataset(key, data=InputsData[key])

        ### Add group for output data (=node)
        pathGroupOutputs = pathGroupDataModel+'/Outputs'
        try:
            groupOutputs = self.h5file.create_group(pathGroupOutputs)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print("\tLooks like this group Outputs",Name,"exist or name is used twice")
            if (overwrite):
                print("\t=> proceed with overwriting it")
                groupOutputs = self.h5file.require_group(pathGroupOutputs)
            else:
                print("ERROR: LISAhdf5:addLISADataModel : Cannot overwrite an existing ouput ",pathGroupOutputs,"!")
                sys.exit(1)
        for i, key in enumerate(OutputsData.keys()):
            if (key in ks):
                if (overwrite):
                    dataset = groupOutputs[key]
                    dataset[...] = OutputsData[key]
            else:
                dataset = groupOutputs.create_dataset(key, data=OutputsData[key])

        self.h5file.close()


    def addPreProcess(self, parameters, overwrite=False, TDIdata=None):
        """
        Add PreProcess (TDI)
        @param parameters is a set of parameters with their units
        @param overwrite is true if we want to replace existing PreProcess with the same name
        """
        ### Open hdf5 file
        self.h5file = h5.File(self.filename, 'a')

        ### Check if pathGroup DataSources exist, if not create it
        pathGroupPreProcess = "/"+self.MainGroupName+"/PreProcess"

        groupPreProcess = self.h5file.require_group(pathGroupPreProcess)

        ### Add parameters
        ks =  list(groupPreProcess.keys())
        for i, key in enumerate(parameters.pars.keys()):
            if (key in ks):
                if (overwrite):
                    dataset = groupPreProcess[key]
                    dataset[...] = parameters.pars[key]
            else:
                dataset = groupPreProcess.create_dataset(key, data=parameters.pars[key])
            dataset.attrs["Units"] = parameters.units[key]

        if ( np.any(TDIdata) ):
            if ('TDIdata' not in ks):
                tdiDataset = groupPreProcess.create_dataset('TDIdata', np.shape(TDIdata), dtype=h5.h5t.IEEE_F64LE)
                tdiDataset[...] = TDIdata
                tdiDataset.attrs['struct'] = 't, value'
            elif (overwrite):
                    tdiDataset = groupPreProcess['TDIdata']
                    tdiDataset[...] = TDIdata
                    tdiDataset.attrs['struct'] = 't, value'
            else:
                print("Nothing to do for TDIdata: TDIdata already exist and not overwrite.")

        self.h5file.close()


    def getPreProcessTDI(self):
        """
        Return the numpy matrix of TDI data from PreProcess
        """
        self.h5file = h5.File(self.filename, 'r')
        try:
            list_of_names = []
            self.h5file.visit(list_of_names.append)
            if 'H5LISA/PreProcess/TDIdata' in list_of_names :
                d = np.copy(self.h5file['H5LISA/PreProcess/TDIdata'].value)
                self.h5file.close()
                return d
            else:
                print("ERROR: LISAhdf5:getPreProcessTDI : No TDIdata (H5LISA/PreProcess/TDIdata) in this LISAhdf5 (see content below)")
                sys.exit(1)
        except:
            self.display()
            self.h5file.close()
            print("ERROR: LISAhdf5:getPreProcessTDI : No TDIdata (H5LISA/PreProcess/TDIdata) in this LISAhdf5 (see content below)")
            sys.exit(1)


    def addSimulation(self, parameters, overwrite=False):
        """
        Add Simulation
        @param parameters is a set of parameters with their units
        @param overwrite is true if we want to replace existing DataSource withthe same name
        """
        ### Open hdf5 file
        self.h5file = h5.File(self.filename, 'a')

        ### Check if pathGroup DataSources exist, if not create it
        pathGroupSimulation = "/"+self.MainGroupName+"/Simulation"

        groupSimulation = self.h5file.require_group(pathGroupSimulation)

        ### Add parameters
        ks =  list(groupSimulation.keys())
        for i, key in enumerate(parameters.pars.keys()):
            if (key in ks):
                dataset = groupSimulation[key]
                if (overwrite):
                    dataset[...] = parameters.pars[key]
            else:
                dataset = groupSimulation.create_dataset(key, data=parameters.pars[key])
            dataset.attrs["Units"] = parameters.units[key]

        self.h5file.close()


    def addUserRequest(self, Name, text, overwrite=True):
        """
        Add the user request for creating the simulated data
        @param Name name (string) to be given to this data source
        @param text User request in form of the string
        """
        ### Format text
        strtext = ""
        if type(text)==str:
            strtext = text
        elif type(text)==list:
            for x in text:
                strtext = strtext + x
        ### Open hdf5 file
        self.h5file = h5.File(self.filename, 'a')

        ### Check if pathGroup DataSources exist, if not create it
        pathGroupUserRequest = "/"+self.MainGroupName+"/UserRequest/"+Name

        try:
            groupUserRequest = self.h5file.create_group(pathGroupUserRequest)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print("\tLooks like this user request ",Name,"exist or name is used twice")
            if (overwrite):
                print("\t=> proceed with overwriting it")
                groupUserRequest = self.h5file.require_group(pathGroupUserRequest)
            else:
                print("ERROR: LISAhdf5:addUserRequest : Cannot overwrite an existing user request",pathGroupUserRequest,"!")
                sys.exit(1)
        ks =  list(groupUserRequest.keys())
        if ("Request" in ks):
            dataset = groupUserRequest["Request"]
            dataset[...] = strtext
        else:
            dataset = groupUserRequest.create_dataset("Request", data=strtext)

        self.h5file.close()


    def addHistory(self, filepath, options, args):
        """
        Add running script informations to the history
        @param filepath is file path : os.path.realpath(__file__)
        @param options is the dictionary of options : vars(options)
        @param args is the list of arguments : args
        """
        NewHist = ''
        try:
            NewHist =  GetStrCodeGitCmd(filepath,options,args)
        except:
            print("WARNING in LISAhdf5.addHistory: cannot recover script informations")
        self.h5file = h5.File(self.filename, 'a')
        pathGroupHistory = "/"+self.MainGroupName+"/History"
        try:
            groupHistory = self.h5file.create_group(pathGroupHistory)
        except:
            groupHistory = self.h5file.require_group(pathGroupHistory)
        ks =  list(groupHistory.keys())
        NewHistName = '%04d'%(len(ks))
        dataset = groupHistory.create_dataset(NewHistName, data=NewHist)
        self.h5file.close()


    def readSources(self):

        self.h5file = h5.File(self.filename, 'r')
        gp = self.h5file.get(self.MainGroupName+"/GWSources")
        Nsrc = len(gp)
        print("Number of sources",  Nsrc)

        keys = list(gp.keys())
        Nsrc = len(keys)
        print(keys)
        print(Nsrc)
        srcs = []
        for ky in keys:
            print("Reading source:", ky)
            src = gp.get(ky)
            nms =  list(src.keys())
            vls = list(src.values())
            #print nms

            source = {}
            att = {}
            source["Source"] = ky
            for i, nm in enumerate(nms):
                #print nm, src.get(nm)[...]
                stas_tmp = np.array(src.get(nm)[...])

                source[nm] = stas_tmp
                if (nm != 'TDIdata'):
                   #att[nm] = src.attrs.get('Units')
                   rec = src.get(nm)
                   #print "attr:", rec.attrs.get('Units')
                   att[nm] = rec.attrs.get('Units')

                #print "stas", stas_tmp

            srcs.append([source, att])
        self.h5file.close()
        return srcs


    def readDataModel(self):

        self.h5file = h5.File(self.filename, 'r')
        gp = self.h5file.get(self.MainGroupName+"/Observatory/DataSources/")

        Ndm = len(gp)
        keys = list(gp.keys())
        Ndm = len(keys)
        print(keys)
        print(Ndm)
        srcs = []
        for ky in keys:
            print("Reading source:", ky)
            src = gp.get(ky)
            nms =  list(src.keys())
            vls = list(src.values())
            #print nms, vls
            source = {}
            att = {}
            source["DataModel"] = ky
            for i, nm in enumerate(nms):
                #print nm, src.get(nm)[...]
                stas_tmp = np.array(src.get(nm)[...])

                source[nm] = stas_tmp
                if (nm != 'data'):
                   #att[nm] = src.attrs.get('Units')
                   rec = src.get(nm)
                   #print "attr:", rec.attrs.get('Units')
                   att[nm] = rec.attrs.get('Units')

                #print "stas", stas_tmp

            srcs.append([source, att])

        self.h5file.close()
        return(srcs)





        """self.h5file = h5.File(self.filename, 'r')
        key =  self.h5file.keys()
        if (key[0] != 'BBH_GW_signals'):
            print "Cannot find the group: BBH_GW_signals "
            self.h5file.close()
            sys.exit(1)
        gp = self.h5file.get('BBH_GW_signals')
        keys = gp.keys()
        Nsrc = len(keys)
        if ('Author' in keys):
            Nsrc = Nsrc -1
        print "found", Nsrc, "sources"
        srcs = []
        for ky in keys:
            if (ky != "Author"):
                #print ky
                src = gp.get(ky)
                nms =  src.keys()
                vls = src.values()
                #print nms

                source = {}
                att = {}
                for i, nm in enumerate(nms):
                    #print nm, src.get(nm)[...]
                    stas_tmp = np.array(src.get(nm)[...])

                    source[nm] = stas_tmp
                    if (nm != 'TDIdata'):
                       #att[nm] = src.attrs.get('Units')
                       rec = src.get(nm)
                       #print "attr:", rec.attrs.get('Units')
                       att[nm] = rec.attrs.get('Units')

                    #print "stas", stas_tmp

                srcs.append([source, att])
        #print srcs
        #print len(srcs)
        s1 = srcs[0]
        #print s1.keys()
        #print s1.values()
        self.h5file.close()
        return srcs
        """
