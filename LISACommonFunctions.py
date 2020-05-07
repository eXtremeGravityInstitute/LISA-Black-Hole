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
#            LISA common functions               #
#                 version 1.0                    #
#                                                #
#         A. Petiteau, ...                       #
#      for the LISA Data Challenge team          #
#                                                #
##################################################


import os, sys, re
import numpy as np
import subprocess


def run(command, disp=False, NoExit=False):
    """
    Run system command
    @param command is a string to run as a system command
    @param disp is true to display the command
    @param is true to continue if the command failed
    """
    commandline = command % globals()
    if disp :
        print("----> %s" % commandline)
    try:
        assert(os.system(commandline) == 0)
    except:
        print('Script %s failed at command "%s".' % (sys.argv[0],commandline))
        if not NoExit :
            sys.exit(1)


def makefromtemplate(output,template,keyChar,**kwargs):
    """
    Create an output file identical to the template file with the string
    between the key character KeyChar and corresponding to a keyword in kwargs
    replaced by the argument of the keyword
    @param output is the output file
    @param template is the template file
    @param keyChar is the key character surrounding keywords
    @param kwargs are keywords and corresponding arguments
    """
    fi = open(template,'r')
    fo = open(output,'w')
    for line in fi:
        repline = line
        for kw in kwargs:
            repline = re.sub(keyChar + kw + keyChar,str(kwargs[kw]),repline)
        print(repline, end=' ', file=fo)

def LoadFileNpStruct(FileName):
    """
    Load txt file of data and return a structured numpy array
    @param filename with one column per record and title of the record on the first line "#rec1 rec1 ..."
    @return structured numpy array
    """
    fIn = open(FileName,'r')
    lines = fIn.readlines()
    fIn.close()
    w = re.split('\s+',lines[0])
    w[0] = w[0][1:]
    dty = []
    for x in w:
        if x != '':
            dty.append((x,np.float))
    d = np.loadtxt(FileName,dtype=dty)
    return d


def GetStrCodeGitCmd(filepath,options,args):
    """
    Get script history and running informations:
    @param filepath is file path : os.path.realpath(__file__)
    @param options is the dictionary of options : vars(options)
    @param args is the list of arguments : args
    @return string with informations: script, git hash, git branch, commandline
    """
    dirCurrent = os.getcwd()+"/"
    dirScript  = os.path.dirname(filepath)+"/"
    #print(dirScript)
    nameScript = os.path.basename(filepath)
    #print(nameScript)
    os.chdir(dirScript)
    tmp = dirCurrent+'/tmp_GetStrCodeGitCmd.txt'
    run('cd '+dirScript+'; git rev-parse HEAD > '+tmp+' ; git branch >> '+tmp+' ; cd '+dirCurrent)
    fIn = open(tmp,'r')
    lines = fIn.readlines()
    fIn.close()
    gitHash = lines[0][:-1]
    for x in lines:
        if x[0]=='*':
            gitBranch = re.split("\s+",x)[1]
    os.chdir(dirCurrent)
    r = "#%s (%s): python3 %s"%(gitHash,gitBranch,filepath)
    for i,k in enumerate(options):
        r = r + " --"+str(k)+"="
        if type(options[k])==str:
            r = r + options[k]
        else:
            r = r + str(options[k])
    for x in args:
        r = r + " "+x
    return r
