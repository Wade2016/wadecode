import sys
from subprocess import Popen, PIPE
import multiprocessing as mp
sys.path.insert(0,'/homea/ias-5/wenping/software/wadecode/')
import mod_io as mio

def gmx_api(cmd):
        process=Popen(cmd,stdout=PIPE)
        (output, err)=process.communicate()
        exit_code=process.wait()
        if err!=None:
                print("Error:", err)
        output=output.decode('utf-8')
        return output

def apbs_api(inp):
        process=Popen(['apbs',inp], stdout=PIPE)
        (output, err)=process.communicate()
        exit_code=process.wait()
        if err!=None:
                print("Error:", err)
        output=output.decode('utf-8')
        E_polar=None
        for l in output.split('\n'):
                if l[:18]=='Solvation energy =':
                        E_polar=float(l.split(' ')[-2])*4.184
        return output, E_polar

def delphi_api(inp):
        process=Popen(['delphi',inp], stdout=PIPE)
        (output, err)=process.communicate()
        exit_code=process.wait()
        if err!=None:
                print("Error:", err)
        output=output.decode('utf-8')
        E_polar=None
        for l in output.split('\n'):
                if l[:56]==' Energy> Corrected reaction field energy               :':
                        E_polar=float(l[56:].split(' ')[-2])*2.479
        return output, E_polar

def creat_index_by_element(pdb,elements):
        import os
        cmd0='gmx editconf -f '+pdb+' -legend -o sys.pdb '
        os.system(cmd0)
        atom=open('sys.pdb').readlines()
        atoms=[a for a in atom if a.split(' ')[0]=='ATOM' and len(a[:-1])==78 and not (a[77] in elements)]
        #print(atoms)
        heavy=[a[5:11] for a in atoms]
        heavy=['[heavy_atoms]\n']+heavy+['\n']
        #print(heavy)
        mio.savetxt(heavy,'heavy.ndx')
        return 'heavy.ndx'
