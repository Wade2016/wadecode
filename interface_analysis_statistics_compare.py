import pdb
import csv
import numpy as np
from itertools import groupby
from natsort import natsorted, ns
import pandas as pd
import sys, os, time
import multiprocessing as mp
import matplotlib.pylab as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


def Combine_Ensemble(filelist, keylist):
	frames=[pd.DataFrame.from_csv(str(fname)) for fname in filelist]
	Ensemble=pd.concat(frames, keys=keylist)
	EFrame=[]
	[[EFrame.append(iD+"_"+str(F)) for F in Ensemble.ix[iD]["Frame"].values] for iD in keylist]
	print(len(EFrame))
	Ensemble["EFrame"]=pd.DataFrame(EFrame, index=Ensemble.index)
	print(Ensemble.shape)
	del frames
	return Ensemble

def Creat_Matrix(chunk):
        Tbegin=time.time()
        # Create a dataframe to store contactID matrix
        NCI_contactID=pd.DataFrame(columns=ContactID)
        for frame in chunk:
                one_frame=NCI_data[NCI_data["EFrame"]==frame]
                #print("Frame:", frame)
                [NCI_contactID.set_value(frame,one_frame.loc[idx]["Cname"], one_frame.loc[idx]["Dist"]) for idx in one_frame.index]
        print(NCI_contactID.shape)
        outid=os.getpid()
        with open(".NCI_contactID_"+str(outid),  'a') as fout:
                NCI_contactID.to_csv(fout, header=False)
        print("LoopTime:", time.time()-Tbegin)
        return outid

def Creat_ContMatrix(NCI_data, factor):

        framelist=natsorted(set(NCI_data["EFrame"].values))
        chunks=Split_data(framelist, factor)
        #threading on
        print("Running:")
        pool = mp.Pool(processes = mp.cpu_count()-1)
        outid=pool.map(Creat_Matrix, chunks)
        pool.close()
        pool.join()
        print("Off")
        outid=list(set(outid))
        ContactMatrix=Merge_data(outid, ".NCI_contactID_", "NCI_contactID.csv")
        print(ContactMatrix.shape)

        return ContactMatrix

def Split_data(testlist, factor):
        # split data set into the number of processors
        numproc = (mp.cpu_count()-1)*factor
        numfile = len(testlist)
        setsize = int(numfile/numproc)
        resset  = numfile%numproc
        chunks  = [testlist[i:i+setsize] for i in range(0, numfile, setsize)]
        print("numproc:",numproc,"numfile:",numfile,"setsize:",setsize,"resset",resset)
        print("chunks:", len(chunks))
        return chunks

def Merge_data(outid, idfile, outfile):
        # Merge files to sigle file
        print("Merge data:")
        for i in natsorted(outid):
                print(i)
                os.system("cat "+ idfile + str(i) + " >> " + outfile)
        os.system("rm " + idfile+ "*")
        sdata=pd.DataFrame.from_csv(outfile)
        sdata=sdata.reindex(index=natsorted(sdata.index))
        sdata.to_csv(outfile)
        print(sdata.shape)
        return sdata

def map_pattern(colname, bty, mty, dataframe, level):
  #Attention, we must re-allocate memory for pattern DataFrame 
  zeros_data=np.zeros([level+5, len(dataframe.columns)])
  pattern=pd.DataFrame(zeros_data, columns=dataframe.columns)
  for ic in range(len(pattern.columns)):
    tp = ord(list(bty[ic])[1])
    for ir in range(level):
      if ir <= level*dataframe.loc[colname][ic]:
         pattern.set_value(pattern.index[ir], pattern.columns[ic] , tp)
    if "ALY" in mty[ic]:
      pattern.set_value(pattern.index[level:level+5], pattern.columns[ic], 5)
    elif "SEP" in mty[ic]:
      pattern.set_value(pattern.index[level:level+5], pattern.columns[ic], 15)
    else:
      pattern.set_value(pattern.index[level:level+5], pattern.columns[ic], 20)

  pattern.to_csv("pattern_"+colname+".csv")
  ptnew=pattern.where(pattern>0)
  ptnew.fillna(value=-100)
  fig, ax = plt.subplots()
  
  #pdb.set_trace()
  #cmap = mpl.cm.get_cmap('PiYG', 10)
  cax = ax.imshow(ptnew)

  #cb  = mpl.colorbar.ColorbarBase(ax2, ticks=[80,78,72,67], boundaries=[60,70,80, 90], format='%1i') 
  #ax1.set_title(colname)
#  cbar=fig.colorbar(cax, ticks=[80,78,72,67], orientation='horizontal')
#  cbar.ax.set_yticklabels(['P','N','H','C'])
  plt.savefig("pattern_"+colname +".jpg")
  #plt.show()
  plt.close()
  return ptnew
  
def creat_pattern(dataframe, level):
  bty=[]
  mty=[]
  for ti in dataframe.columns.values:
    bty.append(str.split(ti, "_")[2])
    mty.append(str.split(ti, "_")[0])
  for colname in ['WT', 'A', 'B', 'C', 'D']:
     print(colname)
     pattern=map_pattern(colname, bty, mty, dataframe, level)
  return 0

def plot_pattern():
   temp_wt=pd.DataFrame.from_csv("pattern_WT.csv")
   temp_a=pd.DataFrame.from_csv("pattern_A.csv")
   temp_b=pd.DataFrame.from_csv("pattern_B.csv")
   temp_c=pd.DataFrame.from_csv("pattern_C.csv")
   temp_d=pd.DataFrame.from_csv("pattern_D.csv")
   temp_wt=temp_wt.where(temp_wt>0)
   temp_a=temp_a.where(temp_a>0)
   temp_b=temp_b.where(temp_b>0)
   temp_c=temp_c.where(temp_c>0)
   temp_d=temp_d.where(temp_d>0)

   # Setting up a colormap that's a simple transtion
   mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red', 'black', 'green'])
   # Using contourf to provide my colorbar info, then clearing the figure
   Z = [[0,0],[0,0]]
   levels = [60, 67, 72, 78, 80, 85]
   CS3 = plt.contourf(Z, levels, cmap=mymap)
   plt.clf()

   fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, sharey=True)

   #cmap= plt.cm.get_cmap('PiYG', 4)

   im1=ax1.imshow(temp_wt, cmap=mymap)
   ax1.set_title('Protein-DNA Interaction Patterns')
   ax2.imshow(temp_a, cmap=mymap)
   ax3.imshow(temp_b, cmap=mymap)
   ax3.set_ylabel('Occupancy Rate (%)')
   ax4.imshow(temp_c, cmap=mymap)
   ax5.imshow(temp_d, cmap=mymap)
   ax1.annotate('WT', xy=(600, 80), xytext=(600, 80) )
   ax2.annotate('A', xy=(600, 80), xytext=(600, 80) )
   ax3.annotate('B', xy=(600, 80), xytext=(600, 80) )
   ax4.annotate('C', xy=(600, 80), xytext=(600, 80) )
   ax5.annotate('D', xy=(600, 80), xytext=(600, 80) )
   #divider = make_axes_locatable(ax6)
   #fig.delaxes(ax6)
   fig.subplots_adjust(hspace=0)
   plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
   plt.gca().invert_yaxis()
   plt.colorbar(CS3) 

#   cax=fig.append_axes('right', size='10%', pad=0.2)
#   fig.colorbar(im1, cax=cax, orientation='vertical')
#   data=["WT", "A", "B", "C", "D"]
#   for i in range(5):
#     temp=pd.DataFrame.from_csv("pattern_"+data[i]+".csv")
#     ("ax"+i).imshow(temp)
   plt.savefig("pattern_wt_iso.jpg")
   plt.show()
   #pdb.set_trace()
   return 0
  
if __name__=='__main__':
	Tbegin=time.time()
	#numxx=int(sys.argv[1])
	#factor=int(sys.argv[2])
	
	fopen1=sys.argv[1]
	fopen2=sys.argv[2]
	fopen3=sys.argv[3]
	fopen4=sys.argv[4]
	fopen5=sys.argv[5]
	fopen6=sys.argv[6]
	fopen7=sys.argv[7]
	fopen8=sys.argv[8]
	fopen9=sys.argv[9]
	fopen10=sys.argv[10]
	fopen11=sys.argv[11]
	fopen12=sys.argv[12]
	fopen13=sys.argv[13]
	fopen14=sys.argv[14]
	fopen15=sys.argv[15]
	factor =int(sys.argv[16])
	filelistA=[fopen1, fopen2, fopen3]
	filelistB=[fopen4, fopen5, fopen6]
	filelistC=[fopen7, fopen8, fopen9]
	filelistD=[fopen10, fopen11, fopen12]
	filelistWT=[fopen13, fopen14, fopen15]
	keyall=["A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3","W1","W2","W3"]
	print(filelistA)
	fileall=filelistA+filelistB+filelistC+filelistD+filelistWT
	
	"""
	NCI_data=Combine_Ensemble(fileall, keyall)
	print(NCI_data.shape)
	"""
	"""
	###################################################
	#NCI_data=pd.DataFrame.from_csv("EnsembleA.csv")
		
	# Creat contacts identifier: P_resid_Presname_Patom_PTM_D_resid_Dresname_Datom
	Contacts=NCI_data["P_resid"].map(str)+"_"+NCI_data["P_resname"]+"_"+NCI_data["P_atom"]+"_"+NCI_data["PTM"]+"_"+NCI_data["D_resid"].map(str)+"_"+NCI_data["D_resname"]+"_"+NCI_data["D_atom"]
	NCI_data["Cname"]=pd.Series(Contacts, index=Contacts.index)
	ContactID=list(natsorted(set(Contacts.values)))
	print("ContactID number is: ", len(ContactID))
	print("ContactID list built up!")
	del Contacts
	# keep low memory
	NCI_data=NCI_data[["EFrame","Dist","Type","Cname"]]
	# Create a new dataframe to store contactID matrix
	NCI_contactID=pd.DataFrame(columns=ContactID)
	NCI_contactID.to_csv("NCI_contactID.csv")
	NCI_data.to_csv("Ensemble_all.csv")

	
	
	keyall=["A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3","W1","W2","W3"]
	##################################################
	NCI_contactID=pd.DataFrame.from_csv("NCI_contactID.csv", index_col=[0,1])
	NCI_data=pd.DataFrame.from_csv("Nucplot_classified.csv.en", index_col=0)
	ContactID=list(NCI_contactID.columns.values)

	##################################################
	# creat contact matrix
	print("Creating contactMatrix Start")
	loopt=time.time()
	
	print("Constructing matrix for ", keyall[numxx])
	ContactMatrix=Creat_ContMatrix(NCI_data, factor)
	print("Time to construct matrix for ", ikey, "is: ", time.time()-loopt)
	print("Creating contactMatrix End")
	"""
  	##################################################	

	print("Running Time:", time.time()-Tbegin)
	print("Done")
