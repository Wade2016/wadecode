import pandas as pd
import numpy as np

def delt_H():
	data={}
	for s in "noptm,a,b,c,d,ptdna".split(","):
		print("sys",s)
		data[s]=pd.read_csv("statistics_H_"+s+'.csv',index_col=0)
		data[s][["coul","psolv"]]*=1
		if s!="ptdna":
			data["comp_"+s]=pd.read_csv("statistics_H_comp_"+s+'.csv',index_col=0)
			data["comp_"+s][["coul","psolv"]]*=1
	print(data.keys())
	deltH={}
	deltH_mean=pd.DataFrame(columns="noele,mmcoul,coul,psolv,npsolv,deltH,deltHmm".split(","))
	deltH_std=pd.DataFrame(columns="noele,mmcoul,coul,psolv,npsolv,deltH,deltHmm".split(","))
	deltHpd={}
	for s in "noptm,a,b,c,d".split(","):
		deltH[s]=data["comp_"+s]-data[s]-data["ptdna"]
		deltH[s]["deltH"]=deltH[s]["noele,coul,psolv,npsolv".split(",")].sum(axis=1)
		deltH[s]["deltHmm"]=deltH[s]["noele,mmcoul,psolv,npsolv".split(",")].sum(axis=1)
		deltHpd[s]=deltH[s]["deltH"]
		#print(deltH[s].mean(axis=0))
		deltH_mean.loc[s]=pd.Series(deltH[s].mean(axis=0).values,index=deltH[s].columns)
		deltH_std.loc[s]=pd.Series(deltH[s].std(axis=0).values,index=deltH[s].columns)
	deltH_mean.to_csv("deltH_mean.csv")
	deltH_std.to_csv("deltH_std.csv")
	print(deltH_mean,'\n',deltH_std)
	return pd.DataFrame.from_dict(deltHpd)

def delt_S():
	data={}
	for s in "noptm,isoa,isob,isoc,isod,ptdna".split(","):
		print("sys",s)
		data[s]=pd.read_csv(s+"_mmpbsa.csv",index_col=0).iloc[:,0]
		if s !="ptdna":
			data["comp_"+s]=pd.read_csv("comp_"+s+"_mmpbsa.csv",index_col=0).iloc[:,0]
	datapd=pd.DataFrame.from_dict(data)
	print(datapd)
	datapd["deltS_noptm"]=datapd.comp_noptm-datapd.ptdna-datapd.noptm
	datapd["deltS_isoa"]=datapd.comp_isoa-datapd.ptdna-datapd.isoa
	datapd["deltS_isob"]=datapd.comp_isob-datapd.ptdna-datapd.isob
	datapd["deltS_isoc"]=datapd.comp_isoc-datapd.ptdna-datapd.isoc
	datapd["deltS_isod"]=datapd.comp_isod-datapd.ptdna-datapd.isod
	print(datapd.mean(axis=0))
	print(datapd.std(axis=0))
	#print(datapd)
	return datapd.iloc[:,-5:]

def bootstraping(slist,n=None, repeat=10):
        if n==None:
                n=len(slist)
        out=[]
        for i in range(repeat):
                resample=np.floor(np.random.rand(n)*len(slist)).astype(int)
                temp=slist[resample]
                out.append(temp)
        return out
	
def delt_G(deltH,deltS):
	lenH=len(deltH.index)
	lenS=len(deltS.index)
	print(deltH)
	print(deltS)
	if lenH>lenS:
		print(lenH, lenS)
		lists=bootstraping(deltH.index, int(lenH/10.0), lenS)
		Hnew={}
		for s in range(lenS):
			Hnew[s]=deltH.loc[lists[s]].mean(axis=0)
		#print(Hnew)
		Hnewpd=pd.DataFrame.from_dict(Hnew,orient='index')
		print(Hnewpd)
	deltG={}
	deltG['noptm']=Hnewpd['noptm']+deltS['deltS_noptm']
	deltG['isoa']=Hnewpd['a']+deltS['deltS_isoa']
	deltG['isob']=Hnewpd['b']+deltS['deltS_isob']
	deltG['isoc']=Hnewpd['c']+deltS['deltS_isoc']
	deltG['isod']=Hnewpd['d']+deltS['deltS_isod']
	print(deltG)	
	deltG=pd.DataFrame.from_dict(deltG)
	print(deltG.mean(axis=0))
	print(deltG.std(axis=0))
	

deltH=delt_H()
deltS=delt_S()
deltG=delt_G(deltH,deltS)
