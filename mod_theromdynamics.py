import sys,os
import numpy as np
sys.path.insert(0,'/homea/ias-5/wenping/software/wadecode/')
import mod_io as mio
import mod_api as mapi
import mdtraj as md
import pandas as pd

class configS(object):
	print("Here we perform configurational entropy claulation for protein-dna association undergoes large conformational changes!")
	gmx="gmx_mpi"
	cutoff=0.2
	b=50000
	e=100000000
	dt=2
	maxwarn=5
	temperature=298.5
	index="index.ndx"
	index_comp="index_comp.ndx"
	pdb="pdb"
	traj="traj"
	trajfull="traj"
	tpr="tpr"
	top="top"
	cluster_cut=None
	mdp="mdp"
	nomp=1
	elements=["H","O","N"]
	heavy=False
	def cluster(self,run):
		if self.heavy==True:
			self.index=mapi.creat_index_by_element(self.pdb,self.elements)
		clustering="cluster,-f,"+self.traj+",-s,"+self.tpr+",-sz,-xvg,none,-cl,-clid,-method,gromos,-cutoff,"+str(self.cutoff)+",-n,"+self.index+",-b,"+str(self.b)+",-e,"+str(self.e)+",-dt,"+str(self.dt)

		print("Clustering using command:",clustering)
		print([self.gmx])
		print(clustering.split(","))
		clustering=[self.gmx]+clustering.split(",")
		print(clustering)
		if run==True:
			out=mapi.gmx_api(clustering)
		else:
			out="Culstering is not run this time!"
		mio.savetxt(out,"cluster_standard_out.log")
		clid=np.loadtxt("clust-id.xvg")
		clidpd=pd.DataFrame(clid)
		sz=np.loadtxt("clust-size.xvg")
		szpd=pd.DataFrame(sz)
		pdbs=open("clusters.pdb").readlines()
		cctime=[float(i[:-1].split(" ")[-1]) for i in pdbs if i[:5]=="TITLE"]
		szpd['center']=pd.Series(cctime,index=szpd.index)
		print(szpd)
		szpd.to_csv("clust-size.csv")
		clidpd.to_csv("clust-id.csv")
		return szpd,clidpd

	def statistics_estimate_S(self, run, pdb, xtc):
		print(print(pdb,xtc))
		traj=md.load(xtc,top=pdb)
		print('n_atoms:',traj.n_atoms)
		nframes=traj.n_frames
		slist=np.array(list(range(1,nframes)))
		resample=mio.bootstraping(slist, int(nframes/2.0), 10)
		i=0
		out=[]
		pwd=os.getcwd()
		for r in resample:
			os.system("mkdir resample_"+str(i))
			#os.chdir(pwd+"/resample_"+str(i))
			temp=traj.slice(r)
			temp.save_xtc('resample_traj_'+str(i)+".xtc")
			path=os.getcwd()
			print(path,pwd)
			self.vibS_tpr(run)
			os.chdir(pwd)
			i+=1
		os.system('rm *#')
		print(out)
		print("AV:",np.mean(np.array(out)),"STD:",np.std(np.array(out)))
		outpd=pd.DataFrame(out,index=range(len(out)),columns=[pdb])
		outpd.to_csv("freeE_-TdS.csv")
		return out

	def vibS_tpr(self, run):
		szpd,clipd=self.cluster(True)
		pdb_cmd=self.gmx+",editconf,-f,"+self.tpr+",-o,whole.pdb"
		mapi.gmx_api(pdb_cmd.split(","))
		print("Loading xtc data...")
		traj=md.load(self.trajfull,top="whole.pdb")
		print("xtc file has been loaded!")
		if self.cluster_cut==None:
			_i=input("Please chose a cutoff for individual configrations...\ncutoff=")
			print("Your choice:",_i)
			cut=int(_i)
		else:
			cut=self.cluster_cut
			print("Using coded cutoff",self.cluster_cut)
		szpd=szpd.where(szpd[1]>=cut).dropna()
		print("Performed cluster_cut by",cut)
		szpd["ratio"]=szpd[1]/szpd[1].sum()
		nfig=len(szpd)
		print(szpd)
		print("The number of individual configuration is",nfig)
		if run==False:
			return 0
		print("Loading xtc data...")
		traj=md.load(self.trajfull,top="whole.pdb")
		print("xtc file has been loaded!")
		outlog=[]
		tprs=[]
		for i in range(len(szpd)):
			c=szpd.center.values[i]
			fid=[t for t in range(traj.n_frames) if traj.time[t]==c]
			out=traj[fid]
			confname="iconfig_"+str(i)
			out.save_pdb(confname+".pdb")
			tpr_cmd=self.gmx+",grompp,-f,"+self.mdp+",-c,"+confname+".pdb,-n,"+self.index_comp+",-p,"+self.top+",-o,"+confname+",-maxwarn,"+str(self.maxwarn)
			tpr_cmd=tpr_cmd.split(",")
			temp=mapi.gmx_api(tpr_cmd)
			outlog.append("Config_"+str(i)+"\n"+temp+"\n")
			tprs.append(confname)
		mio.savetxt(outlog,"Config.log")
		return tprs

	def vibS_cal_mpi(self,run,tprname,trajname,oname,b,e):
		import mod_parallel_mpi as mpm
		_mpi=mpm.Parallel()
		if _mpi.rank==0:
			print("USING",_mpi.size,"MPI RANKS!")
		args=(self.index, tprname, trajname,b,e)
		output=_mpi.run(self.vibS_cal, args)
		if _mpi.rank==0:
			outpd=pd.concat(output)
			outpd.to_csv(oname)

	def vibS_cal(self, patch, args):
		print("Entropy\n",args)
		(index, tprname, trajname,tb,te)=args
		tpr=tprname+str(patch)+".tpr"
		traj=trajname+str(patch)
		eignve="evector"+str(patch)+".trr"
		eignva="evalue"+str(patch)+".xvg"

		# we have to re-define the inputs in this method
		#index=mapi.creat_index_by_element(pdb)
		cmd1='covar,-f,'+traj+',-s,'+tpr+',-n,'+index+',-v,'+eignve+',-o,'+eignva+',-b,'+str(tb)+',-e,'+str(te)
		out1=mapi.gmx_api([self.gmx]+cmd1.split(","))
		cmd2='anaeig,-entropy,-eig,'+eignva+',-v,'+eignve
		cmd2=cmd2.split(",")
		out=mapi.gmx_api([self.gmx]+cmd2)
		fe=None
		for l in out.split("\n"):
			print(l)
			if l[:54]=="The Entropy due to the Quasi Harmonic approximation is":
				entropy=[a for a in l[54:].split(" ") if a!=""][0]
				fe=-1*self.temperature*float(entropy)/1000.0
		mio.savetxt(out,"entropy_output.log"+str(patch))
		print("Quasi Harmonic approximation free energy change",fe,'kJ/mol')
		fepd=pd.DataFrame(fe,columns=['qhE_vib'],index=[patch])
		return fepd
	
	def extS_cal(self):
		pdb=self.pdb
		traj=self.traj
		data=md.load(traj,top=pdb)
		index=
		




