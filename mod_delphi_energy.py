import pandas as pd
import sys, glob

def get_energy(inp,outf):
	data=open(inp).readlines()
	out=pd.DataFrame()
	i=None
	factor=2.479
	for line in data:
		if line[:6]=="Frame:":
			i=int(line[7:-1])
			print("Frame:",i)
		if line[:56]==" Energy> Corrected reaction field energy               :":
			e=line[56:].split(' ')[-2]
			out.set_value(i,"E_psol",float(e)*factor)
		if line[:56]==" Energy> Coulombic energy                              :":
			e=line[56:].split(' ')[-2]
			out.set_value(i,"E_col",float(e)*factor)
		if line[:56]==" Energy> rho*phi/2 term in solution                    :":
			e=line[56:].split(' ')[-2]
			out.set_value(i,"rhophi2",float(e)*factor)
		if line[:56]==" Energy> All required energy terms but grid energy     :":
			e=line[56:].split(' ')[-2]
			out.set_value(i,"E_pball",float(e)*factor)
	print(out.shape)
	if outf!=None:
		out.to_csv(outf)
	return out

def main():
	inp=sys.argv[1]
	outf=sys.argv[2]
	out=[]
	inps=str(inp+"*")
	for log in glob.glob(inps):
		out.append(get_energy(log, None))
		print("File:",log)
	outall=pd.concat(out)
	print(outall.shape)
	print(outall.mean())
	outall.to_csv(outf)



if __name__=='__main__':
	main()
