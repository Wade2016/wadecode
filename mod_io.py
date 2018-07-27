def bootstraping(slist,n=None, repeat=10):
        import numpy as np
        if n==None:
                n=len(slist)
        out=[]
        for i in range(repeat):
                resample=np.floor(np.random.rand(n)*len(slist)).astype(int)
                temp=slist[resample]
                out.append(temp)
        return out

def savetxt(data, fname):
        with open(fname,'w') as f:
                for l in data:
                        f.write(l)
        #print('Saved file:',fname)
        return 0

