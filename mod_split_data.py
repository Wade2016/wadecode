def Split_data(testlist,numproc):
        # split data set into the number of processors
        numfile = len(testlist)
        setsize = int(numfile/numproc)
        resset  = numfile%numproc
        chunks  = [testlist[i:i+setsize] for i in range(0, numfile, setsize)]
        print("numproc:",numproc,"numfile:",numfile,"setsize:",setsize,"resset",resset)
        print("chunks:", len(chunks))
        return chunks

