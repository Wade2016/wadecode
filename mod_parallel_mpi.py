class Parallel(object):
        from mpi4py import MPI
        comm=MPI.COMM_WORLD
        rank=comm.rank
        size=comm.size

        def run(self, target, inputs):
                print("Inputs",inputs)
                output=target(self.rank, inputs)
                outputs=self.comm.gather(output,root=0)
                if self.rank==0:
                        print("MPI jobs are done!")
                        return outputs
