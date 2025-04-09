'''
main.py
Author: D. Younis
    University of Rochester, Department of Physics & Astronomy
Revised: 2/8/2023 at 0051
LEAP version: 0.8.0dev
'''
import os
import sys
import time
import h5py
import subprocess
import numpy as np
#import Primitives as P
from random import randint
from operator import attrgetter
from scipy.io import FortranFile

import leap_ec.ops as ops
from leap_ec.problem import ScalarProblem
from leap_ec.algorithm import generational_ea
from leap_ec import Individual, Representation
from leap_ec.real_rep.ops import mutate_gaussian

import toolz
from leap_ec.distrib import synchronous
from leap_ec.decoder import IdentityDecoder
from leap_ec.distrib.individual import DistributedIndividual

import multiprocessing.popen_spawn_posix
import distributed
from distributed import Client, LocalCluster

# usage: python main.py [serial/parallel] <input_file>.deck <output_file>.h5
args = sys.argv
if (args[2][-5:] != '.deck'): args[2] += '.deck'
if (args[3][-3:] != '.h5'): args[3] += '.h5'

'''Skeleton struct for Interior Minimization Parameters (Nelder-Mead over 11θ/7ϕ).'''
class imp:
    nmins = int(1) # no. interior minimizations to perform
    itmax = int(1) # max iterations per simplex
    ftol = 0. # cutoff fractional tolerance
    da = 0. # simplex leg step-size

'''Skeleton struct for Exterior Maximization Parameters (Genetic Algorithm over x,y,z).'''
class emp:
    # evolution parameters
    pop_size = int(1)
    max_generations = int(1)
    # mutate_gaussian parameters
    std = 0.
    expected_num_mutations = int(1)
    hard_bounds = [(0.,0.) for _ in range(3)]
    # uniform_crossover parameters
    p_swap = 0.
    # initial guess ranges [a,b)
    gx = np.array([0.,0.])
    gy = np.array([0.,0.])
    gz = np.array([0.,0.])

'''Read input deck parameters.'''
with open(args[2]) as fin:
    for ln in fin:
        sln = ln.split()
        if (len(sln) == 0): continue
        if (sln[0] == 'nmins:'): imp.nmins = int(eval(str(sln[1])))
        if (sln[0] == 'itmax:'): imp.itmax = int(eval(str(sln[1])))
        if (sln[0] == 'ftol:'): imp.ftol = float(eval(str(sln[1])))
        if (sln[0] == 'da:'): imp.da = np.pi*float(eval(str(sln[1])))
        if (sln[0] == 'e-ps:'): emp.pop_size = int(eval(str(sln[1])))
        if (sln[0] == 'e-mg:'): emp.max_generations = int(eval(str(sln[1])))
        if (sln[0] == 'e-std:'): emp.std = eval(str(sln[1]))
        if (sln[0] == 'e-expected_num_mutations:'):
            emp.expected_num_mutations = int(eval(str(sln[1])))
        if (sln[0] == 'e-hard_bounds:'):
            emp.hard_bounds = eval(''.join(sln[1:]))
        if (sln[0] == 'e-p_swap:'):
            emp.p_swap = float(sln[1])
        if (sln[0] == 'e-gx:'): emp.gx = eval(''.join(sln[1:]))
        if (sln[0] == 'e-gy:'): emp.gy = eval(''.join(sln[1:]))
        if (sln[0] == 'e-gz:'): emp.gz = eval(''.join(sln[1:]))

'''Output containers.'''
class out:
    # record of best fitness (sigma_12) value each generation
    fitness = np.zeros(emp.max_generations+1, dtype=float)
    # (x,y,z) vector of fittest individual
    X = np.zeros(3, dtype=float)
    # sigma_12 value of fittest individual
    cur_max = -np.inf

def write_data():
    '''Write data to HDF5 file.'''
    with h5py.File(args[3],'w') as hf:
        hf.create_dataset('fitness', data=out.fitness)
        hf.create_dataset('X', data=out.X)
        hf.create_dataset('cur_max', data=out.cur_max)

def init_xyz(gx,gy,gz):
    '''Initialize the (x,y,z) parameters randomly in [low,high).'''
    def create():
        return np.array([
            np.random.uniform(gx[0],gx[1]),
            np.random.uniform(gy[0],gy[1]),
            np.random.uniform(gz[0],gz[1])])
    return create

'''
Evolutionary algorithm problem for maximization of ObjectiveI1 (sigma_12) over (x,y,z).
'''
class MaxI1(ScalarProblem):
    def __init__(self):
        super().__init__(maximize=True)
    def evaluate(self,phenome):
        '''
        Outputs phenome[:3]=(x,y,z) to a Fortran-readable file with a 10-digit
        unique identifier (uid) and runs the Nelder-Mead binary (tetra_mini) on it.
        Upon completion, read the output from feval-{uid}.dat and delete temporary files.
        '''
        uid = ''.join(['{}'.format(randint(0,9)) for _ in range(10)])
        Xeval_name = 'Xeval-{}.dat'.format(uid)
        feval_name = 'feval-{}.dat'.format(uid)

        with FortranFile(Xeval_name,'w') as hf:
            hf.write_record(np.array(phenome))

        cmd_list = ['./tetra_mini', args[2], Xeval_name]
        process = subprocess.Popen(cmd_list, stdout=subprocess.PIPE)
        process.wait()

        if (process.returncode == 0):
            with FortranFile(feval_name,'r') as hf:
                ans = hf.read_record(float)[0]
        else:
            ans = -np.inf
            print('tetra_mini returned {} at X = ({})'
                .format(process.returncode, ', '.join(str(np.round(elt,3)) for elt in phenome)),
                flush=True)

        os.remove(Xeval_name)
        while os.path.isfile(feval_name):
            try:
                os.remove(feval_name)
            except OSError:
                continue

        return ans

#############################
### Main computation loop ###
#############################

tic = time.time()

if (args[1]=='serial'):
    # create evolutionary algorithm task
    ea = generational_ea(max_generations=emp.max_generations, pop_size=emp.pop_size, problem=MaxI1(),
        representation=Representation(individual_cls=Individual, initialize=init_xyz(emp.gx,emp.gy,emp.gz)),
        pipeline=[
            ops.tournament_selection,
            ops.clone,
            mutate_gaussian(
                std=emp.std,
                expected_num_mutations=emp.expected_num_mutations,
                hard_bounds=emp.hard_bounds
                ),
            ops.uniform_crossover(p_swap=emp.p_swap),
            ops.evaluate,
            ops.pool(size=emp.pop_size)
            ])

    # execute task, g = (generation #, Individual)
    n = 0
    for g in ea:
        # store current fitness value
        out.fitness[n] = g[1].fitness
        # output data if an improvement was made
        if (g[1].fitness > out.cur_max):
            out.X = g[1].genome
            out.cur_max = g[1].fitness
            write_data()
        # print info to log file
        if (np.mod(n, 0.02*emp.max_generations) == 0):
            print('obj {} gen {}, {}%'
                .format(str(g[1].fitness), n, 100.*n/emp.max_generations), flush=True)
            print('Xbest = ({})\n'
                .format(', '.join(str(np.round(elt,3)) for elt in g[1].genome)), flush=True)
        n += 1

    # finalize output data
    write_data()
    print('elapsed time: {} s'.format(np.round(time.time()-tic,2)), flush=True)

if (args[1]=='parallel' and __name__=='__main__'):
    with LocalCluster() as cluster, Client(cluster) as client:
        # initialize
        print(client, flush=True)
        parents = DistributedIndividual.create_population(emp.pop_size, problem=MaxI1(),
            initialize=init_xyz(emp.gx,emp.gy,emp.gz), decoder=IdentityDecoder())
        parents = synchronous.eval_population(parents, client=client)

        # execute task
        for n in range(emp.max_generations+1):
            offspring = toolz.pipe(
                parents,
                ops.tournament_selection,
                ops.clone,
                mutate_gaussian(
                    std=emp.std,
                    expected_num_mutations=emp.expected_num_mutations,
                    hard_bounds=emp.hard_bounds
                    ),
                ops.uniform_crossover(p_swap=emp.p_swap),
                synchronous.eval_pool(client=client, size=len(parents)))
            parents = offspring

            # store current fitness value
            this_best = max(offspring, key=attrgetter('fitness'))
            out.fitness[n] = this_best.fitness
            # output data if an improvement was made
            if (this_best.fitness > out.cur_max):
                out.X = this_best.genome
                out.cur_max = this_best.fitness
                write_data()
            # print info to log file
            if (np.mod(n, 0.02*emp.max_generations) == 0):
                print('obj {} gen {}, {}%'
                    .format(str(out.fitness[n]), n, 100.*n/emp.max_generations), flush=True)
                print('Xbest = ({})\n'
                    .format(', '.join(str(np.round(elt,3)) for elt in this_best.genome)), flush=True)

    # finalize output data
    write_data()
    print('elapsed time: {} s'.format(np.round(time.time()-tic,2)), flush=True)
