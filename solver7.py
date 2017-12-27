import argparse
import numpy as np
from collections import defaultdict
from collections import Counter
from itertools import permutations, product
import random
from simanneal import Annealer

"""
======================================================================
  Complete the following function.
======================================================================
"""

def solve(num_wizards, num_constraints, wizards, constraints):
    """
    Write your algorithm here.
    potential_ordering: dict{name : order number}
    Input:
        num_wizards: Number of wizards
        num_constraints: Number of constraints
        wizards: An array of wizard names, in no particular order
        constraints: A 2D-array of constraints, 
                     where constraints[0] may take the form ['A', 'B', 'C']i

    Output:
        An array of wizard names in the ordering your algorithm returns
    """
    
    np.random.shuffle(wizards)
    np.random.shuffle(constraints)
    initialState = wizards[:]
    tsp = WizardProblem(initialState, num_wizards, num_constraints, wizards, constraints)
    auto_schedule = tsp.auto(minutes=20)
    tsp.set_schedule(auto_schedule) 
    itinerary, miles = tsp.anneal()
    return itinerary

def getBrokenConstraints(results, constraints):
    brokenConstraints = []
    for constraint in constraints:
        wiz_a = constraint[0]
        wiz_b = constraint[1]
        wiz_mid = constraint[2]
        if wiz_a in results and wiz_b in results and wiz_mid in results:
            wiz_a_loc = results.index(wiz_a)
            wiz_b_loc = results.index(wiz_b)
            wiz_mid_loc = results.index(wiz_mid)
            if (wiz_a_loc < wiz_mid_loc < wiz_b_loc) or (wiz_b_loc < wiz_mid_loc < wiz_a_loc):
                brokenConstraints = brokenConstraints + [constraint]
    return brokenConstraints

class WizardProblem(Annealer):

        def __init__(self, state, num_wizards, num_constraints, wizards, constraints):
            self.num_wizards = num_wizards
            self.num_constraints = num_constraints
            self.wizards = wizards
            self.constraints = constraints
            self.opt = None
            super(WizardProblem, self).__init__(state)

        def move(self):
            brokenConstraints = getBrokenConstraints(self.state, self.constraints)
            if len(brokenConstraints) == 1:
                rand_num = 0
                constraintToFix = brokenConstraints[rand_num]
                wiz_mid = self.state.index(constraintToFix[2])
                wiz_b = self.state.index(constraintToFix[0])
                self.state[wiz_mid], self.state[wiz_b] = self.state[wiz_b], self.state[wiz_mid]
            
            if len(brokenConstraints) != 0 and len(brokenConstraints) !=1:
                rand_num = random.randint(0,len(brokenConstraints)-1)
                constraintToFix = brokenConstraints[rand_num]
                rand_num = random.random()
                if rand_num >= .5:
                    wiz_mid = self.state.index(constraintToFix[2])
                    wiz_b = self.state.index(constraintToFix[0])
                    wiz_a = self.state.index(constraintToFix[1])
                    self.state[wiz_mid], self.state[wiz_b] = self.state[wiz_b], self.state[wiz_mid]
                else:
                    wiz_mid = self.state.index(constraintToFix[2])
                    wiz_b = self.state.index(constraintToFix[1])
                    wiz_a = self.state.index(constraintToFix[0])
                    self.state[wiz_mid], self.state[wiz_b] = self.state[wiz_b], self.state[wiz_mid]
        
        def energy(self):
            brokenConstraints = getBrokenConstraints(self.state, self.constraints)
            return len(brokenConstraints)

"""
======================================================================
   No need to change any code below this line
======================================================================
"""

def read_input(filename):
    with open(filename) as f:
        num_wizards = int(f.readline())
        num_constraints = int(f.readline())
        constraints = []
        wizards = set()
        for _ in range(num_constraints):
            c = f.readline().split()
            constraints.append(c)
            for w in c:
                wizards.add(w)
                
    wizards = list(wizards)
    return num_wizards, num_constraints, wizards, constraints

def write_output(filename, solution):
    with open(filename, "w") as f:
        for wizard in solution:
            f.write("{0} ".format(wizard))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Constraint Solver.")
    parser.add_argument("input_file", type=str, help = "___.in")
    parser.add_argument("output_file", type=str, help = "___.out")
    args = parser.parse_args()

    num_wizards, num_constraints, wizards, constraints = read_input(args.input_file)
    solution = solve(num_wizards, num_constraints, wizards, constraints)
    write_output(args.output_file, solution)
