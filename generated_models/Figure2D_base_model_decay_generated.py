import sbmltoodepy.modelclasses
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Generated_Figure2D_base_model_decay(sbmltoodepy.modelclasses.Model):

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['ktx__'] = sbmltoodepy.modelclasses.Parameter(0.013, 'ktx__', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['kmat__'] = sbmltoodepy.modelclasses.Parameter(0.00417, 'kmat__', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['kx__'] = sbmltoodepy.modelclasses.Parameter(1e-06, 'kx__', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['krev__'] = sbmltoodepy.modelclasses.Parameter(2.7e-07, 'krev__', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['krep__'] = sbmltoodepy.modelclasses.Parameter(1e-05, 'krep__', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['k_test_decay'] = sbmltoodepy.modelclasses.Parameter(0.0001, 'k_test_decay', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))

		self.c = {} #Dictionary of compartments
		self.c['default'] = sbmltoodepy.modelclasses.Compartment(1e-06, 3, True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("default"))

		self.s = {} #Dictionary of chemical species
		self.s['dna_Input_1'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("dna_Input_1"))
		self.s['rna_Input_1'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rna_Input_1"))
		self.s['dna_Gate_1_2'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("dna_Gate_1_2"))
		self.s['rna_Gate_1_2_pre'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rna_Gate_1_2_pre"))
		self.s['rna_Gate_1_2'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rna_Gate_1_2"))
		self.s['rna_Output_1_2'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rna_Output_1_2"))
		self.s['dna_Signal_2_ss'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("dna_Signal_2_ss"))
		self.s['rna_GateOut_1_1'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rna_GateOut_1_1"))
		self.s['dna_Reporter_2'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("dna_Reporter_2"))
		self.s['rnadna_ReporterOut_2'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rnadna_ReporterOut_2"))

		self.r = {} #Dictionary of reactions
		self.r['r0'] = r0(self)
		self.r['r1'] = r1(self)
		self.r['r2'] = r2(self)
		self.r['r3'] = r3(self)
		self.r['r3rev'] = r3rev(self)
		self.r['r4'] = r4(self)
		self.r['r_decay'] = r_decay(self)

		self.f = {} #Dictionary of function definitions
		self.time = 0

		self.AssignmentRules()



	def AssignmentRules(self):

		return

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['dna_Input_1'].amount, self.s['rna_Input_1'].amount, self.s['dna_Gate_1_2'].amount, self.s['rna_Gate_1_2_pre'].amount, self.s['rna_Gate_1_2'].amount, self.s['rna_Output_1_2'].amount, self.s['dna_Signal_2_ss'].amount, self.s['rna_GateOut_1_1'].amount, self.s['dna_Reporter_2'].amount, self.s['rnadna_ReporterOut_2'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[ 0,0,0,0,0,0,0.],[ 1,0,0,-1,1,0,0.],[ 0,0,0,0,0,0,0.],[ 0,1,-1,0,0,0,0.],[ 0,0,1,-1,1,0,0.],[ 0,0,0,1,-1,-1,-1.],[ 0,0,0,0,0,1,0.],[ 0,0,0,1,-1,0,0.],[ 0,0,0,0,0,-1,0.],[ 0,0,0,0,0,1,0.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['r0'](), self.r['r1'](), self.r['r2'](), self.r['r3'](), self.r['r3rev'](), self.r['r4'](), self.r['r_decay']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['dna_Input_1'].amount, self.s['rna_Input_1'].amount, self.s['dna_Gate_1_2'].amount, self.s['rna_Gate_1_2_pre'].amount, self.s['rna_Gate_1_2'].amount, self.s['rna_Output_1_2'].amount, self.s['dna_Signal_2_ss'].amount, self.s['rna_GateOut_1_1'].amount, self.s['dna_Reporter_2'].amount, self.s['rnadna_ReporterOut_2'].amount], dtype = np.float64)
		self.s['dna_Input_1'].amount, self.s['rna_Input_1'].amount, self.s['dna_Gate_1_2'].amount, self.s['rna_Gate_1_2_pre'].amount, self.s['rna_Gate_1_2'].amount, self.s['rna_Output_1_2'].amount, self.s['dna_Signal_2_ss'].amount, self.s['rna_GateOut_1_1'].amount, self.s['dna_Reporter_2'].amount, self.s['rnadna_ReporterOut_2'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class r0:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r0")

	def __call__(self):
		return self.parent.c['default'].size * self.parent.p['ktx__'].value * self.parent.s['dna_Input_1'].concentration

class r1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r1")

	def __call__(self):
		return self.parent.c['default'].size * self.parent.p['ktx__'].value * self.parent.s['dna_Gate_1_2'].concentration

class r2:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r2")

	def __call__(self):
		return self.parent.c['default'].size * self.parent.p['kmat__'].value * self.parent.s['rna_Gate_1_2_pre'].concentration

class r3:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r3")

	def __call__(self):
		return self.parent.c['default'].size * self.parent.p['kx__'].value * self.parent.s['rna_Gate_1_2'].concentration * self.parent.s['rna_Input_1'].concentration

class r3rev:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r3rev")

	def __call__(self):
		return self.parent.c['default'].size * self.parent.p['krev__'].value * self.parent.s['rna_Output_1_2'].concentration * self.parent.s['rna_GateOut_1_1'].concentration

class r4:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r4")

	def __call__(self):
		return self.parent.c['default'].size * self.parent.p['krep__'].value * self.parent.s['rna_Output_1_2'].concentration * self.parent.s['dna_Reporter_2'].concentration

class r_decay:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r_decay")

	def __call__(self):
		return self.parent.c['default'].size * self.parent.p['k_test_decay'].value * self.parent.s['rna_Output_1_2'].concentration

