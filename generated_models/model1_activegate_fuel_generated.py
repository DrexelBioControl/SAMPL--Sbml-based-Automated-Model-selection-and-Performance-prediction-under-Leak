import sbmltoodepy.modelclasses
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Generated_model1_activegate_fuel(sbmltoodepy.modelclasses.Model):

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['k_txn'] = sbmltoodepy.modelclasses.Parameter(0.00962763757, 'k_txn', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['ksd'] = sbmltoodepy.modelclasses.Parameter(9.97801643e-06, 'ksd', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['kfsd'] = sbmltoodepy.modelclasses.Parameter(2.29972555e-06, 'kfsd', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['krev'] = sbmltoodepy.modelclasses.Parameter(1.42705205e-06, 'krev', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['kf_rep'] = sbmltoodepy.modelclasses.Parameter(4.77116458e-06, 'kf_rep', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['kRz'] = sbmltoodepy.modelclasses.Parameter(0.0406799972, 'kRz', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['klkg'] = sbmltoodepy.modelclasses.Parameter(6.98976431e-05, 'klkg', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.p['basal_frac'] = sbmltoodepy.modelclasses.Parameter(0.0220510273, 'basal_frac', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))

		self.c = {} #Dictionary of compartments
		self.c['default'] = sbmltoodepy.modelclasses.Compartment(1.0, 3, True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))

		self.s = {} #Dictionary of chemical species
		self.s['uRSDg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['RSDg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['IN'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['OUT'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['DRL'] = sbmltoodepy.modelclasses.Species(500.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['ROL'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['I_RSDg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['F'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['F_RSDg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['F_temp'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['IN_temp'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))
		self.s['RSD_temp'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))

		self.r = {} #Dictionary of reactions
		self.r['d_uRSDg_dt'] = d_uRSDg_dt(self)
		self.r['d_RSDg_dt'] = d_RSDg_dt(self)
		self.r['d_IN_dt'] = d_IN_dt(self)
		self.r['d_OUT_dt'] = d_OUT_dt(self)
		self.r['d_DRL_dt'] = d_DRL_dt(self)
		self.r['d_ROL_dt'] = d_ROL_dt(self)
		self.r['d_I_RSDg_dt'] = d_I_RSDg_dt(self)
		self.r['d_F_dt'] = d_F_dt(self)
		self.r['d_F_RSDg_dt'] = d_F_RSDg_dt(self)

		self.f = {} #Dictionary of function definitions
		self.time = 0

		self.AssignmentRules()



	def AssignmentRules(self):

		return

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['uRSDg'].amount, self.s['RSDg'].amount, self.s['IN'].amount, self.s['OUT'].amount, self.s['DRL'].amount, self.s['ROL'].amount, self.s['I_RSDg'].amount, self.s['F'].amount, self.s['F_RSDg'].amount, self.s['F_temp'].amount, self.s['IN_temp'].amount, self.s['RSD_temp'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[1,0,0,0,0,0,0,0,0.],[0,1,0,0,0,0,0,0,0.],[0,0,1,0,0,0,0,0,0.],[0,0,0,1,0,0,0,0,0.],[0,0,0,0,1,0,0,0,0.],[0,0,0,0,0,1,0,0,0.],[0,0,0,0,0,0,1,0,0.],[0,0,0,0,0,0,0,1,0.],[0,0,0,0,0,0,0,0,1.],[0,0,0,0,0,0,0,0,0.],[0,0,0,0,0,0,0,0,0.],[0,0,0,0,0,0,0,0,0.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['d_uRSDg_dt'](), self.r['d_RSDg_dt'](), self.r['d_IN_dt'](), self.r['d_OUT_dt'](), self.r['d_DRL_dt'](), self.r['d_ROL_dt'](), self.r['d_I_RSDg_dt'](), self.r['d_F_dt'](), self.r['d_F_RSDg_dt']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['uRSDg'].amount, self.s['RSDg'].amount, self.s['IN'].amount, self.s['OUT'].amount, self.s['DRL'].amount, self.s['ROL'].amount, self.s['I_RSDg'].amount, self.s['F'].amount, self.s['F_RSDg'].amount, self.s['F_temp'].amount, self.s['IN_temp'].amount, self.s['RSD_temp'].amount], dtype = np.float64)
		self.s['uRSDg'].amount, self.s['RSDg'].amount, self.s['IN'].amount, self.s['OUT'].amount, self.s['DRL'].amount, self.s['ROL'].amount, self.s['I_RSDg'].amount, self.s['F'].amount, self.s['F_RSDg'].amount, self.s['F_temp'].amount, self.s['IN_temp'].amount, self.s['RSD_temp'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class d_uRSDg_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return self.parent.p['k_txn'].value * self.parent.s['RSD_temp'].concentration - self.parent.p['kRz'].value * self.parent.s['uRSDg'].concentration

class d_RSDg_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return (self.parent.p['kRz'].value * self.parent.s['uRSDg'].concentration - self.parent.p['ksd'].value * self.parent.s['RSDg'].concentration * self.parent.s['IN'].concentration) + self.parent.p['krev'].value * self.parent.s['I_RSDg'].concentration * self.parent.s['OUT'].concentration - self.parent.p['klkg'].value * self.parent.s['RSDg'].concentration * self.parent.s['F'].concentration

class d_IN_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return (self.parent.p['k_txn'].value * self.parent.s['IN_temp'].concentration - self.parent.p['ksd'].value * self.parent.s['RSDg'].concentration * self.parent.s['IN'].concentration) + self.parent.p['krev'].value * self.parent.s['I_RSDg'].concentration * self.parent.s['OUT'].concentration + self.parent.p['kfsd'].value * self.parent.s['I_RSDg'].concentration * self.parent.s['F'].concentration

class d_OUT_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return (self.parent.p['ksd'].value * self.parent.s['RSDg'].concentration * self.parent.s['IN'].concentration - self.parent.p['kf_rep'].value * self.parent.s['OUT'].concentration * self.parent.s['DRL'].concentration - self.parent.p['krev'].value * self.parent.s['I_RSDg'].concentration * self.parent.s['OUT'].concentration) + self.parent.p['klkg'].value * self.parent.s['RSDg'].concentration * self.parent.s['F'].concentration + self.parent.p['basal_frac'].value * self.parent.p['k_txn'].value * self.parent.s['RSD_temp'].concentration

class d_DRL_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return -self.parent.p['kf_rep'].value * self.parent.s['OUT'].concentration * self.parent.s['DRL'].concentration

class d_ROL_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return self.parent.p['kf_rep'].value * self.parent.s['OUT'].concentration * self.parent.s['DRL'].concentration

class d_I_RSDg_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return self.parent.p['ksd'].value * self.parent.s['RSDg'].concentration * self.parent.s['IN'].concentration - self.parent.p['krev'].value * self.parent.s['I_RSDg'].concentration * self.parent.s['OUT'].concentration - self.parent.p['kfsd'].value * self.parent.s['I_RSDg'].concentration * self.parent.s['F'].concentration

class d_F_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return (self.parent.p['k_txn'].value * self.parent.s['F_temp'].concentration - self.parent.p['kfsd'].value * self.parent.s['I_RSDg'].concentration * self.parent.s['F'].concentration) + self.parent.p['ksd'].value * self.parent.s['F_RSDg'].concentration * self.parent.s['IN'].concentration - self.parent.p['klkg'].value * self.parent.s['RSDg'].concentration * self.parent.s['F'].concentration

class d_F_RSDg_dt:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")

	def __call__(self):
		return (self.parent.p['kfsd'].value * self.parent.s['I_RSDg'].concentration * self.parent.s['F'].concentration - self.parent.p['ksd'].value * self.parent.s['F_RSDg'].concentration * self.parent.s['IN'].concentration) + self.parent.p['klkg'].value * self.parent.s['RSDg'].concentration * self.parent.s['F'].concentration

