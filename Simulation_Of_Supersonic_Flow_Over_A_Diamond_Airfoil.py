from tkinter import Tk
from tkinter.ttk import Label, Button, Entry, Frame
from numpy import sin, cos, tan, atan, sqrt, radians, ndarray, object_, float64
from scipy.optimize import fsolve
from numpy.linalg import solve
from sympy import sin as s, cos as c, tan as t, sqrt as sq, symbols, rad as r, cot, atan as at

class IsentropicFlow:

	def sat_prop(prop,stat_value,m,gamma):

		mGammaFun = 1 + ((gamma - 1)/2)*(m**2)
		exp = 1

		if prop == 'temperature':
			exp = 1

		elif prop == 'pressure':
			exp = gamma/(gamma - 1)

		elif prop == 'density':
			exp = 1/(gamma - 1)

		else:
			pass

		return stat_value*mGammaFun**exp

	def stat_prop(prop,sat_value,m,gamma):

		mGammaFun = 1 + ((gamma - 1)/2)*(m**2)
		exp = 1

		if prop == 'temperature':
			exp = 1

		elif prop == 'pressure':
			exp = gamma/(gamma - 1)

		elif prop == 'density':
			exp = 1/(gamma - 1)

		else:
			pass

		return sat_value/(mGammaFun**exp)

class NormalShockWave:

	def mach_downstream(m,gamma):

		num = 2/(gamma - 1) + m**2
		den = (2*gamma/(gamma - 1))*(m**2) - 1

		return sqrt(num/den)

	def prop_downstream(prop,upstream_value,m,gamma):

		mDown = NormalShockWave.mach_downstream(m,gamma)
		mUp = m

		if prop == 'temperature':
			result = upstream_value*(1 + ((gamma - 1)/2)*(mUp**2))/(1 + ((gamma - 1)/2)*(mDown**2))

		elif prop == 'pressure':
			result = upstream_value*(1 + gamma*mUp**2)/(1 + gamma*mDown**2)

		else:
			pass

		return result

class ObliqueShockWave:

	def theta_beta_m(beta,theta,m,gamma):

		return tan(radians(theta)) - ((m**2)*sin(radians(2*beta)) - 2/tan(radians(beta)))/(2 + (m**2)*(gamma + cos(radians(2*beta))))

	def mach_downstream(theta,m,gamma):

		beta = fsolve(ObliqueShockWave.theta_beta_m,45,args = (theta,m,gamma))[0]
		mUpNormal = m*sin(radians(beta))
		mDownNormal = NormalShockWave.mach_downstream(mUpNormal,gamma)
		mDown = mDownNormal/sin(radians(beta - theta))

		return mDown

	def prop_downstream(prop,upstream_value,theta,m,gamma):

		beta = fsolve(ObliqueShockWave.theta_beta_m,45,args = (theta,m,gamma))[0]
		mUpNormal = m*sin(radians(beta))
		mDownNormal = NormalShockWave.mach_downstream(mUpNormal,gamma)

		result = NormalShockWave.prop_downstream(prop,upstream_value,mUpNormal,gamma)

		return result

class ExpansionWave:

	def omega(m,gamma):

		expr1 = sqrt((gamma + 1)/(gamma - 1))
		expr2 = 1/expr1
		expr3 = sqrt(m**2 - 1)

		return expr1*atan(expr2*expr3) - atan(expr3)

	def symomega(m,gamma):

		expr1 = sq((gamma + 1)/(gamma - 1))
		expr2 = 1/expr1
		expr3 = sq(m**2 - 1)

		return expr1*at(expr2*expr3) - at(expr3)

	def prandtl_meyer_expansion_function(m2,m1,gamma,theta):

		return radians(theta) - ExpansionWave.omega(m2,gamma) + ExpansionWave.omega(m1,gamma)

class SupersonicFlowSimulation:

	def __init__(self,master):

		def simulate():

			p1 = float(self.p1Entry.get())
			t1 = float(self.t1Entry.get())
			m1 = float(self.m1Entry.get())
			theta = float(self.thetaEntry.get())
			alpha = float(self.alphaEntry.get())
			gamma = float(self.gammaEntry.get())

			if alpha >= 0 and theta > alpha:

				p2 = ObliqueShockWave.prop_downstream('pressure', p1, theta - alpha, m1, gamma)
				t2 = ObliqueShockWave.prop_downstream('temperature', t1, theta - alpha, m1, gamma)
				m2 = ObliqueShockWave.mach_downstream(theta - alpha, m1, gamma)

				p4 = ObliqueShockWave.prop_downstream('pressure', p1, theta + alpha, m1, gamma)
				t4 = ObliqueShockWave.prop_downstream('temperature', t1, theta + alpha, m1, gamma)
				m4 = ObliqueShockWave.mach_downstream(theta + alpha, m1, gamma)

			if alpha >= 0 and theta == alpha:

				p2, t2, m2 = p1, t1, m1

				p4 = ObliqueShockWave.prop_downstream('pressure', p1, theta + alpha, m1, gamma)
				t4 = ObliqueShockWave.prop_downstream('temperature', t1, theta + alpha, m1, gamma)
				m4 = ObliqueShockWave.mach_downstream(theta + alpha, m1, gamma)

			if alpha >= 0 and theta < alpha:

				m2 = fsolve(ExpansionWave.prandtl_meyer_expansion_function,m1,args=(m1,gamma,theta))[0]

				p01 = IsentropicFlow.sat_prop('pressure', p1, m1, gamma)
				p02 = p01
				p2 = IsentropicFlow.stat_prop('pressure', p02, m2, gamma)

				t01 = IsentropicFlow.sat_prop('temperature', t1, m1, gamma)
				t02 = t01
				t2 = IsentropicFlow.stat_prop('temperature', t02, m2, gamma)

				p4 = ObliqueShockWave.prop_downstream('pressure', p1, theta + alpha, m1, gamma)
				t4 = ObliqueShockWave.prop_downstream('temperature', t1, theta + alpha, m1, gamma)
				m4 = ObliqueShockWave.mach_downstream(theta + alpha, m1, gamma)

			if alpha < 0 and abs(theta) > abs(alpha):

				p2 = ObliqueShockWave.prop_downstream('pressure', p1, theta + alpha, m1, gamma)
				t2 = ObliqueShockWave.prop_downstream('temperature', t1, theta + alpha, m1, gamma)
				m2 = ObliqueShockWave.mach_downstream(theta + alpha, m1, gamma)

				p4 = ObliqueShockWave.prop_downstream('pressure', p1, theta - alpha, m1, gamma)
				t4 = ObliqueShockWave.prop_downstream('temperature', t1, theta - alpha, m1, gamma)
				m4 = ObliqueShockWave.mach_downstream(theta - alpha, m1, gamma)

			if alpha < 0 and abs(theta) == abs(alpha):

				p4, t4, m4 = p1, t1, m1

				p2 = ObliqueShockWave.prop_downstream('pressure', p1, theta + alpha, m1, gamma)
				t2 = ObliqueShockWave.prop_downstream('temperature', t1, theta + alpha, m1, gamma)
				m2 = ObliqueShockWave.mach_downstream(theta + alpha, m1, gamma)

			if alpha < 0 and abs(theta) < abs(alpha):

				m4 = fsolve(ExpansionWave.prandtl_meyer_expansion_function,m1,args=(m1,gamma,theta))[0]

				p01 = IsentropicFlow.sat_prop('pressure', p1, m1, gamma)
				p04 = p01
				p4 = IsentropicFlow.stat_prop('pressure', p04, m4, gamma)

				t04 = IsentropicFlow.sat_prop('temperature', t1, m1, gamma)
				t04 = t01
				t4 = IsentropicFlow.stat_prop('temperature', t04, m4, gamma)

				p2 = ObliqueShockWave.prop_downstream('pressure', p1, theta + alpha, m1, gamma)
				t2 = ObliqueShockWave.prop_downstream('temperature', t1, theta + alpha, m1, gamma)
				m2 = ObliqueShockWave.mach_downstream(theta + alpha, m1, gamma)

			m3 = fsolve(ExpansionWave.prandtl_meyer_expansion_function,m2,args=(m2,gamma,theta))[0]

			p02 = IsentropicFlow.sat_prop('pressure', p2, m2, gamma)
			p03 = p02
			p3 = IsentropicFlow.stat_prop('pressure', p03, m3, gamma)

			t02 = IsentropicFlow.sat_prop('temperature', t2, m2, gamma)
			t03 = t02
			t3 = IsentropicFlow.stat_prop('temperature', t03, m3, gamma)

			m5 = fsolve(ExpansionWave.prandtl_meyer_expansion_function,m2,args=(m2,gamma,theta))[0]

			p04 = IsentropicFlow.sat_prop('pressure', p4, m4, gamma)
			p05 = p04
			p5 = IsentropicFlow.stat_prop('pressure', p05, m5 , gamma)

			t04 = IsentropicFlow.sat_prop('temperature', t4, m4, gamma)
			t05 = t04
			t5 = IsentropicFlow.stat_prop('temperature', t05, m5, gamma)

			p6_case1 = ObliqueShockWave.prop_downstream('pressure', p3, 2*theta, m3, gamma)
			p7_case1 = p5

			p6_case2 = p3
			p7_case2 = ObliqueShockWave.prop_downstream('pressure', p5, 2*theta, m5, gamma)

			if p6_case1 == p7_case1:

				p6_final = p6_case1
				t6_final = ObliqueShockWave.prop_downstream('temperature', t3, 2*theta, m3, gamma)
				m6_final = ObliqueShockWave.mach_downstream(2*theta, m3, gamma)

				p7_final = p7_case1
				t7_final = t5
				m7_final = m5

				phi_final = theta
				direction_final = 'Upward'

				self.p2ValueLabel.configure(text = str(float(p2)))
				self.t2ValueLabel.configure(text = str(float(t2)))
				self.m2ValueLabel.configure(text = str(float(m2)))

				self.p3ValueLabel.configure(text = str(float(p3)))
				self.t3ValueLabel.configure(text = str(float(t3)))
				self.m3ValueLabel.configure(text = str(float(m3)))

				self.p4ValueLabel.configure(text = str(float(p4)))
				self.t4ValueLabel.configure(text = str(float(t4)))
				self.m4ValueLabel.configure(text = str(float(m4)))

				self.p5ValueLabel.configure(text = str(float(p5)))
				self.t5ValueLabel.configure(text = str(float(t5)))
				self.m5ValueLabel.configure(text = str(float(m5)))

				self.p6ValueLabel.configure(text = str(float(p6_final)))
				self.t6ValueLabel.configure(text = str(float(t6_final)))
				self.m6ValueLabel.configure(text = str(float(m6_final)))

				self.p7ValueLabel.configure(text = str(float(p7_final)))
				self.t7ValueLabel.configure(text = str(float(t7_final)))
				self.m7ValueLabel.configure(text = str(float(m7_final)))

				self.phiValueLabel.configure(text = str(float(phi_final)))
				self.dirValueLabel.configure(text = str(direction_final))

			if p6_case2 == p7_case2:

				p6_final = p6_case2
				t6_final = t3
				m6_final = m3

				p7_final = ObliqueShockWave.prop_downstream('pressure', p3, 2*theta, m3, gamma)
				t7_final = ObliqueShockWave.prop_downstream('temperature', t3, 2*theta, m3, gamma)
				m7_final = ObliqueShockWave.mach_downstream(2*theta, m5, gamma)

				phi_final = theta
				direction_final = 'Downward'

				self.p2ValueLabel.configure(text = str(float(p2)))
				self.t2ValueLabel.configure(text = str(float(t2)))
				self.m2ValueLabel.configure(text = str(float(m2)))

				self.p3ValueLabel.configure(text = str(float(p3)))
				self.t3ValueLabel.configure(text = str(float(t3)))
				self.m3ValueLabel.configure(text = str(float(m3)))

				self.p4ValueLabel.configure(text = str(float(p4)))
				self.t4ValueLabel.configure(text = str(float(t4)))
				self.m4ValueLabel.configure(text = str(float(m4)))

				self.p5ValueLabel.configure(text = str(float(p5)))
				self.t5ValueLabel.configure(text = str(float(t5)))
				self.m5ValueLabel.configure(text = str(float(m5)))

				self.p6ValueLabel.configure(text = str(float(p6_final)))
				self.t6ValueLabel.configure(text = str(float(t6_final)))
				self.m6ValueLabel.configure(text = str(float(m6_final)))

				self.p7ValueLabel.configure(text = str(float(p7_final)))
				self.t7ValueLabel.configure(text = str(float(t7_final)))
				self.m7ValueLabel.configure(text = str(float(m7_final)))

				self.phiValueLabel.configure(text = str(float(phi_final)))
				self.dirValueLabel.configure(text = str(direction_final))

			if p6_case1 > p7_case1 and p6_case2 < p7_case2:

				def guess():

					beta3Guess = float(beta3GuessEntry.get())
					beta5Guess = float(beta5GuessEntry.get())
					theta3Guess = float(theta3GuessEntry.get())
					theta5Guess = float(theta5GuessEntry.get())

					eqn = ndarray((1,10),dtype = object_)[0]
					syms = ['p6','t6','m6','p7','t7','m7','beta3','beta5','theta3','theta5']
					x0 = [p3,t3,m3,p5,t5,m5,beta3Guess,beta5Guess,theta3Guess,theta5Guess]

					p6,t6,m6,p7,t7,m7,beta3,beta5,theta3,theta5 = symbols(syms)
					symvars = [p6,t6,m6,p7,t7,m7,beta3,beta5,theta3,theta5]

					subs_dict = dict()

					eqn[0] = p6 - p7
					eqn[1] = theta3 + theta5 - 2*theta
					eqn[2] = t(r(theta3)) - ((m3**2)*s(r(2*beta3)) - 2*cot(r(2*beta3)))/(2 + (m3**2)*(gamma + c(r(2*beta3))))
					eqn[3] = t(r(theta5)) - ((m5**2)*s(r(2*beta5)) - 2*cot(r(2*beta5)))/(2 + (m5**2)*(gamma + c(r(2*beta5))))
					eqn[4] = p6 - p3*(1 + gamma*(m3*s(r(beta3)))**2)/(1 + gamma*(m6*s(r(beta3 - theta3)))**2)
					eqn[5] = t6 - t3*(1 + ((gamma - 1)/2)*(m3*s(r(beta3)))**2)/(1 + ((gamma - 1)/2)*(m6*s(r(beta3 - theta3)))**2)
					eqn[6] = p7 - p5*(1 + gamma*(m5*s(r(beta5)))**2)/(1 + gamma*(m7*s(r(beta5 - theta5)))**2)
					eqn[7] = t7 - t5*(1 + ((gamma - 1)/2)*(m5*s(r(beta5)))**2)/(1 + ((gamma - 1)/2)*(m7*s(r(beta3 - theta3)))**2)
					eqn[8] = m6*s(r(beta3 - theta3)) - (((2/(gamma - 1)) + (m3*s(r(beta3)))**2)/(2*gamma/(gamma - 1)*(m3*s(r(beta3)))**2 - 1))**(1/2)
					eqn[9] = m7*s(r(beta5 - theta5)) - (((2/(gamma - 1)) + (m5*s(r(beta5)))**2)/(2*gamma/(gamma - 1)*(m5*s(r(beta5)))**2 - 1))**(1/2)

					j = ndarray((10,10),dtype = object_)

					eqn0 = ndarray((1,10),dtype = float64)[0]

					for iterno in range(100):

						for i in range(len(syms)):subs_dict[syms[i]] = x0[i]

						for k in range(len(eqn)):eqn0[k] = float(eqn[k].evalf(subs = subs_dict))

						for i_ in range(len(eqn)):
							for j_ in range(len(symvars)):
								j[i][j_] = eqn[i].diff(symvars[j_]).evalf(subs = subs_dict)

						y = solve(float64(j),float64(eqn0.reshape((10,1)))).reshape(10)
						x0 = x0 + y

					p6_final = x0[0]
					t6_final = x0[1]
					m6_final = x0[2]

					p7_final = x0[3]
					t7_final = x0[4]
					m7_final = x0[5]

					phi_final = theta - x0[9]

					if phi_final < 0:

						phi_final = abs(phi_final)
						direction_final = 'Downward'

					else:

						phi_final = abs(phi_final)
						direction_final = 'Upward'

					self.p2ValueLabel.configure(text = str(float(p2)))
					self.t2ValueLabel.configure(text = str(float(t2)))
					self.m2ValueLabel.configure(text = str(float(m2)))

					self.p3ValueLabel.configure(text = str(float(p3)))
					self.t3ValueLabel.configure(text = str(float(t3)))
					self.m3ValueLabel.configure(text = str(float(m3)))

					self.p4ValueLabel.configure(text = str(float(p4)))
					self.t4ValueLabel.configure(text = str(float(t4)))
					self.m4ValueLabel.configure(text = str(float(m4)))

					self.p5ValueLabel.configure(text = str(float(p5)))
					self.t5ValueLabel.configure(text = str(float(t5)))
					self.m5ValueLabel.configure(text = str(float(m5)))

					self.p6ValueLabel.configure(text = str(float(p6_final)))
					self.t6ValueLabel.configure(text = str(float(t6_final)))
					self.m6ValueLabel.configure(text = str(float(m6_final)))

					self.p7ValueLabel.configure(text = str(float(p7_final)))
					self.t7ValueLabel.configure(text = str(float(t7_final)))
					self.m7ValueLabel.configure(text = str(float(m7_final)))

					self.phiValueLabel.configure(text = str(float(phi_final)))
					self.dirValueLabel.configure(text = str(direction_final))

				guessWindow = Tk()
				guessWindow.title("Guess")

				guessLabel = Label(guessWindow, text = "Enter Guess Values")
				guessLabel.grid(column = 0, row = 0, columnspan = len(guessLabel.cget("text")))

				beta3GuessLabel = Label(guessWindow, text = chr(946)+"3")
				beta3GuessLabel.grid(column = 0, row = 1, columnspan = len(beta3GuessLabel.cget("text")))

				beta3GuessEntry = Entry(guessWindow)
				beta3GuessEntry.grid(column = len(beta3GuessLabel.cget("text")) + 1, row = 1, columnspan = 20)

				beta5GuessLabel = Label(guessWindow, text = chr(946)+"5")
				beta5GuessLabel.grid(column = 0, row = 2, columnspan = len(beta5GuessLabel.cget("text")))

				beta5GuessEntry = Entry(guessWindow)
				beta5GuessEntry.grid(column = len(beta5GuessLabel.cget("text")) + 1, row = 2, columnspan = 20)

				theta3GuessLabel = Label(guessWindow, text = chr(952)+"3")
				theta3GuessLabel.grid(column = 0, row = 3, columnspan = len(theta3GuessLabel.cget("text")))

				theta3GuessEntry = Entry(guessWindow)
				theta3GuessEntry.grid(column = len(theta3GuessLabel.cget("text")) + 1, row = 3, columnspan = 20)

				theta5GuessLabel = Label(guessWindow, text = chr(952)+"5")
				theta5GuessLabel.grid(column = 0, row = 4, columnspan = len(theta5GuessLabel.cget("text")))

				theta5GuessEntry = Entry(guessWindow)
				theta5GuessEntry.grid(column = len(theta5GuessLabel.cget("text")) + 1, row = 4, columnspan = 20)

				guessEntryButton = Button(guessWindow, text = "Guess", command = guess)
				guessEntryButton.grid(column = 0, row = 5, columnspan = 10)

			if p6_case1 > p7_case1 and p6_case2 > p7_case2:

				def guess():

					exp3Guess = float(exp3GuessEntry.get())
					beta5Guess = float(beta5GuessEntry.get())
					phiGuess = float(phiGuessEntry.get())
					theta5Guess = float(theta5GuessEntry.get())

					eqn = ndarray((1,10),dtype = object_)[0]
					syms = ['p6','t6','m6','p7','t7','m7','exp3','beta5','phi','theta5']
					x0 = [p3,t3,m3,p5,t5,m5,30,45,30,30]

					p6,t6,m6,p7,t7,m7,exp3,beta5,phi,theta5 = symbols(syms)
					symvars = [p6,t6,m6,p7,t7,m7,exp3,beta5,phi,theta5]

					subs_dict = dict()

					eqn[0] = p6 - p7
					eqn[1] = theta5 - phi
					eqn[2] = exp3 - phi
					eqn[3] = t(r(theta5)) - ((m5**2)*s(r(2*beta5)) - 2*cot(r(2*beta5)))/(2 + (m5**2)*(gamma + c(r(2*beta5))))
					eqn[4] = p7 - p5*(1 + gamma*(m5*s(r(beta5)))**2)/(1 + gamma*(m7*s(r(beta5 - theta5)))**2)
					eqn[5] = p6 - p3*(((1 + ((gamma - 1)/2)*(m3**2))/(1 + ((gamma - 1)/2)*(m6**2)))**(gamma/(gamma - 1)))
					eqn[6] = t7 - t5*(1 + ((gamma - 1)/2)*(m5*s(r(beta5)))**2)/(1 + ((gamma - 1)/2)*(m7*s(r(beta5 - theta5)))**2)
					eqn[7] = t6 - t3*((1 + ((gamma - 1)/2)*(m3**2))/(1 + ((gamma - 1)/2)*(m6**2)))
					eqn[8] = m7*s(r(beta5 - theta5)) - sq(((2/(gamma - 1)) + (m5*s(r(beta5)))**2)/((2*gamma/(gamma - 1))*((m5*s(r(beta5)))**2) - 1))
					eqn[9] = r(exp3) - ExpansionWave.symomega(m7,gamma) - ExpansionWave.symomega(m5,gamma)

					j = ndarray((10,10),dtype = object_)

					eqn0 = ndarray((1,10),dtype = float64)[0]

					for iterno in range(100):

						for i in range(len(syms)):
							subs_dict[syms[i]] = x0[i]
						
						for k in range(len(eqn)):eqn0[k] = float(eqn[k].evalf(subs = subs_dict))

						for i_ in range(len(eqn)):
							for j_ in range(len(symvars)):
								j[i_][j_] = eqn[i_].diff(symvars[j_]).evalf(subs = subs_dict)


						y = solve(float64(j),float64(eqn0.reshape((10,1)))).reshape(10)
						x0 = x0 + y

					p6_final = x0[0]
					t6_final = x0[1]
					m6_final = x0[2]

					p7_final = x0[3]
					t7_final = x0[4]
					m7_final = x0[5]

					phi_final = x0[8]
					direction_final = 'Downward'

					self.p2ValueLabel.configure(text = str(float(p2)))
					self.t2ValueLabel.configure(text = str(float(t2)))
					self.m2ValueLabel.configure(text = str(float(m2)))

					self.p3ValueLabel.configure(text = str(float(p3)))
					self.t3ValueLabel.configure(text = str(float(t3)))
					self.m3ValueLabel.configure(text = str(float(m3)))

					self.p4ValueLabel.configure(text = str(float(p4)))
					self.t4ValueLabel.configure(text = str(float(t4)))
					self.m4ValueLabel.configure(text = str(float(m4)))

					self.p5ValueLabel.configure(text = str(float(p5)))
					self.t5ValueLabel.configure(text = str(float(t5)))
					self.m5ValueLabel.configure(text = str(float(m5)))

					self.p6ValueLabel.configure(text = str(float(p6_final)))
					self.t6ValueLabel.configure(text = str(float(t6_final)))
					self.m6ValueLabel.configure(text = str(float(m6_final)))

					self.p7ValueLabel.configure(text = str(float(p7_final)))
					self.t7ValueLabel.configure(text = str(float(t7_final)))
					self.m7ValueLabel.configure(text = str(float(m7_final)))

					self.phiValueLabel.configure(text = str(float(phi_final)))
					self.dirValueLabel.configure(text = str(direction_final))

				guessWindow = Tk()
				guessWindow.title("Guess")

				guessLabel = Label(guessWindow, text = "Enter Guess Values")
				guessLabel.grid(column = 0, row = 0, columnspan = len(guessLabel.cget("text")))

				exp3GuessLabel = Label(guessWindow, text = chr(949)+"3")
				exp3GuessLabel.grid(column = 0, row = 1, columnspan = len(exp3GuessLabel.cget("text")))

				exp3GuessEntry = Entry(guessWindow)
				exp3GuessEntry.grid(column = len(exp3GuessLabel.cget("text")) + 1, row = 1, columnspan = 20)

				beta5GuessLabel = Label(guessWindow, text = chr(946)+"5")
				beta5GuessLabel.grid(column = 0, row = 2, columnspan = len(beta5GuessLabel.cget("text")))

				beta5GuessEntry = Entry(guessWindow)
				beta5GuessEntry.grid(column = len(beta5GuessLabel.cget("text")) + 1, row = 2, columnspan = 20)

				phiGuessLabel = Label(guessWindow, text = chr(981))
				phiGuessLabel.grid(column = 0, row = 3, columnspan = len(phiGuessLabel.cget("text")))

				phiGuessEntry = Entry(guessWindow)
				phiGuessEntry.grid(column = len(phiGuessLabel.cget("text")) + 1, row = 3, columnspan = 20)

				theta5GuessLabel = Label(guessWindow, text = chr(952)+"5")
				theta5GuessLabel.grid(column = 0, row = 4, columnspan = len(theta5GuessLabel.cget("text")))

				theta5GuessEntry = Entry(guessWindow)
				theta5GuessEntry.grid(column = len(theta5GuessLabel.cget("text")) + 1, row = 4, columnspan = 20)

				guessEntryButton = Button(guessWindow, text = "Guess", command = guess)
				guessEntryButton.grid(column = 0, row = 5, columnspan = 10)

			if p6_case1 < p7_case1 and p6_case2 < p7_case2:

				def guess():

					beta3Guess = float(beta3GuessEntry.get())
					exp5Guess = float(exp5GuessEntry.get())
					phiGuess = float(phiGuessEntry.get())

					eqn = ndarray((1,10),dtype = object_)[0]
					syms = ['p6','t6','m6','p7','t7','m7','exp5','beta3','phi','theta3']
					x0 = [p3,t3,m3,p5,t5,m5,30,45,30,30]

					p6,t6,m6,p7,t7,m7,exp5,beta3,phi,theta3 = symbols(syms)
					symvars = [p6,t6,m6,p7,t7,m7,exp5,beta3,phi,theta3]

					subs_dict = dict()

					eqn[0] = p6 - p7
					eqn[1] = theta3 - phi
					eqn[2] = exp5 - phi
					eqn[3] = t(r(theta3)) - ((m3**2)*s(r(2*beta3)) - 2*cot(r(2*beta3)))/(2 + (m3**2)*(gamma + c(r(2*beta3))))
					eqn[4] = p6 - p3*(1 + gamma*(m3*s(r(beta3)))**2)/(1 + gamma*(m6*s(r(beta3 - theta3)))**2)
					eqn[5] = p7 - p5*(((1 + ((gamma - 1)/2)*(m5**2))/(1 + ((gamma - 1)/2)*(m7**2)))**(gamma/(gamma - 1)))
					eqn[6] = t6 - t3*(1 + ((gamma - 1)/2)*(m3*s(r(beta3)))**2)/(1 + ((gamma - 1)/2)*(m6*s(r(beta3 - theta3)))**2)
					eqn[7] = t7 - t5*((1 + ((gamma - 1)/2)*(m5**2))/(1 + ((gamma - 1)/2)*(m7**2)))
					eqn[8] = m6*s(r(beta3 - theta3)) - sq(((2/(gamma - 1)) + (m3*s(r(beta3)))**2)/((2*gamma/(gamma - 1))*((m3*s(r(beta3)))**2) - 1))
					eqn[9] = r(exp5) - ExpansionWave.symomega(m6,gamma) - ExpansionWave.symomega(m3,gamma)

					eqn0 = ndarray((1,10),dtype = float64)[0]

					j = ndarray((10,10),dtype = object_)

					for iterno in range(100):

						for i in range(len(syms)):
							subs_dict[syms[i]] = x0[i]

						for k in range(len(eqn)):eqn0[k] = float(eqn[k].evalf(subs = subs_dict))

						for i_ in range(len(eqn)):
							for j_ in range(len(symvars)):
								j[i_][j_] = eqn[i_].diff(symvars[j_]).evalf(subs = subs_dict)


						y = solve(float64(j),float64(eqn0.reshape((10,1)))).reshape(10)
						x0 = x0 + y

					p6_final = x0[0]
					t6_final = x0[1]
					m6_final = x0[2]

					p7_final = x0[3]
					t7_final = x0[4]
					m7_final = x0[5]

					phi_final = x0[8]
					direction_final = 'Upward'

					self.p2ValueLabel.configure(text = str(float(p2)))
					self.t2ValueLabel.configure(text = str(float(t2)))
					self.m2ValueLabel.configure(text = str(float(m2)))

					self.p3ValueLabel.configure(text = str(float(p3)))
					self.t3ValueLabel.configure(text = str(float(t3)))
					self.m3ValueLabel.configure(text = str(float(m3)))

					self.p4ValueLabel.configure(text = str(float(p4)))
					self.t4ValueLabel.configure(text = str(float(t4)))
					self.m4ValueLabel.configure(text = str(float(m4)))

					self.p5ValueLabel.configure(text = str(float(p5)))
					self.t5ValueLabel.configure(text = str(float(t5)))
					self.m5ValueLabel.configure(text = str(float(m5)))

					self.p6ValueLabel.configure(text = str(float(p6_final)))
					self.t6ValueLabel.configure(text = str(float(t6_final)))
					self.m6ValueLabel.configure(text = str(float(m6_final)))

					self.p7ValueLabel.configure(text = str(float(p7_final)))
					self.t7ValueLabel.configure(text = str(float(t7_final)))
					self.m7ValueLabel.configure(text = str(float(m7_final)))

					self.phiValueLabel.configure(text = str(float(phi_final)))
					self.dirValueLabel.configure(text = str(direction_final))

				guessWindow = Tk()
				guessWindow.title("Guess")

				guessLabel = Label(guessWindow, text = "Enter Guess Values")
				guessLabel.grid(column = 0, row = 0, columnspan = len(guessLabel.cget("text")))

				beta3GuessLabel = Label(guessWindow, text = chr(946)+"3")
				beta3GuessLabel.grid(column = 0, row = 1, columnspan = len(beta3GuessLabel.cget("text")))

				beta3GuessEntry = Entry(guessWindow)
				beta3GuessEntry.grid(column = len(beta3GuessLabel.cget("text")) + 1, row = 1, columnspan = 20)

				exp5GuessLabel = Label(guessWindow, text = chr(949)+"5")
				exp5GuessLabel.grid(column = 0, row = 2, columnspan = len(exp5GuessLabel.cget("text")))

				exp5GuessEntry = Entry(guessWindow)
				exp5GuessEntry.grid(column = len(exp5GuessLabel.cget("text")) + 1, row = 2, columnspan = 20)

				phiGuessLabel = Label(guessWindow, text = chr(981))
				phiGuessLabel.grid(column = 0, row = 3, columnspan = len(phiGuessLabel.cget("text")))

				phiGuessEntry = Entry(guessWindow)
				phiGuessEntry.grid(column = len(phiGuessLabel.cget("text")) + 1, row = 3, columnspan = 20)

				theta3GuessLabel = Label(guessWindow, text = chr(952)+"3")
				theta3GuessLabel.grid(column = 0, row = 4, columnspan = len(theta3GuessLabel.cget("text")))

				theta3GuessEntry = Entry(guessWindow)
				theta3GuessEntry.grid(column = len(theta3GuessLabel.cget("text")) + 1, row = 4, columnspan = 20)

				guessEntryButton = Button(guessWindow, text = "Guess", command = guess)
				guessEntryButton.grid(column = 0, row = 5, columnspan = 10)

		self.inputFrame = Frame(master)
		self.p1Label = Label(self.inputFrame, text = "P1:")
		self.t1Label = Label(self.inputFrame, text = "T1:")
		self.m1Label = Label(self.inputFrame, text = "M1:")
		self.thetaLabel = Label(self.inputFrame, text = chr(952))
		self.alphaLabel = Label(self.inputFrame, text = chr(945))
		self.gammaLabel = Label(self.inputFrame, text = chr(947))

		self.p1Entry = Entry(self.inputFrame)
		self.t1Entry = Entry(self.inputFrame)
		self.m1Entry = Entry(self.inputFrame)
		self.thetaEntry = Entry(self.inputFrame)
		self.alphaEntry = Entry(self.inputFrame)
		self.gammaEntry = Entry(self.inputFrame)

		self.simulateButton = Button(self.inputFrame, text = "Calculate", command = simulate)

		self.region2Frame = Frame(master)
		self.region2Frame.configure(width = 300, height = 200)

		self.p2Label = Label(self.region2Frame, text = "P2:")
		self.t2Label = Label(self.region2Frame, text = "T2:")
		self.m2Label = Label(self.region2Frame, text = "M2:")

		self.p2ValueLabel = Label(self.region2Frame)
		self.t2ValueLabel = Label(self.region2Frame)
		self.m2ValueLabel = Label(self.region2Frame)

		self.region3Frame = Frame(master)
		self.region3Frame.configure(width = 300, height = 200)

		self.p3Label = Label(self.region3Frame, text = "P3:")
		self.t3Label = Label(self.region3Frame, text = "T3:")
		self.m3Label = Label(self.region3Frame, text = "M3:")

		self.p3ValueLabel = Label(self.region3Frame)
		self.t3ValueLabel = Label(self.region3Frame)
		self.m3ValueLabel = Label(self.region3Frame)

		self.region4Frame = Frame(master)
		self.region4Frame.configure(width = 300, height = 200)

		self.p4Label = Label(self.region4Frame, text = "P4:")
		self.t4Label = Label(self.region4Frame, text = "T4:")
		self.m4Label = Label(self.region4Frame, text = "M4:")

		self.p4ValueLabel = Label(self.region4Frame)
		self.t4ValueLabel = Label(self.region4Frame)
		self.m4ValueLabel = Label(self.region4Frame)

		self.region5Frame = Frame(master)
		self.region5Frame.configure(width = 300, height = 200)

		self.p5Label = Label(self.region5Frame, text = "P5:")
		self.t5Label = Label(self.region5Frame, text = "T5:")
		self.m5Label = Label(self.region5Frame, text = "M5:")

		self.p5ValueLabel = Label(self.region5Frame)
		self.t5ValueLabel = Label(self.region5Frame)
		self.m5ValueLabel = Label(self.region5Frame)

		self.region6Frame = Frame(master)
		self.region6Frame.configure(width = 300, height = 200)

		self.p6Label = Label(self.region6Frame, text = "P6:")
		self.t6Label = Label(self.region6Frame, text = "T6:")
		self.m6Label = Label(self.region6Frame, text = "M6:")

		self.p6ValueLabel = Label(self.region6Frame)
		self.t6ValueLabel = Label(self.region6Frame)
		self.m6ValueLabel = Label(self.region6Frame)

		self.region7Frame = Frame(master)
		self.region7Frame.configure(width = 300, height = 200)

		self.p7Label = Label(self.region7Frame, text = "P7:")
		self.t7Label = Label(self.region7Frame, text = "T7:")
		self.m7Label = Label(self.region7Frame, text = "M7:")

		self.p7ValueLabel = Label(self.region7Frame)
		self.t7ValueLabel = Label(self.region7Frame)
		self.m7ValueLabel = Label(self.region7Frame)

		self.miscFrame = Frame(master)
		self.miscFrame.configure(width = 300, height = 200)

		self.phiLabel = Label(self.miscFrame, text = chr(981)+":")
		self.dirLabel = Label(self.miscFrame, text = "Direction:")

		self.phiValueLabel = Label(self.miscFrame)
		self.dirValueLabel = Label(self.miscFrame)

		self.inputFrame.grid(column = 0, row = 0, columnspan = 600)
		self.region2Frame.grid(column = 0, row = 1, columnspan = 300)
		self.region3Frame.grid(column = 301, row = 1, columnspan = 300)
		self.region4Frame.grid(column = 601, row = 1, columnspan = 300)
		self.region5Frame.grid(column = 0, row = 2, columnspan = 300)
		self.region6Frame.grid(column = 301, row = 2, columnspan = 300)
		self.region7Frame.grid(column = 601, row = 2, columnspan = 300)

		self.p1Label.grid(column = 0, row = 0, columnspan = 3)
		self.p1Entry.grid(column = 4, row = 0, columnspan = 20)
		self.t1Label.grid(column = 0, row = 1, columnspan = 3)
		self.t1Entry.grid(column = 4, row = 1, columnspan = 20)
		self.m1Label.grid(column = 0, row = 2, columnspan = 3)
		self.m1Entry.grid(column = 4, row = 2, columnspan = 20)
		self.thetaLabel.grid(column = 0, row = 3, columnspan = 3)
		self.thetaEntry.grid(column = 4, row = 3, columnspan = 20)
		self.alphaLabel.grid(column = 0, row = 4, columnspan = 3)
		self.alphaEntry.grid(column = 4, row = 4, columnspan = 20)
		self.gammaLabel.grid(column = 0, row = 5, columnspan = 3)
		self.gammaEntry.grid(column = 4, row = 5, columnspan = 20)

		self.p2Label.grid(column = 0, row = 0, columnspan = 3)
		self.p2ValueLabel.grid(column = 4, row = 0, columnspan = len(self.p2ValueLabel.cget("text")) + 1)
		self.t2Label.grid(column = 0, row = 1, columnspan = 3)
		self.t2ValueLabel.grid(column = 4, row = 1, columnspan = len(self.t2ValueLabel.cget("text")) + 1)
		self.m2Label.grid(column = 0, row = 2, columnspan = 3)
		self.m2ValueLabel.grid(column = 4, row = 2, columnspan = len(self.m2ValueLabel.cget("text")) + 1)

		self.p3Label.grid(column = 0, row = 0, columnspan = 3)
		self.p3ValueLabel.grid(column = 4, row = 0, columnspan = len(self.p3ValueLabel.cget("text")) + 1)
		self.t3Label.grid(column = 0, row = 1, columnspan = 3)
		self.t3ValueLabel.grid(column = 4, row = 1, columnspan = len(self.t3ValueLabel.cget("text")) + 1)
		self.m3Label.grid(column = 0, row = 2, columnspan = 3)
		self.m3ValueLabel.grid(column = 4, row = 2, columnspan = len(self.m3ValueLabel.cget("text")) + 1)

		self.p4Label.grid(column = 0, row = 0, columnspan = 3)
		self.p4ValueLabel.grid(column = 4, row = 0, columnspan = len(self.p4ValueLabel.cget("text")) + 1)
		self.t4Label.grid(column = 0, row = 1, columnspan = 3)
		self.t4ValueLabel.grid(column = 4, row = 1, columnspan = len(self.t4ValueLabel.cget("text")) + 1)
		self.m4Label.grid(column = 0, row = 2, columnspan = 3)
		self.m4ValueLabel.grid(column = 4, row = 2, columnspan = len(self.m4ValueLabel.cget("text")) + 1)

		self.p5Label.grid(column = 0, row = 0, columnspan = 3)
		self.p5ValueLabel.grid(column = 4, row = 0, columnspan = len(self.p5ValueLabel.cget("text")) + 1)
		self.t5Label.grid(column = 0, row = 1, columnspan = 3)
		self.t5ValueLabel.grid(column = 4, row = 1, columnspan = len(self.t5ValueLabel.cget("text")) + 1)
		self.m5Label.grid(column = 0, row = 2, columnspan = 3)
		self.m5ValueLabel.grid(column = 4, row = 2, columnspan = len(self.m5ValueLabel.cget("text")) + 1)

		self.p6Label.grid(column = 0, row = 0, columnspan = 3)
		self.p6ValueLabel.grid(column = 4, row = 0, columnspan = len(self.p6ValueLabel.cget("text")) + 1)
		self.t6Label.grid(column = 0, row = 1, columnspan = 3)
		self.t6ValueLabel.grid(column = 4, row = 1, columnspan = len(self.t6ValueLabel.cget("text")) + 1)
		self.m6Label.grid(column = 0, row = 2, columnspan = 3)
		self.m6ValueLabel.grid(column = 4, row = 2, columnspan = len(self.m6ValueLabel.cget("text")) + 1)

		self.p7Label.grid(column = 0, row = 0, columnspan = 3)
		self.p7ValueLabel.grid(column = 4, row = 0, columnspan = len(self.p7ValueLabel.cget("text")) + 1)
		self.t7Label.grid(column = 0, row = 1, columnspan = 3)
		self.t7ValueLabel.grid(column = 4, row = 1, columnspan = len(self.t7ValueLabel.cget("text")) + 1)
		self.m7Label.grid(column = 0, row = 2, columnspan = 3)
		self.m7ValueLabel.grid(column = 4, row = 2, columnspan = len(self.m7ValueLabel.cget("text")) + 1)

		self.phiLabel.grid(column = 0, row = 0, columnspan = 10)
		self.dirLabel.grid(column = 0, row = 1, columnspan = 10)

		self.phiValueLabel.grid(column = 10, row = 0, columnspan = 10)
		self.dirValueLabel.grid(column = 10, row = 1, columnspan = 10)

		self.simulateButton.grid(column = 24, row = 5, columnspan = 10)

window = Tk()
SupersonicFlowSimulation(window)
window.mainloop()
