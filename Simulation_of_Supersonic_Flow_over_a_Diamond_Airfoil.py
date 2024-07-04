import math;
import numpy;
import scipy;
import sympy

from tkinter import Tk, Label, Entry, Button, Frame;
from tkinter import messagebox;

class MachValueError(Exception):

        pass;

class GammaValueError(Exception):

        pass;

class PressureValueError(Exception):

        pass;

class TemperatureValueError(Exception):

        pass;

class ThetaValueError(Exception):

        pass;

class AOAError(Exception):

        pass;

class ExpansionRelationErrorRegion5(Exception):

        pass;

class ExpansionRelationErrorRegion3(Exception):

        pass;

class ExpansionRelationErrorRegion6(Exception):

        pass;

class ExpansionRelationErrorRegion7(Exception):

        pass;

class IsentropicFlow:

	def p_sat(p_stat,m_flow,gamma):

		sat_pr = p_stat*((1 + ((gamma - 1)/2)*(m_flow**2))**(gamma/(gamma - 1)));

		return sat_pr;

	def t_sat(t_stat,m_flow,gamma):

		sat_tmp = t_stat*(1 + ((gamma - 1)/2)*(m_flow**2));

		return sat_tmp;

	def d_sat(d_stat,m_flow,gamma):

		sat_d = d_stat*((1 + ((gamma - 1)/2)*(m_flow**2))**(1/(gamma - 1)));

		return sat_d;

class NormalShockWave:

	def m_downstream(m_upstream,gamma):

		exp_1 = (2/(gamma - 1)) + (m_upstream**2);
		exp_2 = ((2*gamma)/(gamma - 1))*(m_upstream**2) - 1;
		mach_downstream = math.sqrt(exp_1/exp_2);

		return mach_downstream;

	def p_ratio(m_upstream,gamma):

		pressure_ratio = (2*gamma/(gamma + 1))*(m_upstream**2) - ((gamma - 1)/(gamma + 1));

		return pressure_ratio;

	def t_ratio(m_upstream,gamma):

		exp_1 = (1 + ((gamma - 1)/2)*(m_upstream**2));
		exp_2 = ((2*gamma)/(gamma - 1))*(m_upstream**2) - 1;
		exp_3 = (((gamma + 1)**2)/(2*(gamma - 1)))*(m_upstream**2);

		temperature_ratio = (exp_1*exp_2)/exp_3;

		return temperature_ratio;

	def d_ratio(pr_ratio, tmp_ratio):

		density_ratio = pr_ratio/tmp_ratio;

		return density_ratio;

def calculate2():

        try:
                calculate()

        except ExpansionRelationErrorRegion5:

                messagebox.showerror("Error!", "The value of Mach Number in Region 4 (M4) is subsonic leading to invalidation of the relation between M4 and M5");

        except ExpansionRelationErrorRegion3:

                messagebox.showerror("Error!", "The value of Mach Number in Region 2 (M2) is subsonic leading to invalidation of the relation between M2 and M3");

        except ExpansionRelationErrorRegion6:

                messagebox.showerror("Error!", "The value of Mach Number in Region 3 (M3) is subsonic leading to invalidation of the relation between M3 and M6");

        except ExpansionRelationErrorRegion7:

                messagebox.showerror("Error!", "The value of Mach Number in Region 5 (M5) is subsonic leading to invalidation of the relation between M5 and M7");

def calculate():

        global zero_case_covered;
        zero_case_covered = 0;

        flag = 0;

        try:
                mach_upstream = float(m_upstream_entry.get());
                gamma = float(gamma_entry.get());
                pressure = float(pressure_entry.get());
                temperature = float(temperature_entry.get());
                theta = float(theta_entry.get());
                aoa = float(aoa_entry.get());

                if mach_upstream <= 1:

                        raise MachValueError;

                if gamma <= 0:

                        raise GammaValueError;

                if pressure < 0:

                        raise PressureValueError;

                if temperature < 0:

                        raise TemperatureValueError;

                if theta < 0:

                        raise ThetaValueError;

                if aoa >= 90 or aoa <= - 90:

                        raise AOAError;

        except MachValueError:

                messagebox.showerror("Error!","Upstream (freestream) Mach Number must be greater than 1");
                flag = 1;

        except GammaValueError:

                messagebox.showerror("Error!","Ratio of specific heats must be greater than zero");
                flag = 1;

        except PressureValueError:

                messagebox.showerror("Error!","Upstream (freestream) pressure value must be non - negative");
                flag = 1;

        except TemperatureValueError:

                messagebox.showerror("Error!","Upstream(freestream) temperature value must be non - negative");
                flag = 1;

        except ThetaValueError:

                messagebox.showerror("Error!","Airfoil Half - Wedge angle must be non - negative");
                flag = 1;

        except AOAError:

                messagebox.showwarning("Warning!", "Angle of Attack must be between 90 and -90 degrees.");
                flag = 1;

        except ValueError:

                messagebox.showerror("Error!","Values entered must be numbers");
                flag = 1;

        finally:

                pass;

        if flag == 0 and (aoa >= 0 and theta > aoa):

                # Region 2 - Upper Leading

                p01 = IsentropicFlow.p_sat(pressure,mach_upstream,gamma);
                t01 = IsentropicFlow.t_sat(temperature,mach_upstream,gamma);

                t02 = t01;

                def theta_beta_m_upper(beta):

                        return ((mach_upstream**2)*math.sin(2*beta*math.pi/180) - 2/math.tan(beta*math.pi/180))/(2 + (mach_upstream**2)*(gamma + math.cos(2*beta*math.pi/180))) - math.tan((theta - aoa)*math.pi/180);

                beta_deg_upper = scipy.optimize.fsolve(theta_beta_m_upper,45)[0];
                mach_upstream_n = mach_upstream*math.sin(math.radians(beta_deg_upper))

                t2 = NormalShockWave.t_ratio(mach_upstream_n,gamma)*temperature;
                m2n = NormalShockWave.m_downstream(mach_upstream_n,gamma);
                m2 = m2n/math.sin((beta_deg_upper - (theta - aoa))*math.pi/180);
                p2 = NormalShockWave.p_ratio(mach_upstream_n,gamma)*pressure;
                p02 = IsentropicFlow.p_sat(p2,m2,gamma);

                # Region 3 - Upper Trailing

                def omega(mach_number):

                        exp_1 = math.sqrt((gamma + 1)/(gamma - 1));
                        exp_2 = math.sqrt((gamma - 1)/(gamma + 1));
                        exp_3 = math.sqrt(mach_number**2 - 1);

                        exp = exp_1*math.atan(exp_2*exp_3) - math.atan(exp_3);

                        return exp;

                def prandtl_meyer_expansion_relation_upper(m):

                        return 2*theta*math.pi/180 + omega(m2) - omega(m)

                if m2 < 1:

                        raise ExpansionRelationErrorRegion3;

                m3 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_upper, m2)[0];
                t03 = t02;

                t3 = t03/(1 + ((gamma - 1)/2)*(m3**2));
                p03 = p02;
                p3 = p03/((1 + ((gamma - 1)/2)*(m3**2))**(gamma/(gamma - 1)));

                # Region 4 - Lower Leading

                def theta_beta_m_lower(beta):

                        return ((mach_upstream**2)*math.sin(2*beta*math.pi/180) - 2/math.tan(beta*math.pi/180))/(2 + (mach_upstream**2)*(gamma + math.cos(2*beta*math.pi/180))) - math.tan((theta + aoa)*math.pi/180);

                beta_deg_lower = scipy.optimize.fsolve(theta_beta_m_lower,45)[0];
                mach_upstream_n = mach_upstream*math.sin(math.radians(beta_deg_lower))

                t4 = NormalShockWave.t_ratio(mach_upstream_n,gamma)*temperature;
                m4n = NormalShockWave.m_downstream(mach_upstream_n,gamma);
                m4 = m4n/math.sin((beta_deg_lower - (theta + aoa))*math.pi/180);
                p4 = NormalShockWave.p_ratio(mach_upstream_n,gamma)*pressure;

                t04 = IsentropicFlow.t_sat(t4,m4,gamma);
                p04 = IsentropicFlow.p_sat(p4,m4,gamma);

                # Region 5 - Lower Trailing

                def prandtl_meyer_expansion_relation_lower(m):

                        return 2*theta*math.pi/180 + omega(m4) - omega(m)

                if m4 < 1:

                        raise ExpansionRelationErrorRegion5;

                m5 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_lower,m4)[0];
                t05 = t04;

                t5 = t05/(1 + ((gamma - 1)/2)*(m5**2));
                p05 = p04;
                p5 = p05/((1 + ((gamma - 1)/2)*(m5**2))**(gamma/(gamma - 1)));

        elif flag == 0 and (aoa >= 0 and theta == aoa):

                t01 = IsentropicFlow.t_sat(temperature,mach_upstream,gamma);
                p01 = IsentropicFlow.p_sat(pressure,mach_upstream,gamma);

                # Region 2 - Upper Leading

                m2 = mach_upstream;
                t02 = t01;
                p02 = p01;

                t2 = t02/(1 + ((gamma - 1)/2)*(m2**2));
                p2 = p02/((1 + ((gamma - 1)/2)*(m2**2))**(gamma/(gamma - 1)));

                # Region 3 - Upper Trailing

                def omega(mach_number):

                        exp_1 = math.sqrt((gamma + 1)/(gamma - 1));
                        exp_2 = math.sqrt((gamma - 1)/(gamma + 1));
                        exp_3 = math.sqrt(mach_number**2 - 1);

                        exp = exp_1*math.atan(exp_2*exp_3) - math.atan(exp_3);

                        return exp;

                def prandtl_meyer_expansion_relation_upper(m):

                        return 2*theta*math.pi/180 + omega(m2) - omega(m)

                if m2 < 1:

                        raise ExpansionRelationErrorRegion3;

                m3 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_upper, m2)[0];
                t03 = t02;

                t3 = t03/(1 + ((gamma - 1)/2)*(m3**2));
                p03 = p02;
                p3 = p03/((1 + ((gamma - 1)/2)*(m3**2))**(gamma/(gamma - 1)));

                # Region 4 - Lower Leading

                def theta_beta_m_lower(beta):

                        return ((mach_upstream**2)*math.sin(2*beta*math.pi/180) - 2/math.tan(beta*math.pi/180))/(2 + (mach_upstream**2)*(gamma + math.cos(2*beta*math.pi/180))) - math.tan((theta + aoa)*math.pi/180);

                beta_deg_lower = scipy.optimize.fsolve(theta_beta_m_lower,45)[0];
                mach_upstream_n = mach_upstream*math.sin(math.radians(beta_deg_lower));

                t4 = NormalShockWave.t_ratio(mach_upstream_n,gamma)*temperature;
                m4n = NormalShockWave.m_downstream(mach_upstream_n,gamma);
                m4 = m4n/math.sin((beta_deg_lower - (theta + aoa))*math.pi/180);
                p4 = NormalShockWave.p_ratio(mach_upstream_n,gamma)*pressure;

                t04 = IsentropicFlow.t_sat(t4,m4,gamma);
                p04 = IsentropicFlow.p_sat(p4,m4,gamma);

                # Region 5 - Lower Trailing

                def prandtl_meyer_expansion_relation_lower(m):

                        return 2*theta*math.pi/180 + omega(m4) - omega(m)

                if m4 < 1:

                        raise ExpansionRelationErrorRegion5;

                m5 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_lower,m4)[0];
                t05 = t04;

                t5 = t05/(1 + ((gamma - 1)/2)*(m5**2));
                p05 = p04;
                p5 = p05/((1 + ((gamma - 1)/2)*(m5**2))**(gamma/(gamma - 1)));

        elif flag == 0 and (aoa >= 0 and theta < aoa):

                t01 = IsentropicFlow.t_sat(temperature,mach_upstream,gamma);
                p01 = IsentropicFlow.p_sat(pressure,mach_upstream,gamma);

                # Region 2 - Upper Leading

                def omega(mach_number):

                        exp_1 = math.sqrt((gamma + 1)/(gamma - 1));
                        exp_2 = math.sqrt((gamma - 1)/(gamma + 1));
                        exp_3 = math.sqrt(mach_number**2 - 1);

                        exp = exp_1*math.atan(exp_2*exp_3) - math.atan(exp_3);

                        return exp;

                def prandtl_meyer_expansion_relation_upper_1(m):

                        return (aoa - theta)*math.pi/180 + omega(mach_upstream) - omega(m);

                m2 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_upper_1,mach_upstream)[0];
                t02 = t01;
                p02 = p01;

                t2 = t02/(1 + ((gamma - 1)/2)*(m2**2));
                p2 = p02/((1 + ((gamma - 1)/2)*(m2**2))**(gamma/(gamma - 1)));
                

                # Region 3 - Upper Trailing

                def prandtl_meyer_expansion_relation_upper_2(m):

                        return 2*(theta)*math.pi/180 + omega(m2) - omega(m);

                if m2 < 1:

                        raise ExpansionRelationErrorRegion3;

                m3 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_upper_2,m2)[0];
                t03 = t02;
                p03 = p02;

                t3 = t03/(1 + ((gamma - 1)/2)*(m3**2));
                p3 = p03/((1 + ((gamma - 1)/2)*(m3**2))**(gamma/(gamma - 1)));

                # Region 4 - Lower Leading

                def theta_beta_m_lower(beta):

                        return ((mach_upstream**2)*math.sin(2*beta*math.pi/180) - 2/math.tan(beta*math.pi/180))/(2 + (mach_upstream**2)*(gamma + math.cos(2*beta*math.pi/180))) - math.tan((theta + aoa)*math.pi/180);

                beta_deg_lower = scipy.optimize.fsolve(theta_beta_m_lower,45)[0];
                mach_upstream_n = mach_upstream*math.sin(math.radians(beta_deg_lower))

                t4 = NormalShockWave.t_ratio(mach_upstream_n,gamma)*temperature;
                m4n = NormalShockWave.m_downstream(mach_upstream_n,gamma);
                m4 = m4n/math.sin((beta_deg_lower - (theta + aoa))*math.pi/180);
                p4 = NormalShockWave.p_ratio(mach_upstream_n,gamma)*pressure;

                t04 = IsentropicFlow.t_sat(t4,m4,gamma);
                p04 = IsentropicFlow.p_sat(p4,m4,gamma);

                # Region 5 - Lower Trailing

                def prandtl_meyer_expansion_relation_lower(m):

                        return 2*theta*math.pi/180 + omega(m4) - omega(m)

                if m4 < 1:

                        raise ExpansionRelationErrorRegion5;

                m5 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_lower,m4)[0];
                t05 = t04;

                t5 = t05/(1 + ((gamma - 1)/2)*(m5**2));
                p05 = p04;
                p5 = p05/((1 + ((gamma - 1)/2)*(m5**2))**(gamma/(gamma - 1)));

        elif flag == 0 and (aoa < 0 and math.fabs(theta) > math.fabs(aoa)):

                # Region 2 - Upper Leading

                p01 = IsentropicFlow.p_sat(pressure,mach_upstream,gamma);
                t01 = IsentropicFlow.t_sat(temperature,mach_upstream,gamma);

                t02 = t01;

                def theta_beta_m_upper(beta):

                        return ((mach_upstream**2)*math.sin(2*beta*math.pi/180) - 2/math.tan(beta*math.pi/180))/(2 + (mach_upstream**2)*(gamma + math.cos(2*beta*math.pi/180))) - math.tan((theta + math.fabs(aoa))*math.pi/180);

                beta_deg_upper = scipy.optimize.fsolve(theta_beta_m_upper,45)[0];
                mach_upstream_n = mach_upstream*math.sin(math.radians(beta_deg_upper))

                t2 = NormalShockWave.t_ratio(mach_upstream_n,gamma)*temperature;
                m2n = NormalShockWave.m_downstream(mach_upstream_n,gamma);
                m2 = m2n/math.sin((beta_deg_upper - (theta - aoa))*math.pi/180);
                p2 = NormalShockWave.p_ratio(mach_upstream_n,gamma)*pressure;
                p02 = IsentropicFlow.p_sat(p2,m2,gamma);

                # Region 3 - Upper Trailing

                def omega(mach_number):

                        exp_1 = math.sqrt((gamma + 1)/(gamma - 1));
                        exp_2 = math.sqrt((gamma - 1)/(gamma + 1));
                        exp_3 = math.sqrt(mach_number**2 - 1);

                        exp = exp_1*math.atan(exp_2*exp_3) - math.atan(exp_3);

                        return exp;

                def prandtl_meyer_expansion_relation_upper(m):

                        return 2*theta*math.pi/180 + omega(m2) - omega(m)

                if m2 < 1:

                        raise ExpansionRelationErrorRegion3;

                m3 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_upper, m2)[0];
                t03 = t02;

                t3 = t03/(1 + ((gamma - 1)/2)*(m3**2));
                p03 = p02;
                p3 = p03/((1 + ((gamma - 1)/2)*(m3**2))**(gamma/(gamma - 1)));

                # Region 4 - Lower Leading

                def theta_beta_m_lower(beta):

                        return ((mach_upstream**2)*math.sin(2*beta*math.pi/180) - 2/math.tan(beta*math.pi/180))/(2 + (mach_upstream**2)*(gamma + math.cos(2*beta*math.pi/180))) - math.tan((theta - math.fabs(aoa))*math.pi/180);

                beta_deg_lower = scipy.optimize.fsolve(theta_beta_m_lower,45)[0];
                mach_upstream_n = mach_upstream*math.sin(math.radians(beta_deg_lower))

                t4 = NormalShockWave.t_ratio(mach_upstream_n,gamma)*temperature;
                m4n = NormalShockWave.m_downstream(mach_upstream_n,gamma);
                m4 = m2n/math.sin((beta_deg_lower - (theta + aoa))*math.pi/180);
                p4 = NormalShockWave.p_ratio(mach_upstream_n,gamma)*pressure;

                t04 = IsentropicFlow.t_sat(t4,m4,gamma);
                p04 = IsentropicFlow.p_sat(p4,m4,gamma);

                # Region 5 - Lower Trailing

                def prandtl_meyer_expansion_relation_lower(m):

                        return 2*theta*math.pi/180 + omega(m4) - omega(m)

                if m4 < 1:

                        raise ExpansionRelationErrorRegion5;

                m5 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_lower,m4)[0];
                t05 = t04;

                t5 = t05/(1 + ((gamma - 1)/2)*(m5**2));
                p05 = p04;
                p5 = p05/((1 + ((gamma - 1)/2)*(m5**2))**(gamma/(gamma - 1)));

        elif flag == 0 and (aoa < 0 and math.fabs(aoa) == math.fabs(theta)):

                # Region 2 - Upper Leading

                p01 = IsentropicFlow.p_sat(pressure,mach_upstream,gamma);
                t01 = IsentropicFlow.t_sat(temperature,mach_upstream,gamma);

                t02 = t01;

                def theta_beta_m_upper(beta):

                        return ((mach_upstream**2)*math.sin(2*beta*math.pi/180) - 2/math.tan(beta*math.pi/180))/(2 + (mach_upstream**2)*(gamma + math.cos(2*beta*math.pi/180))) - math.tan((theta - aoa)*math.pi/180);

                beta_deg_upper = scipy.optimize.fsolve(theta_beta_m_upper,45)[0];
                mach_upstream_n = mach_upstream*math.sin(math.radians(beta_deg_upper))

                t2 = NormalShockWave.t_ratio(mach_upstream_n,gamma)*temperature;
                m2n = NormalShockWave.m_downstream(mach_upstream_n,gamma);
                m2 = m2n/math.sin((beta_deg_upper - (theta - aoa))*math.pi/180);
                p2 = NormalShockWave.p_ratio(mach_upstream_n,gamma)*pressure;
                p02 = IsentropicFlow.p_sat(p2,m2,gamma);

                # Region 3 - Upper Trailing

                def omega(mach_number):

                        exp_1 = math.sqrt((gamma + 1)/(gamma - 1));
                        exp_2 = math.sqrt((gamma - 1)/(gamma + 1));
                        exp_3 = math.sqrt(mach_number**2 - 1);

                        exp = exp_1*math.atan(exp_2*exp_3) - math.atan(exp_3);

                        return exp;

                def prandtl_meyer_expansion_relation_upper(m):

                        return 2*theta*math.pi/180 + omega(m2) - omega(m)

                if m2 < 1:

                        raise ExpansionRelationErrorRegion3;

                m3 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_upper, m2)[0];
                t03 = t02;

                t3 = t03/(1 + ((gamma - 1)/2)*(m3**2));
                p03 = p02;
                p3 = p03/((1 + ((gamma - 1)/2)*(m3**2))**(gamma/(gamma - 1)));

                # Region 4 - Lower Leading

                m4 = mach_upstream;
                p04 = p01;
                t04 = t01;

                t4 = t04/(1 + ((gamma - 1)/2)*(m4**2))
                p4 = p04/((1 + ((gamma - 1)/2)*(m4**2))**(gamma/(gamma - 1)));

                # Region 5 - Lower Trailing

                def prandtl_meyer_expansion_relation_lower(m):

                        return 2*theta*math.pi/180 + omega(m4) - omega(m)

                if m4 < 1:

                        raise ExpansionRelationRegion5;

                m5 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_lower,m4)[0];
                t05 = t04;

                t5 = t05/(1 + ((gamma - 1)/2)*(m5**2));
                p05 = p04;
                p5 = p05/((1 + ((gamma - 1)/2)*(m5**2))**(gamma/(gamma - 1)));

        elif flag == 0 and (aoa < 0 and math.fabs(theta) < math.fabs(aoa)):

                # Region 2 - Upper Leading

                p01 = IsentropicFlow.p_sat(pressure,mach_upstream,gamma);
                t01 = IsentropicFlow.t_sat(temperature,mach_upstream,gamma);

                t02 = t01;

                def theta_beta_m_upper(beta):

                        return ((mach_upstream**2)*math.sin(2*beta*math.pi/180) - 2/math.tan(beta*math.pi/180))/(2 + (mach_upstream**2)*(gamma + math.cos(2*beta*math.pi/180))) - math.tan((theta - aoa)*math.pi/180);

                beta_deg_upper = scipy.optimize.fsolve(theta_beta_m_upper,45)[0];
                mach_upstream_n = mach_upstream*math.sin(math.radians(beta_deg_upper))

                t2 = NormalShockWave.t_ratio(mach_upstream_n,gamma)*temperature;
                m2n = NormalShockWave.m_downstream(mach_upstream_n,gamma);
                m2 = m2n/math.sin((beta_deg_upper - (theta - aoa))*math.pi/180);
                p2 = NormalShockWave.p_ratio(mach_upstream_n,gamma)*pressure;
                p02 = IsentropicFlow.p_sat(p2,m2,gamma);

                # Region 3 - Upper Trailing

                def omega(mach_number):

                        exp_1 = math.sqrt((gamma + 1)/(gamma - 1));
                        exp_2 = math.sqrt((gamma - 1)/(gamma + 1));
                        exp_3 = math.sqrt(mach_number**2 - 1);

                        exp = exp_1*math.atan(exp_2*exp_3) - math.atan(exp_3);

                        return exp;

                def prandtl_meyer_expansion_relation_upper(m):

                        return 2*theta*math.pi/180 + omega(m2) - omega(m)

                if m2 < 1:

                        raise ExpansionRelationErrorRegion3;

                m3 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_upper, m2)[0];
                t03 = t02;

                t3 = t03/(1 + ((gamma - 1)/2)*(m3**2));
                p03 = p02;
                p3 = p03/((1 + ((gamma - 1)/2)*(m3**2))**(gamma/(gamma - 1)));

                # Region 4 - Lower Leading

                def prandtl_meyer_expansion_relation_lower_1(m):

                        return (math.fabs(aoa) - theta)*math.pi/180 + omega(mach_upstream) - omega(m);

                m4 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_lower_1,mach_upstream)[0];
                t04 = t01;
                p04 = p01;

                t4 = t04/(1 + ((gamma - 1)/2)*(m4**2));
                p4 = p04/((1 + ((gamma - 1)/2)*(m4**2))**(gamma/(gamma - 1)));

                # Region 5 - Lower Trailing

                def prandtl_meyer_expansion_relation_lower_2(m):

                        return 2*theta*math.pi/180 + omega(m4) - omega(m);

                if m4 < 1:

                        raise ExpansionRelationErrorRegion5;

                m5 = scipy.optimize.fsolve(prandtl_meyer_expansion_relation_lower_2,m4)[0];
                t05 = t04;
                p05 = p04;

                t5 = t05/(1 + ((gamma - 1)/2)*(m5**2));
                p5 = p05/((1 + ((gamma - 1)/2)*(m5**2))**(gamma/(gamma - 1)));

        else:

                pass;

        p2 = float(p2);
        p3 = float(p3);
        p4 = float(p4);
        p5 = float(p5);

        t2 = float(t2);
        t3 = float(t3);
        t4 = float(t4);
        t5 = float(t5);

        m2 = float(m2);
        m3 = float(m3);
        m4 = float(m4);
        m5 = float(m5);

        if flag == 0:

                if theta == 0 and aoa == 0:

                        p6_final, t6_final, m6_final = p3, t3, m3;
                        p7_final, t7_final, m7_final = p5, t5, m5;

                        direction_final = 'upward';
                        phi_final = 0;

                if theta == 0 and aoa > 0:

                        def omega(m):

                            exp1 = sympy.sqrt((gamma + 1)/(gamma - 1));
                            exp2 = sympy.sqrt((gamma - 1)/(gamma + 1));
                            exp3 = sympy.sqrt(m**2 - 1);

                            exp = exp1*sympy.atan(exp2*exp3) - sympy.atan(exp3);

                            return exp;

                        x0 = [p3,t3,m3,p5,t5,m5,30,45,30,30];

                        p6,t6,m6 = sympy.symbols(('p6','t6','m6'));
                        p7,t7,m7 = sympy.symbols(('p7','t7','m7'));
                        theta3,beta3,exp5,phi = sympy.symbols(('theta3','beta3','exp5','phi'));

                        eqn1 = p6 - p7;
                        eqn2 = theta3 - phi;
                        eqn3 = exp5 - phi;
                        eqn4 = ((m3**2)*sympy.sin(2*beta3*sympy.pi/180) - 2*sympy.cot(beta3*sympy.pi/180))/(2 + (m3**2)*(gamma + sympy.cos(2*beta3*sympy.pi/180))) - sympy.tan(theta3*sympy.pi/180);

                        m3n = m3*sympy.sin(beta3*sympy.pi/180);
                        m6n = m6*sympy.sin((beta3 - theta3)*sympy.pi/180);

                        eqn5 = p6 - p3*((1 + gamma*(m3n**2))/(1 + gamma*(m6n**2)));
                        eqn6 = p7 - p5*(((1 + ((gamma - 1)/2)*(m5**2))/(1 + ((gamma - 1)/2)*(m7**2)))**(gamma/(gamma - 1)));
                        eqn7 = t6 - t3*((1 + ((gamma - 1)/2)*(m3n**2))/(1 + ((gamma - 1)/2)*(m6n**2)));
                        eqn8 = t7 - t5*((1 + ((gamma - 1)/2)*(m5**2))/(1 + ((gamma - 1)/2)*(m7**2)));
                        eqn9 = m6n - sympy.sqrt(((2/(gamma - 1)) + (m3n**2))/((2*gamma/(gamma - 1))*(m3n**2) - 1));
                        eqn10 = omega(m7) - omega(m5) - exp5*sympy.pi/180;

                        for ele in x0:

                                print(type(ele))

                        for iter_n in range(100):

                            subs_dict = {'p6':x0[0],'t6':x0[1],'m6':x0[2],'p7':x0[3],'t7':x0[4],'m7':x0[5],'theta3':x0[6],'beta3':x0[7],'exp5':x0[8],'phi':x0[9]};

                            for key in subs_dict.keys():
                                    subs_dict[key] = float(subs_dict[key])

                            J = numpy.array([[eqn1.diff('p6').evalf(subs = subs_dict),eqn1.diff('t6').evalf(subs = subs_dict),eqn1.diff('m6').evalf(subs = subs_dict),eqn1.diff('p7').evalf(subs = subs_dict),eqn1.diff('t7').evalf(subs = subs_dict),eqn1.diff('m7').evalf(subs = subs_dict),eqn1.diff('theta3').evalf(subs = subs_dict),eqn1.diff('beta3').evalf(subs = subs_dict),eqn1.diff('exp5').evalf(subs = subs_dict),eqn1.diff('phi').evalf(subs = subs_dict)],
                                             [eqn2.diff('p6').evalf(subs = subs_dict),eqn2.diff('t6').evalf(subs = subs_dict),eqn2.diff('m6').evalf(subs = subs_dict),eqn2.diff('p7').evalf(subs = subs_dict),eqn2.diff('t7').evalf(subs = subs_dict),eqn2.diff('m7').evalf(subs = subs_dict),eqn2.diff('theta3').evalf(subs = subs_dict),eqn2.diff('beta3').evalf(subs = subs_dict),eqn2.diff('exp5').evalf(subs = subs_dict),eqn2.diff('phi').evalf(subs = subs_dict)],
                                             [eqn3.diff('p6').evalf(subs = subs_dict),eqn3.diff('t6').evalf(subs = subs_dict),eqn3.diff('m6').evalf(subs = subs_dict),eqn3.diff('p7').evalf(subs = subs_dict),eqn3.diff('t7').evalf(subs = subs_dict),eqn3.diff('m7').evalf(subs = subs_dict),eqn3.diff('theta3').evalf(subs = subs_dict),eqn3.diff('beta3').evalf(subs = subs_dict),eqn3.diff('exp5').evalf(subs = subs_dict),eqn3.diff('phi').evalf(subs = subs_dict)],
                                             [eqn4.diff('p6').evalf(subs = subs_dict),eqn4.diff('t6').evalf(subs = subs_dict),eqn4.diff('m6').evalf(subs = subs_dict),eqn4.diff('p7').evalf(subs = subs_dict),eqn4.diff('t7').evalf(subs = subs_dict),eqn4.diff('m7').evalf(subs = subs_dict),eqn4.diff('theta3').evalf(subs = subs_dict),eqn4.diff('beta3').evalf(subs = subs_dict),eqn4.diff('exp5').evalf(subs = subs_dict),eqn4.diff('phi').evalf(subs = subs_dict)],
                                             [eqn5.diff('p6').evalf(subs = subs_dict),eqn5.diff('t6').evalf(subs = subs_dict),eqn5.diff('m6').evalf(subs = subs_dict),eqn5.diff('p7').evalf(subs = subs_dict),eqn5.diff('t7').evalf(subs = subs_dict),eqn5.diff('m7').evalf(subs = subs_dict),eqn5.diff('theta3').evalf(subs = subs_dict),eqn5.diff('beta3').evalf(subs = subs_dict),eqn5.diff('exp5').evalf(subs = subs_dict),eqn5.diff('phi').evalf(subs = subs_dict)],
                                             [eqn6.diff('p6').evalf(subs = subs_dict),eqn6.diff('t6').evalf(subs = subs_dict),eqn6.diff('m6').evalf(subs = subs_dict),eqn6.diff('p7').evalf(subs = subs_dict),eqn6.diff('t7').evalf(subs = subs_dict),eqn6.diff('m7').evalf(subs = subs_dict),eqn6.diff('theta3').evalf(subs = subs_dict),eqn6.diff('beta3').evalf(subs = subs_dict),eqn6.diff('exp5').evalf(subs = subs_dict),eqn6.diff('phi').evalf(subs = subs_dict)],
                                             [eqn7.diff('p6').evalf(subs = subs_dict),eqn7.diff('t6').evalf(subs = subs_dict),eqn7.diff('m6').evalf(subs = subs_dict),eqn7.diff('p7').evalf(subs = subs_dict),eqn7.diff('t7').evalf(subs = subs_dict),eqn7.diff('m7').evalf(subs = subs_dict),eqn7.diff('theta3').evalf(subs = subs_dict),eqn7.diff('beta3').evalf(subs = subs_dict),eqn7.diff('exp5').evalf(subs = subs_dict),eqn7.diff('phi').evalf(subs = subs_dict)],
                                             [eqn8.diff('p6').evalf(subs = subs_dict),eqn8.diff('t6').evalf(subs = subs_dict),eqn8.diff('m6').evalf(subs = subs_dict),eqn8.diff('p7').evalf(subs = subs_dict),eqn8.diff('t7').evalf(subs = subs_dict),eqn8.diff('m7').evalf(subs = subs_dict),eqn8.diff('theta3').evalf(subs = subs_dict),eqn8.diff('beta3').evalf(subs = subs_dict),eqn8.diff('exp5').evalf(subs = subs_dict),eqn8.diff('phi').evalf(subs = subs_dict)],
                                             [eqn9.diff('p6').evalf(subs = subs_dict),eqn9.diff('t6').evalf(subs = subs_dict),eqn9.diff('m6').evalf(subs = subs_dict),eqn9.diff('p7').evalf(subs = subs_dict),eqn9.diff('t7').evalf(subs = subs_dict),eqn9.diff('m7').evalf(subs = subs_dict),eqn9.diff('theta3').evalf(subs = subs_dict),eqn9.diff('beta3').evalf(subs = subs_dict),eqn9.diff('exp5').evalf(subs = subs_dict),eqn9.diff('phi').evalf(subs = subs_dict)],
                                             [eqn10.diff('p6').evalf(subs = subs_dict),eqn10.diff('t6').evalf(subs = subs_dict),eqn10.diff('m6').evalf(subs = subs_dict),eqn10.diff('p7').evalf(subs = subs_dict),eqn10.diff('t7').evalf(subs = subs_dict),eqn10.diff('m7').evalf(subs = subs_dict),eqn10.diff('theta3').evalf(subs = subs_dict),eqn10.diff('beta3').evalf(subs = subs_dict),eqn10.diff('exp5').evalf(subs = subs_dict),eqn10.diff('phi').evalf(subs = subs_dict)]] ,dtype = numpy.float64);

                            F = numpy.array([[eqn1.evalf(subs = subs_dict)],
                                             [eqn2.evalf(subs = subs_dict)],
                                             [eqn3.evalf(subs = subs_dict)],
                                             [eqn4.evalf(subs = subs_dict)],
                                             [eqn5.evalf(subs = subs_dict)],
                                             [eqn6.evalf(subs = subs_dict)],
                                             [eqn7.evalf(subs = subs_dict)],
                                             [eqn8.evalf(subs = subs_dict)],
                                             [eqn9.evalf(subs = subs_dict)],
                                             [eqn10.evalf(subs = subs_dict)]], dtype = numpy.float64);

                            y = numpy.linalg.solve(J,-F);
                            y = y.transpose()[0];
                            x0 = x0 + y;

                            p6_final = x0[0];
                            t6_final = x0[1];
                            m6_final = x0[2];
                            p7_final = x0[3];
                            t7_final = x0[4];
                            m7_final = x0[5];
                            direction_final = 'upward';
                            phi_final = x0[9];

                if theta == 0 and aoa < 0:

                        def omega(m):

                            exp1 = sympy.sqrt((gamma + 1)/(gamma - 1));
                            exp2 = sympy.sqrt((gamma - 1)/(gamma + 1));
                            exp3 = sympy.sqrt(m**2 - 1);

                            exp = exp1*sympy.atan(exp2*exp3) - sympy.atan(exp3);

                            return exp;

                        x0 = [p3,t3,m3,p5,t5,m5,30,45,30,30];

                        p6,t6,m6 = sympy.symbols(('p6','t6','m6'));
                        p7,t7,m7 = sympy.symbols(('p7','t7','m7'));
                        theta5,beta5,exp3,phi = sympy.symbols(('theta5','beta5','exp3','phi'));

                        eqn1 = p6 - p7;
                        eqn2 = theta5 - phi;
                        eqn3 = exp3 - phi;
                        eqn4 = ((m5**2)*sympy.sin(2*beta5*sympy.pi/180) - 2*sympy.cot(beta5*sympy.pi/180))/(2 + (m5**2)*(gamma + sympy.cos(2*beta5*sympy.pi/180))) - sympy.tan(theta5*sympy.pi/180);

                        m5n = m5*sympy.sin(beta5*sympy.pi/180);
                        m7n = m7*sympy.sin((beta5 - theta5)*sympy.pi/180);

                        eqn5 = p7 - p5*((1 + gamma*(m5n**2))/(1 + gamma*(m7n**2)));
                        eqn6 = p6 - p3*(((1 + ((gamma - 1)/2)*(m3**2))/(1 + ((gamma - 1)/2)*(m6**2)))**(gamma/(gamma - 1)));
                        eqn7 = t7 - t5*((1 + ((gamma - 1)/2)*(m5n**2))/(1 + ((gamma - 1)/2)*(m7n**2)));
                        eqn8 = t6 - t3*((1 + ((gamma - 1)/2)*(m3**2))/(1 + ((gamma - 1)/2)*(m6**2)));
                        eqn9 = m7n - sympy.sqrt(((2/(gamma - 1)) + (m5n**2))/((2*gamma/(gamma - 1))*(m5n**2) - 1));
                        eqn10 = omega(m6) - omega(m3) - exp3*sympy.pi/180;

                        for iter_n in range(100):

                            subs_dict = {'p6':x0[0],'t6':x0[1],'m6':x0[2],'p7':x0[3],'t7':x0[4],'m7':x0[5],'theta5':x0[6],'beta5':x0[7],'exp3':x0[8],'phi':x0[9]};

                            for key in subs_dict.keys():
                                    subs_dict[key] = float(subs_dict[key])

                            J = numpy.array([[eqn1.diff('p6').evalf(subs = subs_dict),eqn1.diff('t6').evalf(subs = subs_dict),eqn1.diff('m6').evalf(subs = subs_dict),eqn1.diff('p7').evalf(subs = subs_dict),eqn1.diff('t7').evalf(subs = subs_dict),eqn1.diff('m7').evalf(subs = subs_dict),eqn1.diff('theta5').evalf(subs = subs_dict),eqn1.diff('beta5').evalf(subs = subs_dict),eqn1.diff('exp3').evalf(subs = subs_dict),eqn1.diff('phi').evalf(subs = subs_dict)],
                                             [eqn2.diff('p6').evalf(subs = subs_dict),eqn2.diff('t6').evalf(subs = subs_dict),eqn2.diff('m6').evalf(subs = subs_dict),eqn2.diff('p7').evalf(subs = subs_dict),eqn2.diff('t7').evalf(subs = subs_dict),eqn2.diff('m7').evalf(subs = subs_dict),eqn2.diff('theta5').evalf(subs = subs_dict),eqn2.diff('beta5').evalf(subs = subs_dict),eqn2.diff('exp3').evalf(subs = subs_dict),eqn2.diff('phi').evalf(subs = subs_dict)],
                                             [eqn3.diff('p6').evalf(subs = subs_dict),eqn3.diff('t6').evalf(subs = subs_dict),eqn3.diff('m6').evalf(subs = subs_dict),eqn3.diff('p7').evalf(subs = subs_dict),eqn3.diff('t7').evalf(subs = subs_dict),eqn3.diff('m7').evalf(subs = subs_dict),eqn3.diff('theta5').evalf(subs = subs_dict),eqn3.diff('beta5').evalf(subs = subs_dict),eqn3.diff('exp3').evalf(subs = subs_dict),eqn3.diff('phi').evalf(subs = subs_dict)],
                                             [eqn4.diff('p6').evalf(subs = subs_dict),eqn4.diff('t6').evalf(subs = subs_dict),eqn4.diff('m6').evalf(subs = subs_dict),eqn4.diff('p7').evalf(subs = subs_dict),eqn4.diff('t7').evalf(subs = subs_dict),eqn4.diff('m7').evalf(subs = subs_dict),eqn4.diff('theta5').evalf(subs = subs_dict),eqn4.diff('beta5').evalf(subs = subs_dict),eqn4.diff('exp3').evalf(subs = subs_dict),eqn4.diff('phi').evalf(subs = subs_dict)],
                                             [eqn5.diff('p6').evalf(subs = subs_dict),eqn5.diff('t6').evalf(subs = subs_dict),eqn5.diff('m6').evalf(subs = subs_dict),eqn5.diff('p7').evalf(subs = subs_dict),eqn5.diff('t7').evalf(subs = subs_dict),eqn5.diff('m7').evalf(subs = subs_dict),eqn5.diff('theta5').evalf(subs = subs_dict),eqn5.diff('beta5').evalf(subs = subs_dict),eqn5.diff('exp3').evalf(subs = subs_dict),eqn5.diff('phi').evalf(subs = subs_dict)],
                                             [eqn6.diff('p6').evalf(subs = subs_dict),eqn6.diff('t6').evalf(subs = subs_dict),eqn6.diff('m6').evalf(subs = subs_dict),eqn6.diff('p7').evalf(subs = subs_dict),eqn6.diff('t7').evalf(subs = subs_dict),eqn6.diff('m7').evalf(subs = subs_dict),eqn6.diff('theta5').evalf(subs = subs_dict),eqn6.diff('beta5').evalf(subs = subs_dict),eqn6.diff('exp3').evalf(subs = subs_dict),eqn6.diff('phi').evalf(subs = subs_dict)],
                                             [eqn7.diff('p6').evalf(subs = subs_dict),eqn7.diff('t6').evalf(subs = subs_dict),eqn7.diff('m6').evalf(subs = subs_dict),eqn7.diff('p7').evalf(subs = subs_dict),eqn7.diff('t7').evalf(subs = subs_dict),eqn7.diff('m7').evalf(subs = subs_dict),eqn7.diff('theta5').evalf(subs = subs_dict),eqn7.diff('beta5').evalf(subs = subs_dict),eqn7.diff('exp3').evalf(subs = subs_dict),eqn7.diff('phi').evalf(subs = subs_dict)],
                                             [eqn8.diff('p6').evalf(subs = subs_dict),eqn8.diff('t6').evalf(subs = subs_dict),eqn8.diff('m6').evalf(subs = subs_dict),eqn8.diff('p7').evalf(subs = subs_dict),eqn8.diff('t7').evalf(subs = subs_dict),eqn8.diff('m7').evalf(subs = subs_dict),eqn8.diff('theta5').evalf(subs = subs_dict),eqn8.diff('beta5').evalf(subs = subs_dict),eqn8.diff('exp3').evalf(subs = subs_dict),eqn8.diff('phi').evalf(subs = subs_dict)],
                                             [eqn9.diff('p6').evalf(subs = subs_dict),eqn9.diff('t6').evalf(subs = subs_dict),eqn9.diff('m6').evalf(subs = subs_dict),eqn9.diff('p7').evalf(subs = subs_dict),eqn9.diff('t7').evalf(subs = subs_dict),eqn9.diff('m7').evalf(subs = subs_dict),eqn9.diff('theta5').evalf(subs = subs_dict),eqn9.diff('beta5').evalf(subs = subs_dict),eqn9.diff('exp3').evalf(subs = subs_dict),eqn9.diff('phi').evalf(subs = subs_dict)],
                                             [eqn10.diff('p6').evalf(subs = subs_dict),eqn10.diff('t6').evalf(subs = subs_dict),eqn10.diff('m6').evalf(subs = subs_dict),eqn10.diff('p7').evalf(subs = subs_dict),eqn10.diff('t7').evalf(subs = subs_dict),eqn10.diff('m7').evalf(subs = subs_dict),eqn10.diff('theta5').evalf(subs = subs_dict),eqn10.diff('beta5').evalf(subs = subs_dict),eqn10.diff('exp3').evalf(subs = subs_dict),eqn10.diff('phi').evalf(subs = subs_dict)]] ,dtype = numpy.float64);

                            F = numpy.array([[eqn1.evalf(subs = subs_dict)],
                                             [eqn2.evalf(subs = subs_dict)],
                                             [eqn3.evalf(subs = subs_dict)],
                                             [eqn4.evalf(subs = subs_dict)],
                                             [eqn5.evalf(subs = subs_dict)],
                                             [eqn6.evalf(subs = subs_dict)],
                                             [eqn7.evalf(subs = subs_dict)],
                                             [eqn8.evalf(subs = subs_dict)],
                                             [eqn9.evalf(subs = subs_dict)],
                                             [eqn10.evalf(subs = subs_dict)]], dtype = numpy.float64);

                            y = numpy.linalg.solve(J,-F);
                            y = y.transpose()[0];
                            x0 = x0 + y;

                            p6_final = x0[0];
                            t6_final = x0[1];
                            m6_final = x0[2];
                            p7_final = x0[3];
                            t7_final = x0[4];
                            m7_final = x0[5];
                            direction_final = 'downward';
                            phi_final = x0[9];

                if theta != 0 and aoa == 0:

                        def theta_beta_m_3_6(beta):

                                num = (m3**2)*math.sin(beta*math.pi/180) - 2/math.tan(beta*math.pi/180);
                                den = 2 + (m3**2)*(gamma + math.cos(2*beta*math.pi/180));

                                exp = (num/den) - math.tan(theta*math.pi/180);

                                return exp;

                        beta_3_6 = scipy.optimize.fsolve(theta_beta_m_3_6,45)[0];
                        m3n = m3*math.sin(beta_3_6*math.pi/180);

                        m6n = NormalShockWave.m_downstream(m3n,gamma);
                        m6 = m6n/math.sin((beta_3_6 - theta)*math.pi/180);
                        p6 = NormalShockWave.p_ratio(m3n,gamma)*p3;
                        t6 = NormalShockWave.t_ratio(m3n,gamma)*t3;

                        p6_final = p6;
                        t6_final = t6;
                        m6_final = m6;
                        p7_final = p6;
                        t7_final = t6;
                        m7_final = m6;
                        direction_final = 'upward';
                        phi_final = 0;

                if theta != 0 and aoa != 0:

                        # Case 1 - Phi = Theta in upward direction

                        def theta_beta_m_3_6(beta):

                                num = (m3**2)*math.sin(beta*math.pi/180) - 2/math.tan(beta*math.pi/180);
                                den = 2 + (m3**2)*(gamma + math.cos(2*beta*math.pi/180));

                                exp = (num/den) - math.tan(2*theta*math.pi/180);

                                return exp;

                        beta_3_6 = scipy.optimize.fsolve(theta_beta_m_3_6,45)[0];
                        m3n = m3*math.sin(beta_3_6*math.pi/180);

                        m6n = NormalShockWave.m_downstream(m3n,gamma);
                        m6_case1 = m6n/math.sin((beta_3_6 - 2*theta)*math.pi/180);
                        p6_case1 = NormalShockWave.p_ratio(m3n,gamma)*p3;
                        t6_case1 = NormalShockWave.t_ratio(m3n,gamma)*t3;

                        p7_case1, t7_case1, m7_case1 = p5, t5, m5;

                        # Case 2 - Phi = Theta in downward direction

                        def theta_beta_m_5_7(beta):

                                num = (m5**2)*math.sin(beta*math.pi/180) - 2/math.tan(beta*math.pi/180);
                                den = 2 + (m5**2)*(gamma + math.cos(2*beta*math.pi/180));

                                exp = (num/den) - math.tan(2*theta*math.pi/180);

                                return exp;

                        beta_5_7 = scipy.optimize.fsolve(theta_beta_m_5_7,45)[0];
                        m5n = m5*math.sin(beta_5_7*math.pi/180);

                        m7n = NormalShockWave.m_downstream(m5n,gamma);
                        m7_case2 = m7n/math.sin((beta_5_7 - 2*theta)*math.pi/180);
                        p7_case2 = NormalShockWave.p_ratio(m5n,gamma)*p5;
                        t7_case2 = NormalShockWave.t_ratio(m5n,gamma)*t5;

                        p6_case2, t6_case2, m6_case2 = p3, t3, m3;

                        if p6_case1 == p7_case1:

                                p6_final, t6_final, m6_final = p6_case1, t6_case1, m6_case1;
                                p7_final, t7_final, m7_final = p7_case1, t7_case1, m7_case1;

                                direction_final = 'upward';
                                phi_final = theta;

                        if p6_case2 == p7_case2:

                                p6_final, t6_final, m6_final = p6_case2, t6_case2, m6_case2;
                                p7_final, t7_final, m7_final = p7_case2, t7_case2, m7_case2;

                                direction_final = 'downward';
                                phi_final = theta;

                        if p6_case1 > p7_case1 and p6_case2 > p7_case2:

                                def omega(m):

                                        exp1 = sympy.sqrt((gamma + 1)/(gamma - 1));
                                        exp2 = sympy.sqrt((gamma - 1)/(gamma + 1));
                                        exp3 = sympy.sqrt(m**2 - 1);

                                        exp = exp1*sympy.atan(exp2*exp3) - sympy.atan(exp3);

                                        return exp;

                                p6,t6,m6 = sympy.symbols(('p6','t6','m6'));
                                p7,t7,m7 = sympy.symbols(('p7','t7','m7'));
                                theta5,beta5,exp3,phi = sympy.symbols(('theta5','beta5','exp3','phi'));

                                eqn1 = p6 - p7;
                                eqn2 = theta5 - phi - theta;
                                eqn3 = exp3 - phi - theta;
                                eqn4 = ((m5**2)*sympy.sin(beta5*sympy.pi/180) - 2*sympy.cot(beta5*sympy.pi/180))/(2 + (m5**2)*(gamma + sympy.cos(2*beta5*sympy.pi/180))) - sympy.tan(theta5*sympy.pi/180)

                                m5n = m5*sympy.sin(beta5*sympy.pi/180);
                                m7n = m7*sympy.sin((beta5 - theta5)*sympy.pi/180);

                                eqn5 = p7 - p5*((1 + gamma*(m5n**2))/(1 + gamma*(m7n**2)));
                                eqn6 = p6 - p3*(((1 + ((gamma - 1)/2)*(m3**2))/(1 + ((gamma - 1)/2)*(m6**2)))**(gamma/(gamma - 1)));
                                eqn7 = t7 - t5*((1 + ((gamma - 1)/2)*(m5n**2))/(1 + ((gamma - 1)/2)*(m7n**2)));
                                eqn8 = t6 - t3*((1 + ((gamma - 1)/2)*(m3**2))/(1 + ((gamma - 1)/2)*(m6**2)));
                                eqn9 = sympy.sqrt(((2/(gamma - 1)) + (m5n**2))/((2*gamma/(gamma - 1))*(m5n**2) - 1));
                                eqn10 = omega(m6) - omega(m3) - exp3*sympy.pi/180;

                                x0 = [p3,t3,m3,p5,t5,m5,30,45,30,30];

                                for iter_n in range(100):

                                        subs_dict = {'p6':x0[0],'t6':x0[1],'m6':x0[2],'p7':x0[3],'t7':x0[4],'m7':x0[5],'theta5':x0[6],'beta5':x0[7],'exp3':x0[8],'phi':x0[9]};

                                        for key in subs_dict.keys():
                                                subs_dict[key] = float(subs_dict[key])

                                        J = numpy.array([[eqn1.diff('p6').evalf(subs = subs_dict),eqn1.diff('t6').evalf(subs = subs_dict),eqn1.diff('m6').evalf(subs = subs_dict),eqn1.diff('p7').evalf(subs = subs_dict),eqn1.diff('t7').evalf(subs = subs_dict),eqn1.diff('m7').evalf(subs = subs_dict),eqn1.diff('theta5').evalf(subs = subs_dict),eqn1.diff('beta5').evalf(subs = subs_dict),eqn1.diff('exp3').evalf(subs = subs_dict),eqn1.diff('phi').evalf(subs = subs_dict)],
                                             [eqn2.diff('p6').evalf(subs = subs_dict),eqn2.diff('t6').evalf(subs = subs_dict),eqn2.diff('m6').evalf(subs = subs_dict),eqn2.diff('p7').evalf(subs = subs_dict),eqn2.diff('t7').evalf(subs = subs_dict),eqn2.diff('m7').evalf(subs = subs_dict),eqn2.diff('theta5').evalf(subs = subs_dict),eqn2.diff('beta5').evalf(subs = subs_dict),eqn2.diff('exp3').evalf(subs = subs_dict),eqn2.diff('phi').evalf(subs = subs_dict)],
                                             [eqn3.diff('p6').evalf(subs = subs_dict),eqn3.diff('t6').evalf(subs = subs_dict),eqn3.diff('m6').evalf(subs = subs_dict),eqn3.diff('p7').evalf(subs = subs_dict),eqn3.diff('t7').evalf(subs = subs_dict),eqn3.diff('m7').evalf(subs = subs_dict),eqn3.diff('theta5').evalf(subs = subs_dict),eqn3.diff('beta5').evalf(subs = subs_dict),eqn3.diff('exp3').evalf(subs = subs_dict),eqn3.diff('phi').evalf(subs = subs_dict)],
                                             [eqn4.diff('p6').evalf(subs = subs_dict),eqn4.diff('t6').evalf(subs = subs_dict),eqn4.diff('m6').evalf(subs = subs_dict),eqn4.diff('p7').evalf(subs = subs_dict),eqn4.diff('t7').evalf(subs = subs_dict),eqn4.diff('m7').evalf(subs = subs_dict),eqn4.diff('theta5').evalf(subs = subs_dict),eqn4.diff('beta5').evalf(subs = subs_dict),eqn4.diff('exp3').evalf(subs = subs_dict),eqn4.diff('phi').evalf(subs = subs_dict)],
                                             [eqn5.diff('p6').evalf(subs = subs_dict),eqn5.diff('t6').evalf(subs = subs_dict),eqn5.diff('m6').evalf(subs = subs_dict),eqn5.diff('p7').evalf(subs = subs_dict),eqn5.diff('t7').evalf(subs = subs_dict),eqn5.diff('m7').evalf(subs = subs_dict),eqn5.diff('theta5').evalf(subs = subs_dict),eqn5.diff('beta5').evalf(subs = subs_dict),eqn5.diff('exp3').evalf(subs = subs_dict),eqn5.diff('phi').evalf(subs = subs_dict)],
                                             [eqn6.diff('p6').evalf(subs = subs_dict),eqn6.diff('t6').evalf(subs = subs_dict),eqn6.diff('m6').evalf(subs = subs_dict),eqn6.diff('p7').evalf(subs = subs_dict),eqn6.diff('t7').evalf(subs = subs_dict),eqn6.diff('m7').evalf(subs = subs_dict),eqn6.diff('theta5').evalf(subs = subs_dict),eqn6.diff('beta5').evalf(subs = subs_dict),eqn6.diff('exp3').evalf(subs = subs_dict),eqn6.diff('phi').evalf(subs = subs_dict)],
                                             [eqn7.diff('p6').evalf(subs = subs_dict),eqn7.diff('t6').evalf(subs = subs_dict),eqn7.diff('m6').evalf(subs = subs_dict),eqn7.diff('p7').evalf(subs = subs_dict),eqn7.diff('t7').evalf(subs = subs_dict),eqn7.diff('m7').evalf(subs = subs_dict),eqn7.diff('theta5').evalf(subs = subs_dict),eqn7.diff('beta5').evalf(subs = subs_dict),eqn7.diff('exp3').evalf(subs = subs_dict),eqn7.diff('phi').evalf(subs = subs_dict)],
                                             [eqn8.diff('p6').evalf(subs = subs_dict),eqn8.diff('t6').evalf(subs = subs_dict),eqn8.diff('m6').evalf(subs = subs_dict),eqn8.diff('p7').evalf(subs = subs_dict),eqn8.diff('t7').evalf(subs = subs_dict),eqn8.diff('m7').evalf(subs = subs_dict),eqn8.diff('theta5').evalf(subs = subs_dict),eqn8.diff('beta5').evalf(subs = subs_dict),eqn8.diff('exp3').evalf(subs = subs_dict),eqn8.diff('phi').evalf(subs = subs_dict)],
                                             [eqn9.diff('p6').evalf(subs = subs_dict),eqn9.diff('t6').evalf(subs = subs_dict),eqn9.diff('m6').evalf(subs = subs_dict),eqn9.diff('p7').evalf(subs = subs_dict),eqn9.diff('t7').evalf(subs = subs_dict),eqn9.diff('m7').evalf(subs = subs_dict),eqn9.diff('theta5').evalf(subs = subs_dict),eqn9.diff('beta5').evalf(subs = subs_dict),eqn9.diff('exp3').evalf(subs = subs_dict),eqn9.diff('phi').evalf(subs = subs_dict)],
                                             [eqn10.diff('p6').evalf(subs = subs_dict),eqn10.diff('t6').evalf(subs = subs_dict),eqn10.diff('m6').evalf(subs = subs_dict),eqn10.diff('p7').evalf(subs = subs_dict),eqn10.diff('t7').evalf(subs = subs_dict),eqn10.diff('m7').evalf(subs = subs_dict),eqn10.diff('theta5').evalf(subs = subs_dict),eqn10.diff('beta5').evalf(subs = subs_dict),eqn10.diff('exp3').evalf(subs = subs_dict),eqn10.diff('phi').evalf(subs = subs_dict)]] ,dtype = numpy.float64);

                                        F = numpy.array([[eqn1.evalf(subs = subs_dict)],
                                                     [eqn2.evalf(subs = subs_dict)],
                                                     [eqn3.evalf(subs = subs_dict)],
                                                     [eqn4.evalf(subs = subs_dict)],
                                                     [eqn5.evalf(subs = subs_dict)],
                                                     [eqn6.evalf(subs = subs_dict)],
                                                     [eqn7.evalf(subs = subs_dict)],
                                                     [eqn8.evalf(subs = subs_dict)],
                                                     [eqn9.evalf(subs = subs_dict)],
                                                     [eqn10.evalf(subs = subs_dict)]], dtype = numpy.float64);

                                        y = numpy.linalg.solve(J,-F);
                                        y = y.transpose()[0];
                                        x0 = x0 + y;

                                p6_final = x0[0];
                                t6_final = x0[1];
                                m6_final = x0[2];
                                p7_final = x0[3];
                                t7_final = x0[4];
                                m7_final = x0[5];
                                direction_final = 'downward';
                                phi_final = x0[9];

                        if p6_case1 > p7_case1 and p6_case2 < p7_case2:

                                p6,t6,m6 = sympy.symbols(('p6','t6','m6'));
                                p7,t7,m7 = sympy.symbols(('p7','t7','m7'));
                                beta3,beta5,theta3,theta5 = sympy.symbols(('beta3','beta5','theta3','theta5'));

                                eqn1 = p6 - p7;
                                eqn2 = theta3 + theta5 - 2*theta;
                                eqn3 = ((m3**2)*sympy.sin(beta3*sympy.pi/180) - 2*sympy.cot(beta3*sympy.pi/180))/(2 + (m3**2)*(gamma + sympy.cos(2*beta3*sympy.pi/180))) - sympy.tan(theta3*sympy.pi/180);
                                eqn4 = ((m5**2)*sympy.sin(beta5*sympy.pi/180) - 2*sympy.cot(beta5*sympy.pi/180))/(2 + (m5**2)*(gamma + sympy.cos(2*beta5*sympy.pi/180))) - sympy.tan(theta5*sympy.pi/180);

                                m3n = m3*sympy.sin(beta3*sympy.pi/180);
                                m5n = m5*sympy.sin(beta5*sympy.pi/180);

                                eqn5 = p6 - p3*((2*gamma/(gamma - 1))*(m3n**2) - (gamma - 1)/(gamma + 1));
                                eqn6 = p7 - p5*((2*gamma/(gamma - 1))*(m5n**2) - (gamma - 1)/(gamma + 1));
                                eqn7 = t6 - t3*(((1 + ((gamma - 1)/2)*(m3n**2))*((2*gamma/(gamma - 1))*(m3n**2)))/((((gamma + 1)**2)/(2*(gamma - 1)))*(m3n**2)));
                                eqn8 = t7 - t5*(((1 + ((gamma - 1)/2)*(m5n**2))*((2*gamma/(gamma - 1))*(m5n**2)))/((((gamma + 1)**2)/(2*(gamma - 1)))*(m5n**2)));
                                eqn9 = m6 - sympy.sqrt(((2/(gamma - 1)) + (m3n**2))/((2*gamma/(gamma - 1))*(m3n**2) - 1));
                                eqn10 = m7 - sympy.sqrt(((2/(gamma - 1)) + (m5n**2))/((2*gamma/(gamma - 1))*(m5n**2) - 1));

                                x0 = [p3,t3,m3,p5,t5,m5,45,45,30,30];

                                for iter_n in range(100):

                                        subs_dict = {'p6':x0[0],'t6':x0[1],'m6':x0[2],'p7':x0[3],'t7':x0[4],'m7':x0[5],'beta3':x0[6],'beta5':x0[7],'theta3':x0[8],'theta5':x0[9]};

                                        for key in subs_dict.keys():
                                                subs_dict[key] = float(subs_dict[key])
                                        
                                        J = numpy.array([[eqn1.diff('p6').evalf(subs = subs_dict),eqn1.diff('t6').evalf(subs = subs_dict),eqn1.diff('m6').evalf(subs = subs_dict),eqn1.diff('p7').evalf(subs = subs_dict),eqn1.diff('t7').evalf(subs = subs_dict),eqn1.diff('m7').evalf(subs = subs_dict),eqn1.diff('beta3').evalf(subs = subs_dict),eqn1.diff('beta5').evalf(subs = subs_dict),eqn1.diff('theta3').evalf(subs = subs_dict),eqn1.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn2.diff('p6').evalf(subs = subs_dict),eqn2.diff('t6').evalf(subs = subs_dict),eqn2.diff('m6').evalf(subs = subs_dict),eqn2.diff('p7').evalf(subs = subs_dict),eqn2.diff('t7').evalf(subs = subs_dict),eqn2.diff('m7').evalf(subs = subs_dict),eqn2.diff('beta3').evalf(subs = subs_dict),eqn2.diff('beta5').evalf(subs = subs_dict),eqn2.diff('theta3').evalf(subs = subs_dict),eqn2.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn3.diff('p6').evalf(subs = subs_dict),eqn3.diff('t6').evalf(subs = subs_dict),eqn3.diff('m6').evalf(subs = subs_dict),eqn3.diff('p7').evalf(subs = subs_dict),eqn3.diff('t7').evalf(subs = subs_dict),eqn3.diff('m7').evalf(subs = subs_dict),eqn3.diff('beta3').evalf(subs = subs_dict),eqn3.diff('beta5').evalf(subs = subs_dict),eqn3.diff('theta3').evalf(subs = subs_dict),eqn3.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn4.diff('p6').evalf(subs = subs_dict),eqn4.diff('t6').evalf(subs = subs_dict),eqn4.diff('m6').evalf(subs = subs_dict),eqn4.diff('p7').evalf(subs = subs_dict),eqn4.diff('t7').evalf(subs = subs_dict),eqn4.diff('m7').evalf(subs = subs_dict),eqn4.diff('beta3').evalf(subs = subs_dict),eqn4.diff('beta5').evalf(subs = subs_dict),eqn4.diff('theta3').evalf(subs = subs_dict),eqn4.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn5.diff('p6').evalf(subs = subs_dict),eqn5.diff('t6').evalf(subs = subs_dict),eqn5.diff('m6').evalf(subs = subs_dict),eqn5.diff('p7').evalf(subs = subs_dict),eqn5.diff('t7').evalf(subs = subs_dict),eqn5.diff('m7').evalf(subs = subs_dict),eqn5.diff('beta3').evalf(subs = subs_dict),eqn5.diff('beta5').evalf(subs = subs_dict),eqn5.diff('theta3').evalf(subs = subs_dict),eqn5.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn6.diff('p6').evalf(subs = subs_dict),eqn6.diff('t6').evalf(subs = subs_dict),eqn6.diff('m6').evalf(subs = subs_dict),eqn6.diff('p7').evalf(subs = subs_dict),eqn6.diff('t7').evalf(subs = subs_dict),eqn6.diff('m7').evalf(subs = subs_dict),eqn6.diff('beta3').evalf(subs = subs_dict),eqn6.diff('beta5').evalf(subs = subs_dict),eqn6.diff('theta3').evalf(subs = subs_dict),eqn6.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn7.diff('p6').evalf(subs = subs_dict),eqn7.diff('t6').evalf(subs = subs_dict),eqn7.diff('m6').evalf(subs = subs_dict),eqn7.diff('p7').evalf(subs = subs_dict),eqn7.diff('t7').evalf(subs = subs_dict),eqn7.diff('m7').evalf(subs = subs_dict),eqn7.diff('beta3').evalf(subs = subs_dict),eqn7.diff('beta5').evalf(subs = subs_dict),eqn7.diff('theta3').evalf(subs = subs_dict),eqn7.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn8.diff('p6').evalf(subs = subs_dict),eqn8.diff('t6').evalf(subs = subs_dict),eqn8.diff('m6').evalf(subs = subs_dict),eqn8.diff('p7').evalf(subs = subs_dict),eqn8.diff('t7').evalf(subs = subs_dict),eqn8.diff('m7').evalf(subs = subs_dict),eqn8.diff('beta3').evalf(subs = subs_dict),eqn8.diff('beta5').evalf(subs = subs_dict),eqn8.diff('theta3').evalf(subs = subs_dict),eqn8.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn9.diff('p6').evalf(subs = subs_dict),eqn9.diff('t6').evalf(subs = subs_dict),eqn9.diff('m6').evalf(subs = subs_dict),eqn9.diff('p7').evalf(subs = subs_dict),eqn9.diff('t7').evalf(subs = subs_dict),eqn9.diff('m7').evalf(subs = subs_dict),eqn9.diff('beta3').evalf(subs = subs_dict),eqn9.diff('beta5').evalf(subs = subs_dict),eqn9.diff('theta3').evalf(subs = subs_dict),eqn9.diff('theta5').evalf(subs = subs_dict)],
                                                         [eqn10.diff('p6').evalf(subs = subs_dict),eqn10.diff('t6').evalf(subs = subs_dict),eqn10.diff('m6').evalf(subs = subs_dict),eqn10.diff('p7').evalf(subs = subs_dict),eqn10.diff('t7').evalf(subs = subs_dict),eqn10.diff('m7').evalf(subs = subs_dict),eqn10.diff('beta3').evalf(subs = subs_dict),eqn10.diff('beta5').evalf(subs = subs_dict),eqn10.diff('theta3').evalf(subs = subs_dict),eqn10.diff('theta5').evalf(subs = subs_dict)]] ,dtype = numpy.float64);

                                        F = numpy.array([[float(eqn1.evalf(subs = subs_dict))],
                                                         [float(eqn2.evalf(subs = subs_dict))],
                                                         [float(eqn3.evalf(subs = subs_dict))],
                                                         [float(eqn4.evalf(subs = subs_dict))],
                                                         [float(eqn5.evalf(subs = subs_dict))],
                                                         [float(eqn6.evalf(subs = subs_dict))],
                                                         [float(eqn7.evalf(subs = subs_dict))],
                                                         [float(eqn8.evalf(subs = subs_dict))],
                                                         [float(eqn9.evalf(subs = subs_dict))],
                                                         [float(eqn10.evalf(subs = subs_dict))]], dtype = numpy.float64);

                                        y = numpy.linalg.solve(J,-F);
                                        y = y.transpose()[0];
                                        x0 = x0 + y;

                                        p6_final = x0[0];
                                        t6_final = x0[1];
                                        m6_final = x0[2];
                                        p7_final = x0[3];
                                        t7_final = x0[4];
                                        m7_final = x0[5];

                                        phi_int = x0[9];

                                        if phi_int < 0:

                                                phi_final = numpy.abs(phi_int);
                                                direction_final = 'downward';

                                        if phi_int >= 0:

                                                phi_final = phi_int;
                                                direction_final = 'upward';

                        if p6_case2 < p7_case2 and p6_case1 < p7_case1:

                                def omega(m):

                                    exp1 = sympy.sqrt((gamma + 1)/(gamma - 1));
                                    exp2 = sympy.sqrt((gamma - 1)/(gamma + 1));
                                    exp3 = sympy.sqrt(m**2 - 1);

                                    exp = exp1*sympy.atan(exp2*exp3) - sympy.atan(exp3);

                                    return exp;

                                x0 = [p3,t3,m3,p5,t5,m5,30,45,30,30];

                                p6,t6,m6 = sympy.symbols(('p6','t6','m6'));
                                p7,t7,m7 = sympy.symbols(('p7','t7','m7'));
                                theta3,beta3,exp5,phi = sympy.symbols(('theta3','beta3','exp5','phi'));

                                eqn1 = p6 - p7;
                                eqn2 = theta3 - phi;
                                eqn3 = exp5 - phi;
                                eqn4 = ((m3**2)*sympy.sin(2*beta3*sympy.pi/180) - 2*sympy.cot(beta3*sympy.pi/180))/(2 + (m3**2)*(gamma + sympy.cos(2*beta3*sympy.pi/180))) - sympy.tan(theta3*sympy.pi/180);

                                m3n = m3*sympy.sin(beta3*sympy.pi/180);
                                m6n = m6*sympy.sin((beta3 - theta3)*sympy.pi/180);

                                eqn5 = p6 - p3*((1 + gamma*(m3n**2))/(1 + gamma*(m6n**2)));
                                eqn6 = p7 - p5*(((1 + ((gamma - 1)/2)*(m5**2))/(1 + ((gamma - 1)/2)*(m7**2)))**(gamma/(gamma - 1)));
                                eqn7 = t6 - t3*((1 + ((gamma - 1)/2)*(m3n**2))/(1 + ((gamma - 1)/2)*(m6n**2)));
                                eqn8 = t7 - t5*((1 + ((gamma - 1)/2)*(m5**2))/(1 + ((gamma - 1)/2)*(m7**2)));
                                eqn9 = m6n - sympy.sqrt(((2/(gamma - 1)) + (m3n**2))/((2*gamma/(gamma - 1))*(m3n**2) - 1));
                                eqn10 = omega(m7) - omega(m5) - exp5*sympy.pi/180;

                                for iter_n in range(100):

                                    subs_dict = {'p6':x0[0],'t6':x0[1],'m6':x0[2],'p7':x0[3],'t7':x0[4],'m7':x0[5],'theta3':x0[6],'beta3':x0[7],'exp5':x0[8],'phi':x0[9]};

                                    for key in subs_dict.keys():
                                            subs_dict[key] = float(subs_dict[key])

                                    J = numpy.array([[eqn1.diff('p6').evalf(subs = subs_dict),eqn1.diff('t6').evalf(subs = subs_dict),eqn1.diff('m6').evalf(subs = subs_dict),eqn1.diff('p7').evalf(subs = subs_dict),eqn1.diff('t7').evalf(subs = subs_dict),eqn1.diff('m7').evalf(subs = subs_dict),eqn1.diff('theta3').evalf(subs = subs_dict),eqn1.diff('beta3').evalf(subs = subs_dict),eqn1.diff('exp5').evalf(subs = subs_dict),eqn1.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn2.diff('p6').evalf(subs = subs_dict),eqn2.diff('t6').evalf(subs = subs_dict),eqn2.diff('m6').evalf(subs = subs_dict),eqn2.diff('p7').evalf(subs = subs_dict),eqn2.diff('t7').evalf(subs = subs_dict),eqn2.diff('m7').evalf(subs = subs_dict),eqn2.diff('theta3').evalf(subs = subs_dict),eqn2.diff('beta3').evalf(subs = subs_dict),eqn2.diff('exp5').evalf(subs = subs_dict),eqn2.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn3.diff('p6').evalf(subs = subs_dict),eqn3.diff('t6').evalf(subs = subs_dict),eqn3.diff('m6').evalf(subs = subs_dict),eqn3.diff('p7').evalf(subs = subs_dict),eqn3.diff('t7').evalf(subs = subs_dict),eqn3.diff('m7').evalf(subs = subs_dict),eqn3.diff('theta3').evalf(subs = subs_dict),eqn3.diff('beta3').evalf(subs = subs_dict),eqn3.diff('exp5').evalf(subs = subs_dict),eqn3.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn4.diff('p6').evalf(subs = subs_dict),eqn4.diff('t6').evalf(subs = subs_dict),eqn4.diff('m6').evalf(subs = subs_dict),eqn4.diff('p7').evalf(subs = subs_dict),eqn4.diff('t7').evalf(subs = subs_dict),eqn4.diff('m7').evalf(subs = subs_dict),eqn4.diff('theta3').evalf(subs = subs_dict),eqn4.diff('beta3').evalf(subs = subs_dict),eqn4.diff('exp5').evalf(subs = subs_dict),eqn4.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn5.diff('p6').evalf(subs = subs_dict),eqn5.diff('t6').evalf(subs = subs_dict),eqn5.diff('m6').evalf(subs = subs_dict),eqn5.diff('p7').evalf(subs = subs_dict),eqn5.diff('t7').evalf(subs = subs_dict),eqn5.diff('m7').evalf(subs = subs_dict),eqn5.diff('theta3').evalf(subs = subs_dict),eqn5.diff('beta3').evalf(subs = subs_dict),eqn5.diff('exp5').evalf(subs = subs_dict),eqn5.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn6.diff('p6').evalf(subs = subs_dict),eqn6.diff('t6').evalf(subs = subs_dict),eqn6.diff('m6').evalf(subs = subs_dict),eqn6.diff('p7').evalf(subs = subs_dict),eqn6.diff('t7').evalf(subs = subs_dict),eqn6.diff('m7').evalf(subs = subs_dict),eqn6.diff('theta3').evalf(subs = subs_dict),eqn6.diff('beta3').evalf(subs = subs_dict),eqn6.diff('exp5').evalf(subs = subs_dict),eqn6.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn7.diff('p6').evalf(subs = subs_dict),eqn7.diff('t6').evalf(subs = subs_dict),eqn7.diff('m6').evalf(subs = subs_dict),eqn7.diff('p7').evalf(subs = subs_dict),eqn7.diff('t7').evalf(subs = subs_dict),eqn7.diff('m7').evalf(subs = subs_dict),eqn7.diff('theta3').evalf(subs = subs_dict),eqn7.diff('beta3').evalf(subs = subs_dict),eqn7.diff('exp5').evalf(subs = subs_dict),eqn7.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn8.diff('p6').evalf(subs = subs_dict),eqn8.diff('t6').evalf(subs = subs_dict),eqn8.diff('m6').evalf(subs = subs_dict),eqn8.diff('p7').evalf(subs = subs_dict),eqn8.diff('t7').evalf(subs = subs_dict),eqn8.diff('m7').evalf(subs = subs_dict),eqn8.diff('theta3').evalf(subs = subs_dict),eqn8.diff('beta3').evalf(subs = subs_dict),eqn8.diff('exp5').evalf(subs = subs_dict),eqn8.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn9.diff('p6').evalf(subs = subs_dict),eqn9.diff('t6').evalf(subs = subs_dict),eqn9.diff('m6').evalf(subs = subs_dict),eqn9.diff('p7').evalf(subs = subs_dict),eqn9.diff('t7').evalf(subs = subs_dict),eqn9.diff('m7').evalf(subs = subs_dict),eqn9.diff('theta3').evalf(subs = subs_dict),eqn9.diff('beta3').evalf(subs = subs_dict),eqn9.diff('exp5').evalf(subs = subs_dict),eqn9.diff('phi').evalf(subs = subs_dict)],
                                                     [eqn10.diff('p6').evalf(subs = subs_dict),eqn10.diff('t6').evalf(subs = subs_dict),eqn10.diff('m6').evalf(subs = subs_dict),eqn10.diff('p7').evalf(subs = subs_dict),eqn10.diff('t7').evalf(subs = subs_dict),eqn10.diff('m7').evalf(subs = subs_dict),eqn10.diff('theta3').evalf(subs = subs_dict),eqn10.diff('beta3').evalf(subs = subs_dict),eqn10.diff('exp5').evalf(subs = subs_dict),eqn10.diff('phi').evalf(subs = subs_dict)]] ,dtype = numpy.float64);

                                    F = numpy.array([[eqn1.evalf(subs = subs_dict)],
                                                     [eqn2.evalf(subs = subs_dict)],
                                                     [eqn3.evalf(subs = subs_dict)],
                                                     [eqn4.evalf(subs = subs_dict)],
                                                     [eqn5.evalf(subs = subs_dict)],
                                                     [eqn6.evalf(subs = subs_dict)],
                                                     [eqn7.evalf(subs = subs_dict)],
                                                     [eqn8.evalf(subs = subs_dict)],
                                                     [eqn9.evalf(subs = subs_dict)],
                                                     [eqn10.evalf(subs = subs_dict)]], dtype = numpy.float64);

                                    y = numpy.linalg.solve(J,-F);
                                    y = y.transpose()[0];
                                    x0 = x0 + y;

                                p6_final = x0[0];
                                t6_final = x0[1];
                                m6_final = x0[2];
                                p7_final = x0[3];
                                t7_final = x0[4];
                                m7_final = x0[5];
                                direction_final = 'upward';
                                phi_final = x0[9];
                 
        if flag == 0:

                def show_slip_line_info():

                        info_text = "Key:- \n\n   Region 6 - Region above the slip line \n   Region 7 - Region below the slip line \n\n Values:- \n\n P6: " + str(p6_final) + " \n T6: " + str(t6_final) + " \n M6 :" + str(m6_final) + "\n\n P7: " + str(p7_final) + " \n T7: " + str(t7_final) + " \n M7 :" + str(m7_final) + "\n\n Slip Line Properties: \n\n Direction (with respect to airfoil chord) : " + str(direction_final) + "\n Angle: " + str(phi_final);
                        messagebox.showinfo("Slip Line Info", info_text);

                window2 = Tk();
                window2.title("Simulation of Supersonic Flow over a Diamond Airfoil");
                window.geometry("600x400");

                table_frame = Frame(window2);

                region_2_label = Label(table_frame, text = "Region 2 - Upper Leading");
                region_3_label = Label(table_frame, text = "Region 3 - Upper Trailing");
                region_4_label = Label(table_frame, text = "Region 4 - Lower Leading");
                region_5_label = Label(table_frame, text = "Region 5 - Lower Trailing");

                p_label = Label(table_frame, text = "P");
                t_label = Label(table_frame, text = "T");
                m_label = Label(table_frame, text = "M");

                p2_label = Label(table_frame, text = str(p2));
                t2_label = Label(table_frame, text = str(t2));
                m2_label = Label(table_frame, text = str(m2));

                p3_label = Label(table_frame, text = str(p3));
                t3_label = Label(table_frame, text = str(t3));
                m3_label = Label(table_frame, text = str(m3));

                p4_label = Label(table_frame, text = str(p4));
                t4_label = Label(table_frame, text = str(t4));
                m4_label = Label(table_frame, text = str(m4));

                p5_label = Label(table_frame, text = str(p5));
                t5_label = Label(table_frame, text = str(t5));
                m5_label = Label(table_frame, text = str(m5));

                region_2_label.grid(column = 0, row = 1);
                region_3_label.grid(column = 0, row = 2);
                region_4_label.grid(column = 0, row = 3);
                region_5_label.grid(column = 0, row = 4);

                p_label.grid(column = 1, row = 0);
                t_label.grid(column = 2, row = 0);
                m_label.grid(column = 3, row = 0);

                p2_label.grid(column = 1, row = 1);
                p3_label.grid(column = 1, row = 2);
                p4_label.grid(column = 1, row = 3);
                p5_label.grid(column = 1, row = 4);

                t2_label.grid(column = 2, row = 1);
                t3_label.grid(column = 2, row = 2);
                t4_label.grid(column = 2, row = 3);
                t5_label.grid(column = 2, row = 4);

                m2_label.grid(column = 3, row = 1);
                m3_label.grid(column = 3, row = 2);
                m4_label.grid(column = 3, row = 3);
                m5_label.grid(column = 3, row = 4);

                table_frame.grid(column = 0, row = 0);

                slip_line_info_button = Button(window2, text = "Slip Line Info", command = show_slip_line_info);
                slip_line_info_button.grid(column = 0, row = 1);

window = Tk();
window.title("Simulation of Supersonic Flow over a Diamond Airfoil");
window.geometry("600x400");

m_upstream_label = Label(window, text = "Upstream Mach Number:");
gamma_label = Label(window,text = "Ratio of specific heats:");
pressure_label = Label(window, text = "Upstream (Freestream) Static Pressure (Pa):");
temperature_label = Label(window, text = "Upstream (Freestream) Static Temperature (K):");
theta_label = Label(window, text = "Airfoil Half - Wedge Angle (deg.):");
aoa_label = Label(window,text = "Angle of Attack (deg.):");

m_upstream_entry = Entry(window);
gamma_entry = Entry(window);
pressure_entry = Entry(window);
temperature_entry = Entry(window);
theta_entry = Entry(window);
aoa_entry = Entry(window);

calculate_button = Button(window, text = "Calculate", command = calculate2);

m_upstream_label.grid(row = 0, column = 0, columnspan = len(m_upstream_label.cget('text')))
gamma_label.grid(row = 1, column = 0, columnspan = len(gamma_label.cget('text')))
pressure_label.grid(row = 2, column = 0, columnspan = len(pressure_label.cget('text')))
temperature_label.grid(row = 3, column = 0, columnspan = len(temperature_label.cget('text')))
theta_label.grid(row = 4, column = 0, columnspan = len(theta_label.cget('text')))
aoa_label.grid(row = 5, column = 0, columnspan = len(aoa_label.cget('text')))

m_upstream_entry.grid(row = 0, column = len(m_upstream_label.cget('text')), columnspan = 20)
gamma_entry.grid(row = 1, column = len(gamma_label.cget('text')), columnspan = 20)
pressure_entry.grid(row = 2, column = len(pressure_label.cget('text')), columnspan = 20)
temperature_entry.grid(row = 3, column = len(temperature_label.cget('text')), columnspan = 20)
theta_entry.grid(row = 4, column = len(theta_label.cget('text')), columnspan = 20)
aoa_entry.grid(row = 5, column = len(aoa_label.cget('text')), columnspan = 20)

calculate_button.grid(row = 6, column = 0, columnspan = len(calculate_button.cget('text')))

window.mainloop()
