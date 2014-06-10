package networkModels
import SI = Modelica.SIunits;
constant Real k_g = 1;
constant Real k_o = 10;
constant Real k_w = -1;
constant Real k_fl = -1000;
constant Real k_gl = -1;

package Interfaces
	import SI = Modelica.SIunits;
	connector pressurePortIn = input SI.Pressure(nominal=1e5, min=0);
	connector pressurePortOut = output SI.Pressure(nominal=1e5, min=0);
	connector flowPortIn = input SI.MassFlowRate(nominal=1, min=0, start=1);
	connector threePhaseFlowPortIn = flowPortIn[3];
	connector flowPortOut = output SI.MassFlowRate(nominal=1, min=0, start=1);
	connector threePhaseFlowPortOut = flowPortOut[3];
	
//	partial block singlePhaseIn 
//		flowPortIn w_in(finalstart=0);		
//	end singlePhaseIn;
//	partial block singlePhaseOut
//		flowPortOut w_out(start=0);
//	end singlePhaseOut;
//	partial blcok threePhaseIn;
//		flowPortIn()
//	partial block pressureIn
//		pressurePortIn p_in(start=1e5, nominal=1e5);
//	end pressureIn;
//	partial block pressureOut
//		pressurePortOut p_out(start=1e5, nominal=1e5);	
//	end pressureOut;
//	
	partial block TPSIBlock 
		threePhaseFlowPortIn w_in;
		pressurePortIn p_in;	 	
	end singlePhase; 

	partial block TPMIBlock 
		parameter Integer nin=1;
		threePhaseFlowPortIn w_in[nin];
		pressurePortIn p_in;	 	
	end singlePhase; 
	
	partial block SPSIBlock 
		flowPortIn w_in;
		pressurePortIn p_in;	 	
	end singlePhase; 

	partial block SPMIBlock 
		parameter Integer nin=1;
		flowPortIn w_in[nin];
		pressurePortIn p_in;	 	
	end singlePhase; 
	
	partial block TPMOBlock
		parameter Integer nout=1;
		threePhaseFlowPortOut w_out[nout];
		pressurePortIn p_out;		
	end singlePhase;
	
	partial block TPSOBlock
		threePhaseFlowPortOut w_out;
		pressurePortIn p_out;		
	end singlePhase;
		
	partial block SPSOBlock
		flowPortOut w_out;
		pressurePortIn p_out;		
	end singlePhase;

	partial block SPMOBlock
		parameter Integer nout=1;
		flowPortOut w_out[nout];
		pressurePortIn p_out;		
	end singlePhase;
	//Complete Interfaces
	partial block SPSISOBlock
		extends SPSIBlock;
		extends SPSOBlock;
	end SPSISOBlock;
	
	partial block TPSISOBlock
		extends TPSIBlock;
		extends TPSOBlock;
	end TPSISOBlock;
	
	partial block SPMISOBlock
		extends SPMIBlock;
		extends SPSOBlock;
	end SPMISOBlock;
	
	partial block TPMISOBlock
		extends TPMIBlock;
		extends TPSOBlock;
	end TPMISOBlock;
	
	partial block SPSITPSOBlock
		extends SPSIBlock;
		extends TPSOBlock;
	end SPSITPSOBlock;
end Interfaces;

package networkComponents
function maxFunc  
	input Real x1, x2; 
	output Real y; 
algorithm  
	y := max(x1,x2);
	//y := if (x2> x1) then x2 else x1; 
end maxFunc; 
function heaviside
	input Real x;
	output Real y;
algorithm
	y:=1/(1+exp(-2*1*x));
end heaviside;
	
record wellParameters
	import SI = Modelica.SIunits;

	parameter SI.Temperature T_a=350;          // [K]	// annular temperature
	parameter SI.Temperature T_t=350;          // [K]	// tubing temperature	   
	parameter SI.MolarMass M_g = 0.0195;       // [kg/mol]	// lift gas molecular weight (Molar Mass)  
	parameter SI.Area A_a=0.02;         // [m^2]	// annular cross sectional area	   
	parameter SI.Area A_t=0.012;        // [m^2]	// tubing cross sectional area	   
	parameter SI.Pressure p_r=25e6;         // [Pa]		// reservoir pressure   
	parameter SI.Volume V_a = 30;               // [m^3]		// Annular volume	   
	parameter SI.Volume V_t=  18;               // [m^3]		// Tubing volume	   
	parameter SI.Length L_w=  400;          // [m]		// tubing length from the injection point to the reservoir	   
	parameter SI.VolumeFraction r_wc= 0.4;          // [-]		// well water cut. Volume or Weight?
	parameter Real r_gor=0.08;     // [-]		// well gor	   
	parameter Real Q_max=0.025;   // [m^3/s]		// Empirical constant representing the theoretical absolute open flow (AOF)	   
	parameter SI.Density rho_w = 1030;       // [kg/m^3]		// water density	   
	parameter SI.Density rho_o = 930;        // [kg/m^3]		// oil density	   
	parameter SI.Density rho_Lw = (rho_w * r_wc + rho_o * (1-r_wc));  	//A.4b - liquid density

	parameter Real C=0.8;            // [-]	
	parameter Real C_iv = 0.00016;   // [m^2]	// gas injection valve constant	   
	parameter Real C_pc = 0.0014;    // [m^2]	// production choke valve constant	   
	parameter Real C_gl = 0.00016;
end wellParameters;


model well "Well Model"
	//******************//***********************//
	//               Parameters					 //
	//******************//***********************//
	    extends Interfaces.SPSITPSOBlock; //Single Phase Singe In - Thre Phase Single Out
		import SI = Modelica.SIunits;
		import R = Modelica.Constants.R;
		import g = Modelica.Constants.g_n;
		
		parameter wellParameters p;	
		parameter Real r_glr=(1-p.r_wc)*p.r_gor;		// Gas-to-liquid ratio in well	
	//Initialization	
		parameter SI.Mass m_ga0(nominal=1e3) = 3629.07517367; 
		parameter SI.Mass m_gt0(nominal=1e3) = 1389.20150433;
		parameter SI.Mass m_lt0(nominal=1e3) = 3352.26319015;	
	//State Reference
		parameter SI.Mass m_ga_ref(nominal=1e3) = 0; 
		parameter SI.Mass m_gt_ref(nominal=1e3) = 0;
		parameter SI.Mass m_lt_ref(nominal=1e3) = 0;
		
	//******************//***********************//
	//               States & Variables			 //
	//******************//***********************//	
	//States
		SI.Mass m_ga(start= m_ga0, fixed=true, nominal=1e3, min = 0); // mass of gas in the well annular
		SI.Mass m_gt(start= m_gt0, fixed=true, nominal=1e3, min = 0); // mass of gas in the well tubing
		SI.Mass m_lt(start= m_lt0, fixed=true, nominal=1e3, min = 0); // mass of liquid in the well tubing
		
	//Variables 
		SI.Pressure p_ai(nominal=1e6, min =1e5);
		SI.Pressure p_p	(nominal=1e6, min =1e5);
		SI.Pressure p_ti(nominal=1e6, min =1e5);
		SI.Pressure p_bh(nominal=1e6, min =1e5);
		SI.Pressure p_ta(nominal=1e6, min =1e5);
		
		SI.Density rho_gi(nominal=1e2, min=0);
		SI.Density rho_p (nominal=1e2, min=0);
		SI.Density rho_gl(nominal=1e2, min=0);

		SI.MassFlowRate w_gi(min=0, start=0);		
		SI.MassFlowRate w_p(min=0, start=0);
		SI.MassFlowRate w_lp(min=0, start=0);
		SI.MassFlowRate w_lr(min=0, start=0);
		SI.MassFlowRate w_gr(min=0, start=0);
		SI.MassFlowRate w_gl_max(min=0, start=0);
		
	//******************//***********************//
	//               Equations					 //
	//******************//***********************//			
	equation 
		p_ai = ((R*p.T_a)/(p.V_a*p.M_g)+g/(2*p.A_a))*m_ga;  //A.3a
		p_p = (R*p.T_t*m_gt)/(p.M_g*p.V_t-p.M_g*(1/p.rho_Lw)*m_lt)-(g*(m_gt))/(2*p.A_t); //A.3b
		p_ti = p_p+(g*(m_lt+m_gt))/p.A_t;  //A.3c
		p_bh = ((1+r_glr+(r_glr*p.M_g*g*p.L_w)/(2*R*p.T_t))*p_ti+p.rho_Lw*g*p.L_w)/(1+r_glr-(r_glr*p.M_g*g*p.L_w)/(2*R*p.T_t)); //A.3d
		p_ta = ((R*p.T_a)/(p.V_a*p.M_g)-g/(2*p.A_a))*m_ga;		//A.3e - Pressure in annulus at gas-lift choke valve
		
		rho_gi = (p.M_g/(R*p.T_a)*p_ai); //A.4a
		rho_p = (p.rho_Lw*p.M_g*p_p*(m_lt+m_gt))/(p.rho_Lw*R*p.T_t*m_lt+p.M_g*p_p*m_gt); //A.4b
		rho_gl = p.M_g*p_in/(R*p.T_a); //A.4d
		
		w_gi = p.C_iv*sqrt(rho_gi*maxFunc(0,p_ai-p_ti)); //A.2a
		w_p = p.C_pc*sqrt(rho_p*maxFunc(0, p_p-p_out)); //A.2b
		w_out[1]= (m_gt/(m_lt+m_gt))*w_p; //A.2c
		w_lp = (m_lt/(m_lt+m_gt))*w_p; //A.2d		
		w_out[2] = (1-p.r_wc)*w_lp; //A.2e
		w_out[3] = p.r_wc*w_lp; //A.2d
		w_lr = p.rho_Lw*p.Q_max*(1-(1-p.C)*(p_bh/p.p_r)-p.C*(p_bh/p.p_r)^2); //A.2g
		w_gr = r_glr*w_lr; //A.2h
		w_gl_max = p.C_gl*sqrt(rho_gl*(p_in - p_ta)); //A.2i
		
		der(m_ga) = w_in-w_gi; //A.1a
		der(m_gt) = w_gr+w_gi-w_out[1]; //A.1b
		der(m_lt) = w_lr-w_lp;  //A.1c
end well;

record pipelineParameters
	import SI = Modelica.SIunits;
	parameter SI.Length L_p = 13000;          // [m]// horizontal pipeline liength
	parameter SI.Radius r_p = 0.10;           // [m] // pipeline radius
	parameter SI.Diameter D_p = 2*r_p;        // [m]	// pipeline diameter	
	parameter SI.Area A_p = Modelica.Constants.pi*r_p^2;     // [m^2]	// pipeline cross sectional area	   
	parameter SI.Volume V_p = L_p*A_p;    // [m^3]	// horizontal pipeline volume	   
	parameter SI.Temperature T_p = 330;            // [K]	// pipeline temperatur	   	
	parameter SI.Temperature T_r = 330;            // [K]	// riser temperature	   	
	parameter SI.DynamicViscosity mu_w = 8.94e-4;       // [Pa*s]	// water dynamic viscosity	   	
	parameter SI.DynamicViscosity mu_o = 1e-4;          // [Pa*s]	// oil dynamic viscosity	   	
	parameter SI.Length L_r = 600;            // [m]	// Riser length	   	
	parameter SI.Radius r_r = 0.10;           // [m]	// Riser radius	   	
	parameter SI.Diameter D_r = 2*r_r;        // [m]	// Riser diameter	   	
	parameter SI.Area A_r = Modelica.Constants.pi*r_r^2;     // [m^2] 	// Riser sectional area	   
	parameter SI.Volume V_r = A_r*L_r;    // [m^3]	// Riser volume	   	
	parameter SI.Length epsilon = 2.8e-5;     // [m]	// Pipeline roughness
	parameter Real K_gp = 0.03;         // [-]	// Orifice equation parameters	   
	parameter Real K_lp = 1;            // [-]	
	parameter Real K_gr = 0.0034;	           // [m^2]	// Valve coefficient for gas/liquid flow through the choke valve at top of riser 
	parameter Real K_lr = 0.0014;	           // [m^2]	
	parameter SI.MolarMass M_g=0.0195;       // [kg/mol]	// gas molar    	
	parameter SI.Density rho_w = 1030;       // [kg/m^3]	// water density	   
	parameter SI.Density rho_o = 930;        // [kg/m^3]	// oil density	   
end pipelineParameters;

model pipeline
	extends Interfaces.TPSISOBlock;
	//******************//***********************//
	//               Parameters					 //
	//******************//***********************//
		parameter pipelineParameters p;	
		import SI = Modelica.SIunits;
		import pi = Modelica.Constants.pi;
		import R = Modelica.Constants.R;
		import g = Modelica.Constants.g_n;
	// Iniital Conditions
		parameter SI.Mass m_gp0(nominal=1e3) = 18632.85129716; // mass of gas in the horizontal pipeline
		parameter SI.Mass m_op0(nominal=1e3) = 8663.20501932; // mass of oil in the horizontal pipeline
		parameter SI.Mass m_wp0(nominal=1e3) = 5496.32160965; // mass of water in the horizontal pipeline
		parameter SI.Mass m_gr0(nominal=1e3) =  563.90194572; // mass of gas in the riser
		parameter SI.Mass m_or0(nominal=1e3) = 1991.90216978; // mass of oil in the riser
		parameter SI.Mass m_wr0(nominal=1e3) = 1330; // mass of water in the riser
		
	//******************//***********************//
	//               States & Variables			 //
	//******************//***********************//	
	//States
	SI.Mass m_gp(start=m_gp0, fixed=true, nominal=1e3, min=0); // mass of gas in the horizontal pipeline
	SI.Mass m_op(start=m_op0, fixed=true, nominal=1e3, min=0); // mass of oil in the horizontal pipeline
	SI.Mass m_wp(start=m_wp0, fixed=true, nominal=1e3, min=0); // mass of water in the horizontal pipeline

	SI.Mass m_gr(start=m_gr0, fixed=true, nominal=1e3, min=0); // mass of gas in the riser
	SI.Mass m_or(start=m_or0, fixed=true, nominal=1e3, min=0); // mass of oil in the riser
	SI.Mass m_wr(start=m_wr0, fixed=true, nominal=1e3, min=0); // mass of water in the riser
	
	// Reference for states
	parameter SI.Mass m_gp_ref = 0; // mass of gas in the horizontal pipeline
	parameter SI.Mass m_op_ref = 0; // mass of oil in the horizontal pipeline
	parameter SI.Mass m_wp_ref = 0; // mass of water in the horizontal pipeline

	parameter SI.Mass m_gr_ref = 0; // mass of gas in the riser
	parameter SI.Mass m_or_ref = 0; // mass of oil in the riser
	parameter SI.Mass m_wr_ref = 0; // mass of water in the riser
	//******************//***********************//
	//         	Auxiliary Variables				 //
	//******************//***********************//	
	// Pipeline
	SI.Volume 					V_gp(nominal=1e6, min=0);	
	SI.Density 					rho_gp(nominal=1e2, min=0);
	SI.Density 					rho_lp(nominal=1e3, min=0);
	SI.VolumeFraction 			alpha_lp(nominal=1e-2, min=0);
	SI.Velocity 				U_slp(nominal=0.1, min=0);
	SI.DynamicViscosity 		mu_lp;
	SI.ReynoldsNumber			Re_p;
	SI.CoefficientOfFriction 	f_p(nominal=1e-2);
	SI.PressureDifference 		DP_fp(nominal=1e3);
	
	// Riser
	SI.Volume 					V_gr(nominal=10, min=0);
	SI.Pressure 				P_r(nominal=1e6, min=0);
	SI.Density 					rho_gr(nominal=10, min=0);
	SI.Density 					rho_lr(nominal=1e3, min=0);
	SI.Density 					rho_mr(nominal=1e2, min=0);
	SI.Density 					rho_t(nominal=1e2, min=0);
	SI.Velocity 				U_slr(nominal=0.1, min=0);
	SI.Velocity 				U_sgr(nominal=1, min=0);
	SI.Velocity 				U_mr(nominal=1, min=0);
	SI.DynamicViscosity 		mu_lr(nominal=1e-4, min=0);
	SI.ReynoldsNumber 			Re_r(nominal=1e5, min=0);
	SI.CoefficientOfFriction 	f_r(nominal=1e-2);
	SI.VolumeFraction 			alpha_lr(nominal=1e-1, min=0);
	SI.PressureDifference		DP_fr(nominal=1e3);
	SI.CoefficientOfFriction	alpha_lt(nominal=1e-1);
	
	// Relations
	SI.Pressure 				p_lp(nominal=1e6, min=0);
	SI.Pressure 				p_lr(nominal=1e6, min=0);
	SI.Area 					A_gp(nominal=1e-2, min=0);
	SI.Area 					A_lp(nominal=1e-3, min=0);
	
	SI.MassFlowRate 			w_glp(nominal=10, min=0);
	SI.MassFlowRate 			w_olp(nominal=10, min=0);
	SI.MassFlowRate 			w_wlp(nominal=10, min=0);
	//******************//***********************//
	//               	Equations				 //
	//******************//***********************//	
	equation
		//// Horizontal Pipeline Equations
		V_gp = p.V_p-m_op/p.rho_o-m_wp/p.rho_w; 								// A.7a - Volume of gas in pipeline
		p_in = (m_gp*R*p.T_p)/(V_gp*p.M_g);										// A.7b - Pressure at pipeline inlet (manifold pressure)
		rho_gp   = m_gp/V_gp;													// A.7c - Density of gas in pipeline
		rho_lp = ((p.rho_o*p.rho_w)/(p.rho_w*m_op+p.rho_o*m_wp))*(m_op+m_wp); 	// A.7d - Density of liquid in pipeline
		alpha_lp = (p.rho_w*m_op+p.rho_o*m_wp)/(p.rho_o*p.rho_w*p.V_p); 		// A.7e - Average liquid volume fraction pipeline
		U_slp = (p.rho_w*w_in[2]+p.rho_o*w_in[3])/(p.rho_o*p.rho_w*pi*p.r_p^2); // A.7f - Superficial velocity of liquid in pipeline
		mu_lp = (m_op/(m_wp+m_op))*p.mu_o + (m_wp/(m_wp+m_op))*p.mu_w;			// A.7g - Liquid dynamic viscosity
		Re_p = (2*rho_lp*U_slp*p.r_p)/mu_lp;									// A.7h - Reynolds number
		f_p = (-1.8*log((p.epsilon/(3.7*p.D_p))^1.11+6.9/(Re_p))/log(10))^(-2); // A.7i - Friction coefficient pipeline
		DP_fp = (alpha_lp*p.L_p*rho_lp*f_p*U_slp^2)/(4*p.r_p); 					// A.7j - Pressure drop due to friction in pipeline
		
		// Riser
		V_gr = p.V_r-m_or/p.rho_o-m_wr/p.rho_w; 								// A.8a - Volume of gas in riser
		P_r = (m_gr*R*p.T_r)/(p.M_g*V_gr); 										// A.8b - Pressure at top of riser
		rho_gr = m_gr/V_gr; 													// A.8c - Density of gas at top of riser
		rho_lr = ((p.rho_o*p.rho_w)/(p.rho_w*m_or+p.rho_o*m_wr))*(m_or+m_wr);	// A.8d - Density of liquid in riser
		rho_mr = (m_gr+m_or+m_wr)/p.V_r; 										// A.8e - Average density of mix in riser
		alpha_lr = (p.rho_w*m_or+p.rho_o*m_wr)/(p.rho_w*p.rho_o*p.V_r); 		// A.8f - Average liquid volume fraction in riser
		U_slr = (p.rho_w*w_in[2]+p.rho_o*w_in[3])/(p.rho_w*p.rho_o*p.A_r); 		// A.8g - Average superficial velocity of liquid in riser
		U_sgr = (w_in[1])/(rho_gr*p.A_r);										// A.8h - Average superficial velocity of gas in riser
		U_mr = U_slr+U_sgr; 													// A.8i - Average superficial velocity of fluid in riser
		mu_lr = (m_wr/(m_wr+m_or))*p.mu_w+(m_or/(m_wr+m_or))*p.mu_o; 			// A.8j - Viscosity of liquid in riser
		Re_r = (2*rho_mr*U_mr*p.r_r)/mu_lr; 									// A.8k - Reynolds number of mixed fluid in riser
		f_r = (-1.8*log((p.epsilon/(3.7*p.D_r))^1.11+6.9/Re_r)/log(10))^(-2);  	// A.8l - Friction coefficient riser
		DP_fr = (alpha_lr*f_r*rho_mr*(U_mr^2)*p.L_r)/(4*p.r_r);					// A.8m - Pressure drop due to friction in riser
		
		// Pipeline-Riser Relation
		p_lp = p_in-DP_fp;														// A.9a - Pressure low point pipeline
		p_lr = P_r+DP_fr+rho_mr*g*p.L_r;										// A.9b - Pressure low point riser		
		A_gp = (V_gp/p.V_p)*p.A_p; 												// A.9c - Cross-sectional area of gas at low point
		A_lp = p.A_p-A_gp;														// A.9d - Cross-sectional area of liqud at low point
		alpha_lt = (2*(m_or+m_wr))/(p.V_r*rho_lr)-A_lp/(p.A_p);					// A.9e - Liquid volume fraction at top of riser
		rho_t = alpha_lt*rho_lr+(1-alpha_lt)*rho_gr;							// A9.f - Density of mixed fluid at top of riser

		// Mass flow of gas, oil and water through low point
	    w_glp = p.K_gp*A_gp*sqrt(rho_gp*maxFunc(p_lp-p_lr, 0.01));								//A.9g - 
	    w_olp = p.K_lp*A_lp*sqrt(p.rho_o*maxFunc(p_lp-p_lr, 0.01))*(m_op/(m_gp+m_op+m_wp));		//A.9h - 
	    w_wlp = p.K_lp*A_lp*sqrt(p.rho_w*maxFunc(p_lp-p_lr, 0.01))*(m_wp/(m_gp+m_op+m_wp));		//A.9i - 

	    // Mass flow of mixed fluid at riser outlet
	    w_out[1] = p.K_gr*sqrt(rho_t*maxFunc(P_r-p_out, 0.01))*(m_gr/(m_gr+m_or+m_wr));			//A.9j
	    w_out[2] = p.K_lr*sqrt(rho_t*maxFunc(P_r-p_out, 0.01))*(m_or/(m_gr+m_or+m_wr));			//A.9k
	    w_out[3] = p.K_lr*sqrt(rho_t*maxFunc(P_r-p_out, 0.01))*(m_wr/(m_gr+m_or+m_wr));			//A.9l

		//State Equaitons   
		der(m_gp) = w_in[1]-w_glp;
	    der(m_op) = w_in[2]-w_olp;
	    der(m_wp) = w_in[3]-w_wlp;
	    der(m_gr) = w_glp-w_out[1];
	    der(m_or) = w_olp-w_out[2];
	    der(m_wr) = w_wlp-w_out[3];
	end pipeline;

model separator
	extends Interfaces.TPSISOBlock;
	import SI = Modelica.SIunits;
	import pi = Modelica.Constants.pi;
	import g = Modelica.Constants.g_n;
	//******************//***********************//
	//               Parameters					 //
	//******************//***********************//	
	//Physical parameters
	parameter SI.Length 			L(nominal=1)= 2.62128 ;  			// Seaparator Length
	parameter SI.Radius 			R(nominal=1)=0.73152;				// Separator Radius 
	parameter SI.Volume 			V_sep(nominal=1) = pi*R^2*L;		// Separator Volume
	parameter SI.Density 			rho_w(nominal=1e3)=1030; 			// Mass denisty of water
	parameter SI.Density 			rho_o(nominal=1e3)=930; 			// Mass denisty of oil
	parameter SI.Diameter 			d_m(nominal=1e-4) = 500*10^(-6); 	// Water droplet average diameter
	parameter SI.DynamicViscosity 	mu_w(nominal=1e-4) = 8.94e-4;       // Water dynamic viscosity[Pa*s] 
	parameter SI.Temperature 		T(nominal=1e2) = 298.15;			// Separator Temperatur
	parameter SI.MolarMass 			M_g(nominal=1e-2) = 0.0195;       	// Gas molar mass [kg/mol]
	constant  SI.Pressure 			P_v(nominal=1e7) = 25834655.6; 		// Stimated Vapor pressure

	//User defined
	parameter SI.Height 			h(nominal=1, start=R) = R;			//Level of water/mixture
	parameter SI.Height 			h_o = 1.4*R;						//Level of Oil
	
	// Geometric Calculations
	parameter SI.Angle 				theta_o = 2*acos((h_o-R)/R);		//G
	parameter SI.Angle 				theta_w = 2*acos((h)/(2*R));
	parameter SI.Volume 			V_g = V_sep -V_w -V_o;
	parameter SI.Volume 			V_o = V_sep - V_w - (R^2*L)/2*(theta_o-sin(theta_o));
	parameter SI.Volume 			V_w = (R^2*L)/2*(theta_o-sin(theta_w));
	parameter SI.Angle 				Phi =atan(h/L); // 0.272145891; //

	parameter SI.Mass 				m_w(nominal=1e3)= V_w*(rho_w);
	parameter SI.Angle 				theta = acos((R-h)/R);					//A.11d

	parameter SI.Velocity			v_v = -(rho_o-rho_w)*g*d_m^2/(18*mu_w); //A.11a
	
	//******************//***********************//
	//         	Auxiliary Variables				 //
	//******************//***********************//	
	SI.Velocity 					v_h(min=0);
	SI.Angle 						Phi_1(min=0, max=pi);
	
	SI.Length 						L_1(start = L, min = L);
	SI.Height 						h_1(start = R, max = h, min=0);
	SI.Angle 						theta_1(start=theta, min=0, max=pi);
	SI.Volume 						V_s1(start=1, min=0);
	SI.Volume 						V_s2(start=1, min=0);
	SI.VolumeFraction 				epsilon(start=1, min =0, max =1);
	
	//SI.MassFlowRate 				w_o_wout(start=1, min=0);
	//SI.MassFlowRate 				w_g_wout(start=0);
	SI.MassFlowRate 				w_g_oout(start=0);
	Real 							x(start=0);
	//******************//***********************//
	//               	Equations				 //
	//******************//***********************//	
equation
	p_out = p_in;
	v_h = (L/(m_w/w_out[3]));					//A.11b and .11c
	Phi_1 = atan(v_v/v_h);						//A.11f
	L_1 = h/tan(-maxFunc(-Phi_1,-Phi));			//A.11g		
	h_1 = L*tan(-maxFunc(-Phi_1,-Phi));  		//A.11h
	theta_1 = acos(1-h_1/R);					//A.11i

	V_s1 = R^2*L_1*(theta-0.5*sin(2*theta)-(3*sin(theta)-3*theta*cos(theta)-sin(theta)^3)/(3*(1-cos(theta)))); 				//A.11j
	V_s2 = R^2*L*(theta-0.5*sin(2*theta) -(3*sin(theta_1)-3*theta_1*cos(theta_1)-sin(theta_1)^3)/(3*(1-cos(theta_1) ))); 	//A.11k

	//Water Phase
	epsilon  = V_s2/V_s1;						//A.12a
	//w_o_wout = (1-epsilon)*w_in[1];				//A.12b
	//w_g_wout = (1-epsilon)*w_in[2];				//A.12c
	w_out[3] = w_in[3];							//A.12d
	
	//Oil Phase
	x = p_in/P_v;								//A.13a
	w_g_oout = x*epsilon*w_in[1];				//A.13b
	w_out[2] = epsilon*w_in[2];					//A.13c
	
	//Gas Phase 
	w_out[1] = (1-x)*epsilon*w_in[1];			//A.14a
end separator;
model compressor "Compressor M odel"
	extends Interfaces.SPSISOBlock;
	import SI = Modelica.SIunits;
	import R = Modelica.Constants.R;

	Interfaces.flowPortIn w_sl;
	
	parameter SI.Temperature T = 325; // [K]	
	parameter SI.MolarMass M_g(nominal=1e-2) = 0.0195;       // [kg/mol]
	   
	parameter SI.Pressure p_map = 30e5; //Pressure that the map was made
	parameter SI.Temperature T_map(nominal=1e5) = 318.15; //Tempearature that the map was made

	Real f_c = 0.1;
	parameter Real au = -0.2354, bu =1.3593, cu = 3.2619;  //Quadratic parameters of upper bound
	parameter Real ad = -0.1484, bd =0.6094, cd = 2.0621;  //Quadratic parameters of 'down' bound
	parameter Real al = 1.1284, bl = -3.2777, cl = 4.2217;  //Quadratic parameters of left bound 
	parameter Real ar = 0.2774, br = -1.9135, cr = 5.1650;  //Quadratic parameters of right bound


	SI.Pressure nu_u(nominal=1e5); 	// Up Slack Variables
	SI.Pressure nu_d(nominal=1e5);	// Up Slack Variables
	SI.Pressure nu_l(nominal=1e5);	// Up Slack Variables
	SI.Pressure nu_r(nominal=1e5);	// UpSlack Variables

	Real  p_ratio;// p_down/p_up converted to the test conditions
	SI.VolumeFlowRate q_vol; //Volume converted to test conditions
equation	
	w_out = w_in; //A.14a
	q_vol = (w_in+w_sl)*R*T/(p_in*M_g)*((T_map/T)*(p_in/p_map)); //A.14d and A.14e with A.14b
	
	p_ratio = (p_out/p_in)*(T/T_map); // A.15a and A.15b
	
	nu_u = p_ratio-(au*q_vol^2+bu*q_vol+cu); //A.16a and .17a
	nu_d = (ad*q_vol^2+bd*q_vol+cd)-p_ratio; //A.16a and .17b
	nu_l = f_c*(-bl+sqrt(bl^2-4*al*(cl-p_ratio)))/(2*al)-q_vol; //A.16b and .17c
	nu_r = q_vol-f_c*(-br+sqrt(br^2-4*ar*(cr-p_ratio)))/(2*ar); //A.16c and .17d
end compressor;
end networkComponents;

package networkComponentsApproximation
	extends networkComponents;

	function maxFuncApproximation
		input Real x1, x2;
		output Real y;
	algorithm
		y := ((x1-x2)^2+Modelica.Constants.eps^2)^0.5/2+(x1+x2)/2;
	end maxFuncApproximation;

	redeclare function maxFunc =  maxFuncApproximation;
end networkComponentsApproximation;



model network
	import SI = Modelica.SIunits;
	replaceable package components = networkComponents;
	
	parameter Integer n_wells =8;
	parameter Integer n_pipelines =2;
	parameter Integer n_separators = 2;
	parameter Integer n_compressors = 2;
	parameter SI.Pressure p_export(nominal=1e5) = 200e5;
	//SI.Pressure p_s[n_separators](each nominal=1e5);

	
	//Network Input
	Interfaces.flowPortIn w_gl[n_wells](each min=1, each max=3, each free=true);	
	Interfaces.flowPortIn w_fl[n_compressors](each min=0, each start=0, each free=true);
	Interfaces.flowPortIn w_sl[n_compressors](each min=0);
	Interfaces.pressurePortIn p_s[n_separators](each nominal=1e5, each min=50e5, each max=60e5, each start=50e5);
	
	//Network Output
	SI.MassFlowRate w_out[3](each min=0, each start =0);
	
	//Network Internal Variables
	SI.MassFlowRate w_gm(min=0, start = 0);
	
	model compressorModel = components.compressor;
	//System Components
	components.well w[n_wells];
	components.pipeline p[n_pipelines];
	components.separator s[n_separators];
	compressorModel c[n_compressors];
	equation 
		//p_s[1] = 50e5+components.heaviside(time-15000)*10e5;
		//p_s[2] = 60e5+components.heaviside(time-15000)*5e5;

		// Well -> Manifold
		for i in 1:4 loop 
			connect(w[i].w_in, w_gl[i]);
			connect(w[i].p_out, p[1].p_in);
			w[i].p_in = p_export; 
		end for;
		
		for i in 5:8 loop 
			connect(w[i].w_in, w_gl[i]);
			connect(w[i].p_out, p[2].p_in);
			w[i].p_in = p_export;
		end for;
		

		p[1].w_in[:] = sum(w[i].w_out[:] for i in 1:4);
		p[2].w_in[:] = sum(w[i].w_out[:] for i in 5:8);
		
		// pipeline separator
		for i in 1:n_separators loop
			s[i].w_in[:] = p[i].w_out[:];
			connect(p[i].p_out, s[i].p_in);
			s[i].p_out = p_s[i];
		end for;
		// Separator -> Compressor
		connect(c[:].w_sl, w_sl[:]);
		
		connect(s[1].p_out, c [1].p_in);
		c[1].w_in = sum(s[i].w_out[1] for i in 1:1)- w_fl[1];		
		c[1].p_out = p_export;
		
		connect(s[2].p_out, c[2].p_in);
		c[2].w_in = sum(s[i].w_out[1] for i in 2:2) - w_fl[2];		
		c[2].p_out = p_export;
		
		w_gm = sum(w_gl[i] for i in 1:n_wells);
		//Network Outputs
		w_out[1] = sum(c[i].w_out for i in 1:n_compressors)-w_gm;
		w_out[2] = sum(s[i].w_out[2] for i in 1:n_separators);
		w_out[3] = sum(s[i].w_out[3] for i in 1:n_separators);


end network;
							
model networkApproximation
	extends network;
	redeclare package components = networkComponentsApproximation;			
end networkApproximation;
													
end networkModels;						
															
																
																	
																		
																			
																				
																					
																						