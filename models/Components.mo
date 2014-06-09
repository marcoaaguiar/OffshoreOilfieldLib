within Network;

package Components 
	import SI = Modelica.SIunits;
	import BC = BoundaryConstants;
	replaceable function maxFunc = Functions.maxFunc;

	constant Real SqrtEpsilon = 1e-1; 
	constant SI.Mass minMass = 1;
	constant SI.Volume minVolume = 1e-1;
	constant SI.Pressure minPressure = 1e5;
	constant SI.Density minDensity = 1;
	constant SI.Velocity minVelocity = 1e-3;
	constant SI.Angle angleEpsilon = 1e-3;
	
	model ProductionManifold
	  extends Interfaces.ProductionManifoldInterface;
	  //Real rounting_table[:];
	  equation  
	    for i in 1:n_inputs loop
	      connect(outputConnector.pressure, inputConnector[i].pressure);
	    end for;
	    outputConnector.massFlow = sum(inputConnector[i].massFlow for i in 1:n_inputs);
	end ProductionManifold;
	
	model Flare
	    extends Interfaces.FlareInterface;
	    input SI.MassFlowRate w_fl(min = 0); 
	    equation
	      outputConnector.pressure = inputConnector.pressure;
	      outputConnector.massFlow[1] = inputConnector.massFlow[1] - w_fl;
	      outputConnector.massFlow[2] = inputConnector.massFlow[2];
	      outputConnector.massFlow[3] = inputConnector.massFlow[3];
	end Flare;
	
	model Valve
		parameter Real C = 0.0014;
	
		SI.Density rho;
		input SI.Pressure p_upstream, p_downstream;
		output SI.MassFlowRate w;
		
		equation
			w = C*sqrt(maxFunc(0.1, rho*(p_upstream-p_downstream)));
	end Valve;
	
	model Well "Well Model"
		//******************//***********************//
		//               Parameters					 //
		//******************//***********************//
			extends Interfaces.WellInterface; 
		  	  
			import SI = Modelica.SIunits;
			import R = Modelica.Constants.R;
			import g = Modelica.Constants.g_n;
			
			parameter ComponentsProperties.MediaProperty mediaP;
			parameter ComponentsProperties.WellProperty wellP;
			parameter ComponentsProperties.ReservoirProperty reservoirP;
			parameter ComponentsProperties.WellValvesParameters valveP;
			
			parameter SI.Density rho_Lw = (mediaP.rho_w * reservoirP.r_wc + mediaP.rho_o * (1-reservoirP.r_wc));  	//A.4b - liquid density
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
			SI.Mass m_ga(start= m_ga0, fixed=true, nominal=1e3, min = minMass); // mass of gas in the well annular
			SI.Mass m_gt(start= m_gt0, fixed=true, nominal=1e3, min = minMass); // mass of gas in the well tubing
			SI.Mass m_lt(start= m_lt0, fixed=true, nominal=1e3, min = minMass); // mass of liquid in the well tubing
			
		//Variables 
			SI.Pressure p_ai(nominal=1e6, min = minPressure);
			SI.Pressure p_p	(nominal=1e6, min = minPressure);
			SI.Pressure p_ti(nominal=1e6, min = minPressure);
			SI.Pressure p_bh(nominal=1e6, min = minPressure);
			SI.Pressure p_ta(nominal=1e6, min = minPressure);
			
			SI.Density rho_gi(nominal=1e2, min = minDensity);
			SI.Density rho_p (nominal=1e2, min = minDensity);
			SI.Density rho_gl(nominal=1e2, min = minDensity);
	
			SI.MassFlowRate w_gi(min=0, start=0);		
			SI.MassFlowRate w_p(min=0, start=0);
			SI.MassFlowRate w_lp(min=0, start=0);
			SI.MassFlowRate w_lr(min=0, start=0);
			SI.MassFlowRate w_gr(min=0, start=0);
			SI.MassFlowRate w_gl_max(min=0, start=0);
			
			SI.MassFlowRate w_in, w_out[3];
			SI.Pressure p_in, p_out;
			
			
			Valve GasLiftChoke(C = valveP.C_gl, rho = rho_gl, p_upstream= p_in, p_downstream = p_ta, w = w_gl_max);
			Valve InjectionValve(C = valveP.C_iv, rho = rho_gi, p_upstream= p_ai, p_downstream = p_ti, w = w_gi);
			Valve ProductionChoke(C = valveP.C_pc, rho = rho_p, p_upstream = p_p, p_downstream = p_out, w = w_p);
		//******************//***********************//
		//               Equations					 //
		//******************//***********************//			
		equation 
			w_in = inputConnector.massFlow[1];
			w_out = outputConnector.massFlow;
			p_in = inputConnector.pressure;
			p_out = outputConnector.pressure;
			
			p_ai = ((R*wellP.T_a)/(wellP.V_a*mediaP.M_g) + g/(2*wellP.A_a))*m_ga;  //A.3a
			p_p = (R*wellP.T_t*m_gt)/(mediaP.M_g*wellP.V_t-mediaP.M_g*(1/rho_Lw)*m_lt)-(g*(m_gt))/(2*wellP.A_t); //A.3b
			p_ti = p_p+(g*(m_lt+m_gt))/wellP.A_t;  //A.3c
			p_bh = ((1+reservoirP.r_glr+(reservoirP.r_glr*mediaP.M_g*g*wellP.L_w)/(2*R*wellP.T_t))*p_ti + rho_Lw*g*wellP.L_w)/(1+reservoirP.r_glr-(reservoirP.r_glr*mediaP.M_g*g*wellP.L_w)/(2*R*wellP.T_t)); //A.3d
			p_ta = ((R*wellP.T_a)/(wellP.V_a*mediaP.M_g)-g/(2*wellP.A_a))*m_ga;		//A.3e - Pressure in annulus at gas-lift choke valve
			
			rho_gi = (mediaP.M_g/(R*wellP.T_a)*p_ai); //A.4a
			rho_p = (rho_Lw*mediaP.M_g*p_p*(m_lt+m_gt))/(rho_Lw*R*wellP.T_t*m_lt+mediaP.M_g*p_p*m_gt); //A.4b
			rho_gl = mediaP.M_g*p_in/(R*wellP.T_a); //A.4d
			
			//w_gi = valveP.C_iv*sqrt(maxFunc(Components.SqrtEpsilon, rho_gi*(p_ai-p_ti))); //A.2a
			//w_p = valveP.C_pc*sqrt(maxFunc(Components.SqrtEpsilon, rho_p*(p_p-p_out))); //A.2b
			w_lp = (m_lt/(m_lt+m_gt))*w_p; //A.2d		
			w_lr = rho_Lw*reservoirP.Q_max*(1-(1-reservoirP.C)*(p_bh/reservoirP.p_r)-reservoirP.C*(p_bh/reservoirP.p_r)^2); //A.2g
			w_gr = reservoirP.r_glr*w_lr; //A.2h
			//w_gl_max = valveP.C_gl*sqrt((rho_gl*(p_in- p_ta))); //A.2i

			w_out[1]= (m_gt/(m_lt+m_gt))*w_p; //A.2c
			w_out[2] = (1-reservoirP.r_wc)*w_lp; //A.2e
			w_out[3] = reservoirP.r_wc*w_lp; //A.2d
			
			der(m_ga) = w_in-w_gi; //A.1a
			der(m_gt) = w_gr+w_gi-w_out[1]; //A.1b
			der(m_lt) = w_lr-w_lp;  //A.1c
	end Well;
	
	model Pipeline
		extends Interfaces.PipelineRiserInterface;
		//******************//***********************//
		//               Parameters					 //
		//******************//***********************//
			parameter ComponentsProperties.pipelineParameters p;	
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
		SI.Mass m_gp(start=m_gp0, fixed=true, nominal=1e3, min=minMass); // mass of gas in the horizontal pipeline
		SI.Mass m_op(start=m_op0, fixed=true, nominal=1e3, min=minMass); // mass of oil in the horizontal pipeline
		SI.Mass m_wp(start=m_wp0, fixed=true, nominal=1e3, min=minMass); // mass of water in the horizontal pipeline
	
		SI.Mass m_gr(start=m_gr0, fixed=true, nominal=1e3, min=minMass); // mass of gas in the riser
		SI.Mass m_or(start=m_or0, fixed=true, nominal=1e3, min=minMass); // mass of oil in the riser
		SI.Mass m_wr(start=m_wr0, fixed=true, nominal=1e3, min=minMass); // mass of water in the riser
		
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
		SI.Volume 					V_gp	(nominal=1e6, min=minVolume);	
		SI.Density 					rho_gp	(nominal=1e2, min=minDensity);
		SI.Density 					rho_lp	(nominal=1e3, min=minDensity);
		SI.VolumeFraction 			alpha_lp(nominal=1e-2, min=0);
		SI.Velocity 				U_slp	(nominal=0.1, min=0);
		SI.DynamicViscosity 		mu_lp	(nominal=1e-2);
		SI.ReynoldsNumber			Re_p	(nominal=1e5);
		SI.CoefficientOfFriction 	f_p		(nominal=1e-2);
		SI.PressureDifference 		DP_fp	(nominal=1e3);
		
		// Riser
		SI.Volume 					V_gr	(nominal=10, min=minVolume);
		SI.Pressure 				P_r		(nominal=1e6, min=minPressure);
		SI.Density 					rho_gr	(nominal=10, min=minDensity);
		SI.Density 					rho_lr	(nominal=1e3, min=minDensity);
		SI.Density 					rho_mr	(nominal=1e2, min=minDensity);
		SI.Density 					rho_t	(nominal=1e2, min=minDensity);
		SI.Velocity 				U_slr	(nominal=0.1, min=0);
		SI.Velocity 				U_sgr	(nominal=1, min=0);
		SI.Velocity 				U_mr	(nominal=1, min=0);
		SI.DynamicViscosity 		mu_lr	(nominal=1e-4, min=0);
		SI.ReynoldsNumber 			Re_r	(nominal=1e5, min=0);
		SI.CoefficientOfFriction 	f_r		(nominal=1e-2);
		SI.VolumeFraction 			alpha_lr(nominal=1e-1, min=0);
		SI.PressureDifference		DP_fr	(nominal=1e3);
		SI.CoefficientOfFriction	alpha_lt(nominal=1e-1);
		
		// Relations
		SI.Pressure 				p_lp(nominal=1e6, min=minPressure);
		SI.Pressure 				p_lr(nominal=1e6, min=minPressure);
		SI.Area 					A_gp(nominal=1e-2, min=0);
		SI.Area 					A_lp(nominal=1e-3, min=0);
		
		SI.MassFlowRate 			w_glp(nominal=10, min=0);
		SI.MassFlowRate 			w_olp(nominal=10, min=0);
		SI.MassFlowRate 			w_wlp(nominal=10, min=0);
		SI.MassFlowRate w_in[3], w_out[3];
		SI.Pressure p_in, p_out;
		//******************//***********************//
		//               	Equations				 //
		//******************//***********************//	
		equation
			w_in = inputConnector.massFlow;
			w_out = outputConnector.massFlow;
			p_in = inputConnector.pressure;
			p_out = outputConnector.pressure;
			
			//// Horizontal Pipeline Equations
			V_gp = p.V_p-m_op/p.rho_o - m_wp/p.rho_w; 								// A.7a - Volume of gas in pipeline
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
			
/* 			w_out[1] = (1-alpha_lt)*p.K_gr*sqrt(rho_t*maxFunc(P_r-p_out, 0.01));			//A.9j
			w_out[2] = alpha_lt*p.K_lr*sqrt(rho_t*maxFunc(P_r-p_out, 0.01))*(m_or/(m_or+m_wr));			//A.9k
			w_out[3] = alpha_lt*p.K_lr*sqrt(rho_t*maxFunc(P_r-p_out, 0.01))*(m_wr/(m_or+m_wr));			//A.9l */

			//State Equaitons   
			der(m_gp) = w_in[1]-w_glp;
			der(m_op) = w_in[2]-w_olp;
			der(m_wp) = w_in[3]-w_wlp;
			der(m_gr) = w_glp-w_out[1];
			der(m_or) = w_olp-w_out[2];
			der(m_wr) = w_wlp-w_out[3];
		end Pipeline;
	
	model Separator
		extends Interfaces.SeparatorInterface;
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
		parameter SI.Angle 				Phi = 0.272145891;//atan(h/L);
	
		parameter SI.Mass 				m_w(nominal=1e3)= V_w*(rho_w);
		parameter SI.Angle 				theta = acos((R-h)/R);					//A.11d
	
		parameter SI.Velocity			v_v = -(rho_o-rho_w)*g*d_m^2/(18*mu_w); //A.11a
		
		//******************//***********************//
		//         	Auxiliary Variables				 //
		//******************//***********************//	
		SI.Velocity 					v_h(min=0);
		SI.Angle 						Phi_1(min=0+angleEpsilon, max=pi-angleEpsilon);
		
		SI.Length 						L_1(start = L, min = L);
		SI.Height 						h_1(start = R, max = h, min=0);
		SI.Angle 						theta_1(start=theta, min=0+angleEpsilon, max=pi-angleEpsilon);
		SI.Volume 						V_s1(start=1, min=0);
		SI.Volume 						V_s2(start=1, min=0);
		SI.VolumeFraction 				epsilon(start=1, min =0, max =1);

		Real 							x(start=0, min=0);
		
		SI.MassFlowRate					w_in[3](each min = 0);
		SI.Pressure						p_in; 
		//******************//***********************//
		//               	Equations				 //
		//******************//***********************//	
	equation
		w_in = inputConnector.massFlow;
		p_in = inputConnector.pressure;
		
		outputConnector[1].pressure = p_in;
		outputConnector[2].pressure = p_in;
		outputConnector[3].pressure = p_in;
		
		v_h = (L/(m_w/outputConnector[3].massFlow[3]));					//A.11b and .11c
		//v_h = L*(outputConnector[3].massFlow[3]/m_w);					//A.11b and .11c
		Phi_1 = atan(v_v/v_h);						//A.11f
		//tan(Phi_1) = (v_v/v_h);						//A.11f
		L_1 = h/tan(-maxFunc(-Phi_1,-Phi));			//A.11g		
		h_1 = L*tan(-maxFunc(-Phi_1,-Phi));  		//A.11h
		theta_1 = acos(1-h_1/R);					//A.11i
	
		V_s1 = R^2*L_1*(theta-0.5*sin(2*theta)-(3*sin(theta)-3*theta*cos(theta)-sin(theta)^3)/(3*(1-cos(theta)))); 				//A.11j
		V_s2 = R^2*L*(theta-0.5*sin(2*theta) -(3*sin(theta_1)-3*theta_1*cos(theta_1)-sin(theta_1)^3)/(3*(1-cos(theta_1) ))); 	//A.11k
	
		//Gas Phase 
		outputConnector[1].massFlow[1] = (1-x)*epsilon*w_in[1];			//A.14a
		outputConnector[1].massFlow[2] = 0; // no oil through the gas outlet
		outputConnector[1].massFlow[3] = 0; // no water through the gas outlet
		

		
		//Oil Phase
		x = p_in/P_v;								//A.13a
		outputConnector[2].massFlow[1] = x*epsilon*w_in[1];				//A.13b
		outputConnector[2].massFlow[2] = epsilon*w_in[2];					//A.13c
		outputConnector[2].massFlow[3] = 0; // no water through the gas outlet

		//Water Phase
		epsilon  = V_s2/V_s1;						//A.12a
		outputConnector[3].massFlow[1] = (1-epsilon)*w_in[1];				//A.12b
		outputConnector[3].massFlow[2] = (1-epsilon)*w_in[2];				//A.12c
		outputConnector[3].massFlow[3] = w_in[3];		//A.12d		

	end Separator;
	
	model Compressor "Compressor Model"
		extends Interfaces.CompressorInterface;
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
	
		Real p_ratio;// p_down/p_up converted to the test conditions
		SI.VolumeFlowRate q_vol; //Volume converted to test conditions
		
		SI.MassFlowRate w_in, w_out;
		SI.Pressure p_in, p_out;
	equation	
		w_in = inputConnector.massFlow[1]; //A.14a
		w_out = w_in;
		outputConnector.massFlow[1] = w_out;
		outputConnector.massFlow[2] = 0;
		outputConnector.massFlow[3] = 0;
		
		p_in = inputConnector.pressure;
		p_out = outputConnector.pressure;
		
		q_vol = (outputConnector.massFlow[1]+w_sl)*R*T/(p_in*M_g)*((T_map/T)*(p_in/p_map)); //A.14d and A.14e with A.14b
		
		p_ratio = (p_out/p_in)*(T/T_map); // A.15a and A.15b
		
		nu_u = p_ratio-(au*q_vol^2+bu*q_vol+cu); //A.16a and .17a
		nu_d = (ad*q_vol^2+bd*q_vol+cd)-p_ratio; //A.16a and .17b
		nu_l = f_c*(-bl+sqrt(bl^2-4*al*(cl-p_ratio)))/(2*al)-q_vol; //A.16b and .17c
		nu_r = q_vol-f_c*(-br+sqrt(br^2-4*ar*(cr-p_ratio)))/(2*ar); //A.16c and .17d
	end Compressor;
end Components;