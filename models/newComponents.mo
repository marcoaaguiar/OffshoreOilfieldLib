within Network;

package newComponents 
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
		
	model Valve
		parameter Real C = 0.0014;
	
		SI.Density rho(start = 1000);
		input SI.Pressure p_upstream(start = 100e5), p_downstream(start = 50e5);
		input Real u(min = 0, max = 1);
		output SI.MassFlowRate w;
		
		equation
			w = C*sqrt(maxFunc(0.1, rho*(p_upstream-p_downstream)))*u;
			//w = C*(sqrt(maxFunc(0, rho*(p_upstream-p_downstream))+0.1) - sqrt(0.1))*u;
	end Valve;
	
	model Annulus "New Annulus II"
		extends Interfaces.WellInterface; 
		  
		import SI = Modelica.SIunits;
		import R = Modelica.Constants.R;
		import g = Modelica.Constants.g_n;
		
		parameter ComponentsProperties.MediaProperty mediaP;
		parameter ComponentsProperties.WellProperty wellP;
		parameter ComponentsProperties.WellValvesParameters valveP;

		parameter SI.Mass m_ga0(nominal=1e3) = 3629.07517367; 

		//States
		SI.Mass m_ga(start= m_ga0, fixed=true, nominal=1e3, min = 1); // mass of gas in the well annular
		
		//Variables 
		SI.Pressure p_ai(nominal=1e6, min = 1e5);
		SI.Pressure p_ta(nominal=1e6, min = 1e5);
		SI.Density rho_gi(nominal=1e2, min = 1);

		SI.MassFlowRate w_gi(min=0, start=1);		

		// Valve InjectionValve(C = valveP.C_iv, rho = rho_gi, p_upstream= p_ai, p_downstream = p_out, w = w_gi, u = 1);
		
		SI.MassFlowRate w_in; 
		SI.MassFlowRate w_out[3] = outputConnector.massFlow;
		SI.Pressure p_in = inputConnector.pressure;
		SI.Pressure	p_out = outputConnector.pressure;
		
		SI.Velocity U_ga;
		SI.ReynoldsNumber Re;
		SI.CoefficientOfFriction f_D;
		Real C,C1,C2;
		equation
			inputConnector.massFlow[1] = w_in;

			w_out[1] = w_gi;
			w_out[2] = 0;
			w_out[3] = 0;
		
			U_ga = w_in*wellP.L_a/m_ga;
			Re = w_in*wellP.D_a/(mediaP.mu_g*wellP.A_a);
			f_D = (-1.8*log10((wellP.epsilon_a/(3.7*wellP.D_a))^(1.11) + 6.9/Re))^(-2);
			C1 = mediaP.M_g/(R*wellP.T_a)*(-g);
			C2 = mediaP.M_g/(R*wellP.T_a)*(+f_D*U_ga^2/(2*wellP.D_a));
			C = C1+C2;
			
			p_ai = (g*m_ga/wellP.A_a - f_D*wellP.L_a*m_ga*U_ga^2/(2*wellP.V_a*wellP.D_a))/(1-exp(C*wellP.L_a));
			p_ta = p_ai - (g*m_ga)/(wellP.A_a) + f_D*m_ga*U_ga^2/(2*wellP.D_a*wellP.A_a);
			rho_gi = (mediaP.M_g/(R*wellP.T_a)*p_ai);
			w_gi = valveP.C_iv*sqrt(maxFunc(Components.SqrtEpsilon, rho_gi*(p_ai-p_out)));
			
			der(m_ga) = inputConnector.massFlow[1]-w_gi;			
	end Annulus;
	
	model Tubing "New Tubing I"
		extends Interfaces.WellInterface; 
		  
		import SI = Modelica.SIunits;
		import R = Modelica.Constants.R;
		import g = Modelica.Constants.g_n;
		
		parameter ComponentsProperties.MediaProperty mediaP;
		parameter ComponentsProperties.WellProperty wellP;
		parameter ComponentsProperties.WellValvesParameters valveP;
		parameter ComponentsProperties.ReservoirProperty reservoirP;
		
		parameter SI.Density rho_l = (mediaP.rho_w * reservoirP.r_wc + mediaP.rho_o * (1-reservoirP.r_wc));
		
		//Initialization	
		parameter SI.Mass m_gt0(nominal=1e3) = 1389.20150433;
		parameter SI.Mass m_lt0(nominal=1e3) = 3352.26319015;	

		//******************//***********************//
		//               States & Variables			 //
		//******************//***********************//	
		
		SI.Mass m_gt(start= m_gt0, fixed=true, nominal=1e3, min = minMass); 
		SI.Mass m_lt(start= m_lt0, fixed=true, nominal=1e3, min = minMass); 
		Real x (min = 0, max = 1, start = 0.3); 
		Real alpha_wh (min = 0, max = 1, start = 0.5);
		Real alpha_ip (min = 0, max = 1, start = 0.5);
		
		SI.Pressure p_p		(nominal=1e6, min = minPressure, start = 100e5);
		SI.Pressure p_ip	(nominal=1e6, min = minPressure);
		SI.Pressure p_bh	(nominal=1e6, min = minPressure, start = 150e5);

		SI.Density rho_p 	(nominal=1e2, min = minDensity);

		SI.MassFlowRate w_gi(min=0, start=0);		
		SI.MassFlowRate w_p	(min=0, start=0);
		SI.MassFlowRate w_lp(min=0, start=0);
		SI.MassFlowRate w_gp(min=0, start=0);
		SI.MassFlowRate w_lr(min=0, start=0);
		SI.MassFlowRate w_gr(min=0, start=0);
		
		SI.ReynoldsNumber Re;
		SI.DynamicViscosity mu_m;
		
		SI.MassFlowRate w_in = inputConnector.massFlow[1]; 
		SI.MassFlowRate w_out[3] = outputConnector.massFlow;
		SI.Pressure p_in = inputConnector.pressure;
		SI.Pressure	p_out = outputConnector.pressure;
		
		Valve ProductionChoke(C = valveP.C_pc, rho = rho_p, p_upstream = p_p, p_downstream = p_out, u = 1);
	equation
		w_gi = w_in;
		p_in = p_ip;
		
		mu_m = x*(mediaP.mu_g) + (1-x)*((1-reservoirP.r_wc)*mediaP.mu_o + reservoirP.r_wc*mediaP.mu_w);
		
		//x = (w_gr +w_gi)/(w_gr + w_gi + w_lr);
		//x = w_gp/(w_gp + w_lp);
		der(x)  = (-x + (0.5*(w_gi + w_gr)/(w_gi + w_gr +w_lr) + 0.5*w_gp/w_p))/1000;
		alpha_wh = x*rho_l/((1-x)*p_p*mediaP.M_g/(R*wellP.T_t) + x*rho_l);
		alpha_ip = x*rho_l/((1-x)*p_ip*mediaP.M_g/(R*wellP.T_t) + x*rho_l);
		
		//p_p = p_ip*exp((-wellP.L_t*g*mediaP.M_g*rho_l + (1-x)*mediaP.M_g*g*(m_gt+m_lt)/wellP.A_t)/(x*rho_l*R*wellP.T_t));
		//p_ip = p_p + g*(m_gt+m_lt)/wellP.A_t;
		
		Re = w_p*wellP.D_t/(mu_m*wellP.A_t);
		
		p_ip = g*(m_gt+m_lt)/wellP.A_t/(1-exp((-wellP.L_t*g*mediaP.M_g*rho_l + (1-x)*mediaP.M_g*g*(m_gt+m_lt)/wellP.A_t)/(x*rho_l*R*wellP.T_t)));
		p_p = p_ip - g*(m_gt+m_lt)/wellP.A_t;
		
		p_bh = p_ip;
		
		rho_p = p_p*mediaP.M_g*rho_l/((1-x)*p_p*mediaP.M_g + x*rho_l*R*wellP.T_t);
		
		w_p = ProductionChoke.w;
		w_lr = maxFunc(0.1,rho_l*reservoirP.Q_max*(1-(1-reservoirP.C)*(p_bh/reservoirP.p_r)-reservoirP.C*(p_bh/reservoirP.p_r)^2)); //A.2g
		w_gr = reservoirP.r_glr*w_lr; //A.2h		
		
		//w_lp = (1-alpha_wh)*w_p; //(m_lt/(m_lt+m_gt))*w_p; //A.2d		
		w_lp = (m_lt/(m_lt+m_gt))*w_p; //A.2d		
		//w_gp = (alpha_wh)*w_p; //(m_lt/(m_lt+m_gt))*w_p; //A.2d		
		w_gp = (m_gt/(m_lt+m_gt))*w_p; //A.2d		
		
		w_out[1] = w_gp;//(m_gt/(m_lt+m_gt))*w_p; //A.2c
		w_out[2] = (1-reservoirP.r_wc)*w_lp; //A.2e
		w_out[3] = reservoirP.r_wc*w_lp; //A.2d
			
		der(m_gt) = w_gr + w_gi - w_gp; //A.1b
		der(m_lt) = w_lr - w_lp;  //A.1c
		
	end Tubing;
	
	model Tubing2 "New Tubing II"
		extends Interfaces.WellInterface; 
		  
		import SI = Modelica.SIunits;
		import R = Modelica.Constants.R;
		import g = Modelica.Constants.g_n;
		
		parameter ComponentsProperties.MediaProperty mediaP;
		parameter ComponentsProperties.WellProperty wellP;
		parameter ComponentsProperties.WellValvesParameters valveP;
		parameter ComponentsProperties.ReservoirProperty reservoirP;
		
		parameter SI.Density rho_l = (mediaP.rho_w * reservoirP.r_wc + mediaP.rho_o * (1-reservoirP.r_wc));
		
		//Initialization	
		parameter SI.Mass m_gt0(nominal=1e3) = 1389.20150433;
		parameter SI.Mass m_lt0(nominal=1e3) = 3352.26319015;	

		//******************//***********************//
		//               States & Variables			 //
		//******************//***********************//	
		
		SI.Mass m_gt(start= m_gt0, fixed=true, nominal=1e3, min = minMass); 
		SI.Mass m_lt(start= m_lt0, fixed=true, nominal=1e3, min = minMass); 
		Real x (min = 0, max = 1, start = 0.3); 
		Real alpha_wh (min = 0, max = 1, start = 0.5);
		Real alpha_ip (min = 0, max = 1, start = 0.5);
		
		SI.Pressure p_p		(nominal=1e6, min = minPressure, start = 100e5);
		SI.Pressure p_ip	(nominal=1e6, min = minPressure);
		SI.Pressure p_bh	(nominal=1e6, min = minPressure, start = 150e5);

		SI.Density rho_p 	(nominal=1e2, min = minDensity);

		SI.MassFlowRate w_gi(min=0, start=0);		
		SI.MassFlowRate w_p	(min=0, start=0);
		SI.MassFlowRate w_lp(min=0, start=0);
		SI.MassFlowRate w_gp(min=0, start=0);
		SI.MassFlowRate w_lr(min=0, start=0);
		SI.MassFlowRate w_gr(min=0, start=0);
		

		
		SI.MassFlowRate w_in = inputConnector.massFlow[1]; 
		SI.MassFlowRate w_out[3] = outputConnector.massFlow;
		SI.Pressure p_in = inputConnector.pressure;
		SI.Pressure	p_out = outputConnector.pressure;
		
		Valve ProductionChoke(C = valveP.C_pc, rho = rho_p, p_upstream = p_p, p_downstream = p_out, u = 1);
		
		Real C(start = -180);
		SI.Pressure Delta_mt, Delta_ut(start =3e5);
		SI.Velocity U_m(start = 2);
		SI.ReynoldsNumber Re(start = 1e7);
		SI.CoefficientOfFriction f_D(start = 0.014);
		SI.DynamicViscosity mu_m(start = 1e4);
	equation
		w_gi = w_in;
		p_in = p_ip;
		// Mass Pressure Drop 

		
		//x = (w_gr +w_gi)/(w_gr + w_gi + w_lr);
		//x = w_gp/(w_gp + w_lp);
		der(x)  = (-x + (0.5*(w_gi + w_gr)/(w_gi + w_gr +w_lr) + 0.5*w_gp/w_p))/1000;
		alpha_wh = x*rho_l/((1-x)*p_p*mediaP.M_g/(R*wellP.T_t) + x*rho_l);
		alpha_ip = x*rho_l/((1-x)*p_ip*mediaP.M_g/(R*wellP.T_t) + x*rho_l);
		
		Delta_mt = g*(m_gt+m_lt)/wellP.A_t;
		
		// Flow Pressure Drop
		mu_m = m_gt/(m_gt + m_lt)*(mediaP.mu_g) + (m_lt)/(m_gt + m_lt)*((1-reservoirP.r_wc)*mediaP.mu_o + reservoirP.r_wc*mediaP.mu_w);
		
		U_m = (w_lr + w_lp)*wellP.L_t/(2*m_lt) + (w_gr + w_gi + w_gp)*wellP.L_t/(2*m_gt);
		Re = (m_lt + m_gt)*wellP.D_t*U_m/(wellP.V_a*mu_m*wellP.A_t);
		f_D = 2*(-1.8*log10((wellP.epsilon_t/(3.7*wellP.D_t))^(1.11) + 6.9/Re))^(-2);
		
		Delta_ut = f_D*(m_gt + m_lt)*U_m^2/(2*wellP.A_t*wellP.D_t);
		
		
		// Pressures
		C =  -mediaP.M_g*rho_l*(g +f_D*U_m^2/(2*wellP.D_t));
		
		
		p_p = p_ip*exp((C*wellP.L_t + (1-x)*mediaP.M_g*(Delta_mt+Delta_ut))/(x*rho_l*R*wellP.T_t));
		p_p = p_ip - Delta_mt - Delta_ut;
		p_bh = p_ip;
		
		//p_ip = g*(m_gt+m_lt)/wellP.A_t/(1-exp((-wellP.L_t*g*mediaP.M_g*rho_l + (1-x)*mediaP.M_g*g*(m_gt+m_lt)/wellP.A_t)/(x*rho_l*R*wellP.T_t)));
		//p_p = p_ip - g*(m_gt+m_lt)/wellP.A_t;
		
		// Reservoir production
		w_lr = maxFunc(0.1,rho_l*reservoirP.Q_max*(1-(1-reservoirP.C)*(p_bh/reservoirP.p_r)-reservoirP.C*(p_bh/reservoirP.p_r)^2)); //A.2g
		w_gr = reservoirP.r_glr*w_lr; //A.2h	
		
		// Production Choke
		rho_p = p_p*mediaP.M_g*rho_l/((1-x)*p_p*mediaP.M_g + x*rho_l*R*wellP.T_t);
		w_p = ProductionChoke.w;		
		w_lp = (m_lt/(m_lt+m_gt))*w_p; //A.2d		
		w_gp = (m_gt/(m_lt+m_gt))*w_p; //A.2d		
		//w_lp = (1-alpha_wh)*w_p; //(m_lt/(m_lt+m_gt))*w_p; //A.2d		
		//w_gp = (alpha_wh)*w_p; //(m_lt/(m_lt+m_gt))*w_p; //A.2d		
		
		w_out[1] = w_gp;//(m_gt/(m_lt+m_gt))*w_p; //A.2c
		w_out[2] = (1-reservoirP.r_wc)*w_lp; //A.2e
		w_out[3] = reservoirP.r_wc*w_lp; //A.2d
			
		der(m_gt) = w_gr + w_gi - w_gp; //A.1b
		der(m_lt) = w_lr - w_lp;  //A.1c
		
	end Tubing2;
	
	model NewWell
		extends Interfaces.WellInterface; 
		parameter ComponentsProperties.MediaProperty mediaP;
		parameter ComponentsProperties.WellProperty wellP;
		parameter ComponentsProperties.ReservoirProperty reservoirP;
		parameter ComponentsProperties.WellValvesParameters valveP;
		
		Annulus ann(mediaP = mediaP, wellP = wellP, valveP = valveP);
		Tubing2 tub(mediaP = mediaP, wellP = wellP, reservoirP = reservoirP, valveP = valveP);
		
		SI.MassFlowRate w_in = inputConnector.massFlow[1]; 
		SI.MassFlowRate w_out[3] = outputConnector.massFlow;
		SI.Pressure p_in = inputConnector.pressure;
		SI.Pressure	p_out = outputConnector.pressure;
		
		equation
			connect(inputConnector, ann.inputConnector);
			connect(ann.outputConnector, tub.inputConnector);
			connect(tub.outputConnector, outputConnector);
	end NewWell;
	
	model Well "Well Model"
		//******************//***********************//
		//               Parameters					 //
		//******************//***********************//
			extends Interfaces.WellInterface; 
		  	  
			import SI = Modelica.SIunits;
			import R = Modelica.Constants.R;
			import g = Modelica.Constants.g_n;
			
			input Real u(min = 0, max = 0);
			
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
			//SI.MassFlowRate w_gl_max(min=0, start=0);
			
			SI.MassFlowRate w_in, w_out[3];
			SI.Pressure p_in, p_out;
			
			
			Valve GasLiftChoke(C = valveP.C_gl, rho = rho_gl, p_upstream= p_in, p_downstream = p_ta, w = w_in, u = u);
			Valve InjectionValve(C = valveP.C_iv, rho = rho_gi, p_upstream= p_ai, p_downstream = p_ti, w = w_gi, u = 1);
			Valve ProductionChoke(C = valveP.C_pc, rho = rho_p, p_upstream = p_p, p_downstream = p_out, w = w_p, u = 1);
		//******************//***********************//
		//               Equations					 //
		//******************//***********************//			
		equation 
			inputConnector.massFlow[1] = w_in;
			inputConnector.massFlow[2] = 0;
			inputConnector.massFlow[3] = 0;
			
			w_out = outputConnector.massFlow;
			p_in = inputConnector.pressure;
			p_out = outputConnector.pressure;
			
			
			
			p_ai = ((R*wellP.T_a)/(wellP.V_a*mediaP.M_g) + g/(2*wellP.A_a))*m_ga;  //A.3a
			p_p = (R*wellP.T_t*m_gt)/(mediaP.M_g*wellP.V_t-mediaP.M_g*(1/rho_Lw)*m_lt) - (g*(m_gt+m_lt))/(2*wellP.A_t); //A.3b (3.4.18 - Binder)
			p_ti = p_p + (g*(m_lt+m_gt))/wellP.A_t;  //A.3c
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
			
			der(m_ga) = w_in - w_gi; //A.1a
			der(m_gt) = w_gr + w_gi-w_out[1]; //A.1b
			der(m_lt) = w_lr - w_lp;  //A.1c
	end Well;
	
end newComponents;