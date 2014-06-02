within Network;

package ComponentsProperties
import SI = Modelica.SIunits;
  
	record NetworkProperty "Network Configuration and Properties"
	  parameter Integer n_wells = 8;  // Number of Wells
		parameter Integer n_manifolds = 2; //Number of Manifolds
		parameter Integer n_pipelines = 2; //Number of Pipeline-risers
		parameter Integer n_separators = 2; //Number of Separators
		parameter Integer n_compressors = 2; // Number of Compressors
		
		//parameter Integer routingTable[n_manifolds,n_wells] = [1,1,1,1,0,0,0,0;	 0,0,0,0,1,1,1,1];
	end NetworkProperty;

	record miniNetworkProperty "Network Configuration and Properties"
	  parameter Integer n_wells = 4;  // Number of Wells
		parameter Integer n_manifolds = 1; //Number of Manifolds
		parameter Integer n_pipelines = 1; //Number of Pipeline-risers
		parameter Integer n_separators = 1; //Number of Separators
		parameter Integer n_compressors = 1; // Number of Compressors
		
		//parameter Integer routingTable[n_manifolds,n_wells] = [1,1,1,1,0,0,0,0;	 0,0,0,0,1,1,1,1];
	end miniNetworkProperty;
	record MediaProperty "Properties of the flowing media"
		parameter SI.Density rho_w = 1030;       // [kg/m^3]		// water density	   
		parameter SI.Density rho_o = 930;        // [kg/m^3]		// oil density	   
		parameter SI.MolarMass M_g = 0.0195;       // [kg/mol]	// lift gas molecular weight (Molar Mass)  
	end MediaProperty;
	
	record WellProperty
	  parameter SI.Temperature T_a=350;          // [K]	// annular temperature
		parameter SI.Temperature T_t=350;          // [K]	// tubing temperature	   
		parameter SI.Area A_a=0.02;         // [m^2]	// annular cross sectional area	   
		parameter SI.Area A_t=0.012;        // [m^2]	// tubing cross sectional area	   
		parameter SI.Volume V_a = 30;               // [m^3]		// Annular volume	  
		parameter SI.Length L_a = V_a/A_a; 
		parameter SI.Volume V_t=  18;               // [m^3]		// Tubing volume
		parameter SI.Length L_t = V_t/A_t; 
		parameter SI.Length L_w=  400;          // [m]		// tubing length from the injection point to the reservoir	  
	end WellProperty;
		
	record ReservoirProperty
	  parameter SI.Pressure p_r(nominal=1e5)=25e6;         // [Pa]		// reservoir pressure   
		parameter SI.VolumeFraction r_wc= 0.4;          // [-]		// well water cut. Volume or Weight?
		parameter Real r_gor=0.08;     // [-]		// well gor	   
		parameter Real r_glr=(1-r_wc)*r_gor;		// Gas-to-liquid ratio in well	
		parameter Real C=0.8;            // [-]	 // Quadratic term of the Volger equation
		parameter Real Q_max=0.025;   // [m^3/s]		// Empirical constant representing the theoretical absolute open flow (AOF)	  
	end ReservoirProperty;
		
		
	record WellValvesParameters
		parameter Real C_iv = 0.00016;   // [m^2]	// gas injection valve constant	   
		parameter Real C_pc = 0.0014;    // [m^2]	// production choke valve constant	   
		parameter Real C_gl = 0.00016;
	end WellValvesParameters;
	
	record pipelineParameters
		import SI = Modelica.SIunits;
		
		parameter SI.Length L_p = 13000;          // [m]// horizontal pipeline liength
		parameter SI.Radius r_p = 0.10;           // [m] // pipeline radius
		parameter SI.Diameter D_p = 2*r_p;        // [m]	// pipeline diameter	
		parameter SI.Area A_p = Modelica.Constants.pi*r_p^2;     // [m^2]	// pipeline cross sectional area	   
		parameter SI.Volume V_p = L_p*A_p;    // [m^3]	// horizontal pipeline volume	   
		parameter SI.Temperature T_p = 330;            // [K]	// pipeline temperatur	 
		     	
		parameter SI.Temperature T_r = 330;            // [K]	// riser temperature	   	
		parameter SI.Length L_r = 600;            // [m]	// Riser length	   	
		parameter SI.Radius r_r = 0.10;           // [m]	// Riser radius	   	
		parameter SI.Diameter D_r = 2*r_r;        // [m]	// Riser diameter	   	
		parameter SI.Area A_r = Modelica.Constants.pi*r_r^2;     // [m^2] 	// Riser sectional area	   
		parameter SI.Volume V_r = A_r*L_r;    // [m^3]	// Riser volume	   	
		parameter SI.Length epsilon = 2.8e-5;     // [m]	// Pipeline roughness
		   	   
		parameter SI.MolarMass M_g=0.0195;       // [kg/mol]	// gas molar    	
		parameter SI.Density rho_w = 1030;       // [kg/m^3]	// water density	   
		parameter SI.Density rho_o = 930;        // [kg/m^3]	// oil density	  
		 
		parameter SI.DynamicViscosity mu_w = 8.94e-4;       // [Pa*s]	// water dynamic viscosity	   	
		parameter SI.DynamicViscosity mu_o = 1e-4;          // [Pa*s]	// oil dynamic viscosity	  
		
		parameter Real K_gp = 0.03;         // [-]	// Orifice equation parameters	   
		parameter Real K_lp = 1;            // [-]	
		parameter Real K_gr = 0.0034;	           // [m^2]	// Valve coefficient for gas/liquid flow through the choke valve at top of riser 
		parameter Real K_lr = 0.0014;	           // [m^2]	   
	end pipelineParameters;
end ComponentsProperties;
