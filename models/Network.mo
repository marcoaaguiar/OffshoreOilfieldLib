package Network
	import SI = Modelica.SIunits;
	constant Real k_g = 1;
	constant Real k_o = 10;
	constant Real k_w = -1;
	constant Real k_fl = -1000;
	constant Real k_gl = -1;
	
	model network
		replaceable package components = Components;
		
		parameter SI.Pressure p_export(nominal=1e5) = 200e5;
		ComponentsProperties.NetworkProperty networkP;
		
		//Network Input
		input SI.MassFlowRate w_gl[networkP.n_wells](each min=-1/1000, each max=1/1000, each free=true);	
		input SI.MassFlowRate w_fl[networkP.n_compressors](each min=0, each start=0, each free=true);
		input SI.MassFlowRate w_sl[networkP.n_compressors](each min=0);
		input SI.Pressure p_s[networkP.n_separators](each nominal=1e5, each min=50e5, each max=60e5, each start=50e5);
		
		Real int_w_gl[networkP.n_wells] (each min=1, each max=3, start = 1.5, fixed = true);	
		//Network Output
		SI.MassFlowRate w_out[3](each min=0, each start =0);
		
		//Network Internal Variables
		SI.MassFlowRate w_gm(min=0, start = 0);
		
		//System Components
		components.Well w[networkP.n_wells];
		components.ProductionManifold m[networkP.n_pipelines] (each n_inputs=4);
		components.Pipeline p[networkP.n_pipelines];
		components.Separator s[networkP.n_separators];
		components.Flare fl[networkP.n_separators];
		components.Compressor c[networkP.n_compressors];
		equation 
			// Well -> Manifold
			for i in 1:4 loop
				der(int_w_gl[i]) = w_gl[i];
				connect(w[i].outputConnector, m[1].inputConnector[i]);
				w[i].inputConnector.massFlow[1] = int_w_gl[i];
				w[i].inputConnector.massFlow[2] = 0;
				w[i].inputConnector.massFlow[3] = 0;
				w[i].inputConnector.pressure = p_export;
			end for;
			
			for i in 5:8 loop
				der(int_w_gl[i]) = w_gl[i];
				connect(w[i].outputConnector, m[2].inputConnector[i-4]);
				w[i].inputConnector.massFlow[1] = int_w_gl[i];
				w[i].inputConnector.massFlow[2] = 0;
				w[i].inputConnector.massFlow[3] = 0;
				w[i].inputConnector.pressure = p_export;
			end for;

			// Manifold -> Pipeline-Riser		    
			for i in 1:networkP.n_pipelines loop
				connect(m[i].outputConnector, p[i].inputConnector);
			end for;
		  
			// Pipeline-Riser -> Separator
			for i in 1:networkP.n_separators loop
				connect(p[i].outputConnector, s[i].inputConnector);
				s[i].inputConnector.pressure = p_s[i];
			end for;
			
			// Separator -> Flare
			for i in 1:networkP.n_separators loop
			    connect(s[i].outputConnector[1], fl[i].inputConnector);
			    fl[i].w_fl = w_fl[i];
			end for;
			
			// Flare -> Compressor
			for i in 1:networkP.n_compressors loop
			   connect(fl[i].outputConnector, c[i].inputConnector);
			   c[i].w_sl = w_sl[i];
			   c[i].outputConnector.pressure = p_export;
			end for;
					
			w_gm = sum(int_w_gl[i] for i in 1:networkP.n_wells);
			//Network Outputs
			w_out[1] = sum(c[i].outputConnector.massFlow[1] for i in 1:networkP.n_compressors)- w_gm;  //-sum (w_sl[i] for i in 1:networkP.n_compressors);
			w_out[2] = sum(s[i].outputConnector[2].massFlow[2] for i in 1:networkP.n_separators);
			w_out[3] = sum(s[i].outputConnector[3].massFlow[3] for i in 1:networkP.n_separators);
			//s[1].outputConnector[2].massFlow[2] = w_out[2];
			//s[1].outputConnector[3].massFlow[3] = w_out[3];
	end network;
								
	model networkApproximation
		extends network;
		redeclare package components = ComponentsApproximation;			
	end networkApproximation;
						
	model miniNetwork
		ComponentsProperties.NetworkProperty networkP(n_wells = 2, n_pipelines=1, n_separators =1, n_compressors = 1);

		Components.Well w[networkP.n_wells];
		Components.ProductionManifold m(n_inputs = networkP.n_wells);
		Components.Pipeline p;
		//Components.Separator s;
		//Components.Flare fl;
		//Components.Compressor c;
		
		input Real w_gl[networkP.n_wells] (each min =0);
		//input Real w_fl (min =0);
		//input Real w_sl (min =0);
		equation
		for i in 1:networkP.n_wells loop
			w[i].inputConnector.pressure =200e5;
			w[i].inputConnector.massFlow[1] = w_gl[i];
			w[i].inputConnector.massFlow[2] =0;
			w[i].inputConnector.massFlow[3] =0;
			connect(w[i].outputConnector, m.inputConnector[i]);
	    end for;
	    connect(m.outputConnector, p.inputConnector);
	    //connect(p.outputConnector, s.inputConnector);
	    //connect(s.outputConnector, fl.inputConnector);
	    //connect(fl.outputConnector, c.inputConnector);
	    
	    //c.outputConnector.pressure =200e5;
	    p.outputConnector.pressure = 60e5;

	    //fl.w_fl = w_fl;
	    //c.w_sl = w_sl;
	  
	end miniNetwork;
	
	model miniNetworkApproximation
		extends miniNetwork;
		redeclare package components = ComponentsApproximation;			
	end miniNetworkApproximation;						
end Network;						
															
																
																	
																		
																			
																				
																					
																						