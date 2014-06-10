model oldSmallNetworkSimulationModel
	parameter Integer n_wells = 2;

	networkModels.networkComponentsApproximation.well well[n_wells]; 
	networkModels.networkComponentsApproximation.pipeline pipeline; 
	networkModels.networkComponentsApproximation.separator separator;
	networkModels.networkComponentsApproximation.compressor compressor;
	
	
	parameter Modelica.SIunits.Pressure p_exp = 200e5; //Exportation Line Pressure
	
	input Modelica.SIunits.Pressure p_sep(start = 50e5, max = 60e5, min = 50e5); //Separator Pressure
	input Modelica.SIunits.MassFlowRate w_gl[n_wells](each min = 0, each max = 10, each start = 1);
	input Modelica.SIunits.MassFlowRate w_sl(start = 10, min =0);
	
	Real cost(start =0, fixed = true, min = 0);
	equation
		
		
		for i in 1:n_wells loop
			well[i].p_in = p_exp;
			well[i].w_in = w_gl[i];
			well[i].p_out = pipeline.p_in;
		end for;
		pipeline.w_in = sum(well[i].w_out for i in 1:n_wells);
		
		pipeline.w_out = separator.w_in;
		pipeline.p_out = separator.p_in;
		
		separator.p_in = p_sep;
		separator.w_out[1] = compressor.w_in;
		separator.p_out = compressor.p_in;
		
		compressor.w_sl = w_sl;
		compressor.p_out = p_exp;
		der(cost) = (separator.w_out[2]-25)^2;
end oldSmallNetworkSimulationModel;
