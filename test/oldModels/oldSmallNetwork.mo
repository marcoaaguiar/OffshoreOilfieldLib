model oldSmallNetwork
	networkModels.networkComponents.well well[4];
	networkModels.networkComponents.pipeline pipeline;
	networkModels.networkComponents.separator separator;
	networkModels.networkComponents.compressor compressor;
		
	Modelica.SIunits.Pressure p_sep = 50e5;
	Modelica.SIunits.Pressure p_exp = 200e5;
	Modelica.SIunits.MassFlowRate w_gl = 1;
	Modelica.SIunits.MassFlowRate w_sl = 10;

	equation
		for i in 1:4 loop
			well[i].p_in = p_exp;
			well[i].w_in = w_gl;
			well[i].p_out = pipeline.p_in;
		end for;
		pipeline.w_in = sum(well[i].w_out for i in 1:4);
		
		pipeline.w_out = separator.w_in;
		pipeline.p_out = separator.p_in;
		
		separator.w_out[1] = compressor.w_in;
		separator.p_out = compressor.p_in;
		separator.p_out = p_sep;
		
		compressor.p_out = p_exp;
		compressor.w_sl = w_sl;

end oldSmallNetwork;
	