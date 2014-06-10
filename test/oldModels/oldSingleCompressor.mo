model oldSingleCompressor
	networkModels.networkComponents.compressor compressor;
	
	Modelica.SIunits.Pressure p_in = 50e5;
	Modelica.SIunits.Pressure p_out = 200e5;
	
	Modelica.SIunits.MassFlowRate w_g = 24;
	Modelica.SIunits.MassFlowRate w_sl = 10;
	
	equation
		compressor.p_in = p_in;
		compressor.p_out = p_out;
		compressor.w_in = w_g;
		
		compressor.w_sl = w_sl;

end oldSingleCompressor;
	