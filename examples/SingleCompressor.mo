model SingleCompressor
	Components.Compressor compressor;
	Sources.GasPressureSource gSource;
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_in = 50e5;
	Modelica.SIunits.Pressure p_out = 200e5;
	
	Modelica.SIunits.MassFlowRate w_g = 24;
	Modelica.SIunits.MassFlowRate w_sl = 10;
	
	equation
		connect(gSource.outputConnector, compressor.inputConnector);
		connect(compressor.outputConnector, pSink.inputConnector);

		gSource.gasMassFlow = w_g;
		gSource.pressure = p_in;

		compressor.w_sl = w_sl;
		
		pSink.pressure = p_out;
end SingleCompressor;
	