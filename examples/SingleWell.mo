model SingleWell
	Components.Well well;
	Sources.GasPressureSource gSource;
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_in = 200e5;
	Modelica.SIunits.Pressure p_out = 90e5;
	Modelica.SIunits.MassFlowRate w_gl = 1.5;
	
	equation
		connect(gSource.outputConnector, well.inputConnector);
		connect(well.outputConnector, pSink.inputConnector);
		gSource.pressure = p_in;
		gSource.gasMassFlow = w_gl;
		pSink.pressure = p_out;
end SingleWell;
	