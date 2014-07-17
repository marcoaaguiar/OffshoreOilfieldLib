model newTubingTest
	newComponents.Tubing2 tubing;
	Sources.GasSource gSource;
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_out = 1.5158e7-0.5e7;
	Modelica.SIunits.MassFlowRate w_gl = 2;
	equation
		connect(gSource.outputConnector, tubing.inputConnector);
		connect(tubing.outputConnector, pSink.inputConnector);
		gSource.gasMassFlow = w_gl;
		pSink.pressure = p_out;
end newTubingTest;
	