model newAnnulusTest
	newComponents.Annulus annulus;
	Sources.GasPressureSource gpSource;
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_in = 200e5;
	Modelica.SIunits.Pressure p_out = 1.5158e7;
	Modelica.SIunits.MassFlowRate w_gl = 2;
	equation
		connect(gpSource.outputConnector, annulus.inputConnector);
		connect(annulus.outputConnector, pSink.inputConnector);
		gpSource.pressure = p_in;
		gpSource.gasMassFlow = w_gl;
		pSink.pressure = p_out;
end newAnnulusTest;
	