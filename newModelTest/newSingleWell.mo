model newSingleWell
	newComponents.Well well;
	Sources.PressureSource pSource;
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_in = 200e5;
	Modelica.SIunits.Pressure p_out = 50e5;
	Modelica.SIunits.MassFlowRate w_gl = 1.5;
	Modelica.Blocks.Sources.Step step(offset = 0.5, height = -0.4, startTime = 7000);
	Modelica.Blocks.Sources.Step step2(offset = 0, height = -0.09, startTime = 14000);
	equation
		connect(pSource.outputConnector, well.inputConnector);
		connect(well.outputConnector, pSink.inputConnector);
		pSource.pressure = p_in;
		pSink.pressure = p_out;
		well.u = step.y + step2.y;
end newSingleWell;
	