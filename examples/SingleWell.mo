model SingleWell
	Components.Well well(valveP = ComponentsProperties.WellValvesParameters(C_pc = 0.0016));
	Sources.GasPressureSource gSource;
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_in = 200e5;
	Modelica.SIunits.Pressure p_out = 90e5;
	//Modelica.SIunits.MassFlowRate w_gl;
	Modelica.Blocks.Sources.Step step(offset = 1, height = -0.5, startTime = 5000);
	
	equation
		connect(gSource.outputConnector, well.inputConnector);
		connect(well.outputConnector, pSink.inputConnector);
		gSource.pressure = p_in;
		gSource.gasMassFlow = step.y;
		pSink.pressure = p_out;
end SingleWell;
	