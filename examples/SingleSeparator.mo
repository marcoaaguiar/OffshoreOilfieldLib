model SingleSeparator
	Components.Separator separator;
	Sources.MultiPhaseSource mpSource;
	Sinks.PressureSink pSink[3];
	
	Modelica.SIunits.Pressure p_out = 50e5;
	Modelica.SIunits.MassFlowRate w_g = 24;
	Modelica.SIunits.MassFlowRate w_o = 10;
	Modelica.SIunits.MassFlowRate w_w = 16;
	
	equation
		connect(mpSource.outputConnector, separator.inputConnector);
		connect(separator.outputConnector[1], pSink[1].inputConnector);
		connect(separator.outputConnector[2], pSink[2].inputConnector);
		connect(separator.outputConnector[3], pSink[3].inputConnector);
		
		mpSource.gasMassFlow = w_g;
		mpSource.oilMassFlow = w_o;
		mpSource.waterMassFlow = w_w;
		
		pSink[1].pressure = p_out;
		pSink[2].pressure = p_out;
		pSink[3].pressure = p_out;
end SingleSeparator;
	