model SinglePipelineRiser
	Components.Pipeline pipeline;
	Sources.MultiPhaseSource mpSource;
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_out = 50e5;
	Modelica.SIunits.MassFlowRate w_g = 24;
	Modelica.SIunits.MassFlowRate w_o = 10;
	Modelica.SIunits.MassFlowRate w_w = 16;
	
	equation
		connect(mpSource.outputConnector, pipeline.inputConnector);
		connect(pipeline.outputConnector, pSink.inputConnector);
		
		mpSource.gasMassFlow = w_g;
		mpSource.oilMassFlow = w_o;
		mpSource.waterMassFlow = w_w;
		pSink.pressure = p_out;
end SinglePipelineRiser;
	