model FourWellOnePipelineNetwork
	Components.Well well[4];
	Components.ProductionManifold manifold;
	Components.Pipeline pipeline;
	
	Sources.GasPressureSource gSource[4];
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_in = 200e5;
	Modelica.SIunits.Pressure p_out = 50e5;
	Modelica.SIunits.MassFlowRate w_gl = 0.5;
	
	equation
		connect(gSource[1:4].outputConnector, well[1:4].inputConnector);
		connect(well[1:4].outputConnector, manifold.inputConnector[1:4]);
		connect(manifold.outputConnector, pipeline.inputConnector);
		connect(pipeline.outputConnector, pSink.inputConnector);
		
		for i in 1:4 loop
			gSource[i].pressure = p_in;
			gSource[i].gasMassFlow = w_gl+0.5*i;
		end for;
		pSink.pressure = p_out;
end FourWellOnePipelineNetwork;
	