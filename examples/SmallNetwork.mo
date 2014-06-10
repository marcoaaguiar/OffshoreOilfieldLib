model SmallNetwork
	parameter Integer n_wells = 4;
	
	ComponentsApproximation.Well well[n_wells];
	ComponentsApproximation.ProductionManifold manifold(n_inputs = n_wells);
	ComponentsApproximation.Pipeline pipeline;
	ComponentsApproximation.Separator separator;
	ComponentsApproximation.Compressor compressor;
	
	Sources.GasPressureSource gSource[n_wells];
	Sinks.PressureSink pSink[3];
	
	Modelica.SIunits.Pressure p_sep = 50e5; //Separator Pressure
	Modelica.SIunits.Pressure p_exp = 200e5; //Exportation Line Pressure
	Modelica.SIunits.MassFlowRate w_gl = 1;
	Modelica.SIunits.MassFlowRate w_sl = 10;
	
	equation
		connect(gSource[1:n_wells].outputConnector, well[1:n_wells].inputConnector);
		connect(well[1:n_wells].outputConnector, manifold.inputConnector[1:n_wells]);
		connect(manifold.outputConnector, pipeline.inputConnector);
		connect(pipeline.outputConnector, separator.inputConnector);
		connect(separator.outputConnector[1], compressor.inputConnector);
		connect(separator.outputConnector[2], pSink[2].inputConnector);
		connect(separator.outputConnector[3], pSink[3].inputConnector);
		connect(compressor.outputConnector, pSink[1].inputConnector);
		
		
		for i in 1:n_wells loop
			gSource[i].pressure = p_exp;
			gSource[i].gasMassFlow = w_gl;
		end for;
		
		pSink[1].pressure = p_exp;
		pSink[2].pressure = p_sep;
		pSink[3].pressure = p_sep;
		
		compressor.w_sl = w_sl;
end SmallNetwork;
	