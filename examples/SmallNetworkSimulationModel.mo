model SmallNetworkSimulationModel
	parameter Integer n_wells = 2;

	ComponentsApproximation.Well well[n_wells]; 
	ComponentsApproximation.ProductionManifold manifold(n_inputs = n_wells);
	ComponentsApproximation.Pipeline pipeline; 
	ComponentsApproximation.Separator separator;
	ComponentsApproximation.Compressor compressor;
	
	Sources.GasPressureSource gSource[n_wells];
	Sinks.PressureSink pSink[3];
	
	parameter Modelica.SIunits.Pressure p_exp = 200e5; //Exportation Line Pressure
	
	input Modelica.SIunits.Pressure p_sep(start = 50e5, max = 60e5, min = 50e5); //Separator Pressure
	input Modelica.SIunits.MassFlowRate w_gl[n_wells](each min = 0, each max = 10, each start = 1);
	input Modelica.SIunits.MassFlowRate w_sl(start = 10, min =0);
	
	Real cost(start =0, fixed = true, min = 0);
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
			gSource[i].gasMassFlow = w_gl[i];
		end for;
		
		pSink[1].pressure = p_exp;
		pSink[2].pressure = p_sep;
		pSink[3].pressure = p_sep;
		
		compressor.w_sl = w_sl;
		der(cost) = (pSink[2].inputConnector.massFlow[2]-25)^2;
end SmallNetworkSimulationModel;
