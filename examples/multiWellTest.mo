model multiWellTest
	parameter Integer n_wells = 100;
	Components.Well well[n_wells];
	Sources.GasPressureSource gSource[n_wells];
	Sinks.PressureSink pSink[n_wells];
	
	Modelica.SIunits.Pressure p_in = 200e5;
	Modelica.SIunits.Pressure p_out = 90e5;
	Modelica.SIunits.MassFlowRate w_gl_step = 0.2;
	Modelica.SIunits.MassFlowRate w_gl_base = 0.1;
	
	equation
		for w in 1:n_wells loop
			connect(gSource[w].outputConnector, well[w].inputConnector);
			connect(well[w].outputConnector, pSink[w].inputConnector);
			gSource[w].pressure = p_in;
			gSource[w].gasMassFlow = w_gl_base+w_gl_step*(w-1);
			pSink[w].pressure = p_out;
		end for;
end multiWellTest;
	