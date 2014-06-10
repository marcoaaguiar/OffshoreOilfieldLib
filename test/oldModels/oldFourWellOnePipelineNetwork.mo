model oldFourWellOnePipelineNetwork
	networkModels.networkComponents.well well[4];
	networkModels.networkComponents.pipeline pipeline;
		
	Modelica.SIunits.Pressure p_in = 200e5;
	Modelica.SIunits.Pressure p_out = 50e5;
	Modelica.SIunits.MassFlowRate w_gl = 0.5;
	
	equation
		for i in 1:4 loop
			well[i].p_in = p_in;
			well[i].w_in = w_gl + 0.5*i;
			well[i].p_out = pipeline.p_in;
		end for;
		pipeline.w_in = sum(well[i].w_out for i in 1:4);
		
		pipeline.p_out = p_out;
end oldFourWellOnePipelineNetwork;
	