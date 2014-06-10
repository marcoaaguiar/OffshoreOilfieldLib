model oldSinglePipelineRiser
	networkModels.networkComponents.pipeline pipeline;
	
	Modelica.SIunits.Pressure p_out = 50e5;
	Modelica.SIunits.MassFlowRate w_g = 24;
	Modelica.SIunits.MassFlowRate w_o = 10;
	Modelica.SIunits.MassFlowRate w_w = 16;
	
	equation
		pipeline.p_out = p_out;
		pipeline.w_in[1] = w_g;
		pipeline.w_in[2] = w_o;
		pipeline.w_in[3] = w_w;
end oldSinglePipelineRiser;
	