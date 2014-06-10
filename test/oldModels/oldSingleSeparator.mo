model oldSingleSeparator
	networkModels.networkComponents.separator separator;
	
	Modelica.SIunits.Pressure p_out = 50e5;
	Modelica.SIunits.MassFlowRate w_g = 24;
	Modelica.SIunits.MassFlowRate w_o = 10;
	Modelica.SIunits.MassFlowRate w_w = 16;
	
	equation
		separator.p_out = p_out;
		separator.w_in[1] = w_g;
		separator.w_in[2] = w_o;
		separator.w_in[3] = w_w;
end oldSingleSeparator;
	