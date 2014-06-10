model oldSingleWellNetwork
	networkModels.networkComponents.well well;
	
	Modelica.SIunits.Pressure p_in = 200e5;
	Modelica.SIunits.Pressure p_out = 90e5;
	Modelica.SIunits.MassFlowRate w_gl = 1.5;
	equation
		well.w_in = w_gl;
		well.p_in = p_in;
		well.p_out = p_out;
end oldSingleWellNetwork;
	