model newWellTest
	newComponents.NewWell nWell(wellP.T_a = 300,
								wellP.T_t = 300,
								wellP.L_w = 0,
								wellP.A_t = 0.0094,
								wellP.A_a = 0.0257,
								wellP.V_a = 64.25,
								wellP.V_t = 23.5,
								mediaP.rho_o = 921.8241,
								mediaP.rho_w = 1000,
								mediaP.M_g = 0.016,
								reservoirP.r_wc = 0.15,
								reservoirP.p_r = 18e6,
								reservoirP.r_gor = 0.0176,
								valveP.C_iv = 0.000135,
								valveP.C_pc = 0.0007,
								ann.m_ga0 = 7000
								);
	Sources.GasPressureSource gpSource;
	Sinks.PressureSink pSink;
	
	Modelica.SIunits.Pressure p_in = 190e5;//200e5;
	Modelica.SIunits.Pressure p_out = 10e5;//1.5158e7-0.5e7;
	Modelica.SIunits.MassFlowRate w_gl = 0.3;
	equation
		connect(gpSource.outputConnector, nWell.inputConnector);
		connect(nWell.outputConnector, pSink.inputConnector);
		gpSource.pressure = p_in;
		gpSource.gasMassFlow = w_gl;
		pSink.pressure = p_out;
end newWellTest;
	