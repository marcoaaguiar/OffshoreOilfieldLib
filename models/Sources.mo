within Network;

package Sources
	import SI = Modelica.SIunits;
	block GasPressureSource "Source with fixed gas and pressure"
		extends Interfaces.SourceInterface;
		SI.Pressure pressure;
		SI.MassFlowRate	gasMassFlow;
		equation
			outputConnector.pressure = pressure;
			outputConnector.massFlow[1] = gasMassFlow;
			outputConnector.massFlow[2] = 0;
			outputConnector.massFlow[3] = 0;
	end GasSource;
	
	block GasSource "Source with fixed gas and pressure"
		extends Interfaces.SourceInterface;
		SI.MassFlowRate	gasMassFlow;
		equation
			outputConnector.massFlow[1] = gasMassFlow;
			outputConnector.massFlow[2] = 0;
			outputConnector.massFlow[3] = 0;
	end GasSource;
	
	block PressureSource "Source with fixed gas and pressure"
		extends Interfaces.SourceInterface;
		SI.Pressure pressure;
		equation
			outputConnector.pressure = pressure;
	end GasSource;
	
	block MultiPhaseSource "Source with fixed gas, oil, and water"
		extends Interfaces.SourceInterface;
		SI.MassFlowRate	gasMassFlow;
		SI.MassFlowRate	oilMassFlow;
		SI.MassFlowRate	waterMassFlow;
		equation
			outputConnector.massFlow[1] = gasMassFlow;
			outputConnector.massFlow[2] = oilMassFlow;
			outputConnector.massFlow[3] = waterMassFlow;
	end GasSource;
end Sources;
