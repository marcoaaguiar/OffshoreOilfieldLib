within Network;

package Sinks
	import SI = Modelica.SIunits;
	block PressureSink
		extends Interfaces.SinkInterface;
		SI.Pressure pressure;
		equation
			inputConnector.pressure = pressure;
	end GasSource;
end Sources;
