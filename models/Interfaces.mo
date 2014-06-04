within Network;
package Interfaces
	import SI = Modelica.SIunits;

	connector PipeConnector
		SI.Pressure pressure(nominal=1e5, min=0, start= 1e5);
		SI.MassFlowRate[3] massFlow(each min=0, each start=1);
	end PipeConnector;
	
	partial block BaseInterface 
		input PipeConnector inputConnector;
		output PipeConnector outputConnector;
	end BaseInterface;

	partial block SourceInterface
		output PipeConnector outputConnector;
	end SourceInterface;	
	
	partial block SinkInterface
		input PipeConnector inputConnector;
	end SinkInterface;
	
	partial block valve = BaseInterface;
		   
	partial model WellInterface = BaseInterface;
	
	partial model PipelineRiserInterface = BaseInterface;
	
	partial model SeparatorInterface 
		input PipeConnector inputConnector;
		output PipeConnector outputConnector[3]; //One for each outlet
	end SeparatorInterface;
	
	partial model CompressorInterface = BaseInterface;
	
	partial model FlareInterface = BaseInterface;
	
	partial model ProductionManifoldInterface
	  parameter Integer n_inputs = 4;
		input PipeConnector inputConnector[n_inputs];
		output PipeConnector outputConnector;
	end ProductionManifoldInterface;
	/// To be deprecated

	connector pressurePortIn = input SI.Pressure(nominal=1e5, min=0);
	connector pressurePortOut = output SI.Pressure(nominal=1e5, min=0);
	connector flowPortIn = input SI.MassFlowRate(nominal=1, min=0, start=1);
	connector threePhaseFlowPortIn = flowPortIn[3];
	connector flowPortOut = output SI.MassFlowRate(nominal=1, min=0, start=1);
	connector threePhaseFlowPortOut = flowPortOut[3];

end Interfaces;