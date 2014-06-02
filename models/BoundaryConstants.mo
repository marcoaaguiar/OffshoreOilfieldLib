within Network

package BoundaryConstants "Constants defined to avoid instabilities in Simulations and Jacobians/Hessians"
	constant Real SqrtEpsilon = 1e-1; 
	constant SI.Mass minMass = 1;
	constant SI.Volume minVolume = 1e-1;
	constant SI.Pressure minPressure = 1e5;
	constant SI.Density minDensity = 1;
	constant SI.Velocity minVelocity = 1e-3;
	constant SI.Angle angleEpsilon = 1e-3;
	
end BoundaryConstants;