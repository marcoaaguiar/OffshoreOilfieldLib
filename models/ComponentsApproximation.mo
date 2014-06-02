within Network;

package ComponentsApproximation
	extends Components; 

	function maxFuncApproximation
		input Real x1, x2;
		output Real y;
	algorithm
		y := ((x1-x2)^2+Modelica.Constants.eps^2)^0.5/2+(x1+x2)/2;
		//y := ((x1-x2)^2+0.1^2)^0.5/2+(x1+x2)/2;
	end maxFuncApproximation;

	redeclare function maxFunc =  maxFuncApproximation;
end ComponentsApproximation;
