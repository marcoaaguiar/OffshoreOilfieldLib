within Network;
  
package Functions
  function maxFunc "Mock function for the built-in max() function"
			input Real x1, x2; 
			output Real y; 
		algorithm  
			y := max(x1,x2);
	end maxFunc; 
	
	function heaviside "Heaviside (Step) Function"
			input Real x;
			output Real y;
		algorithm
			y:=1/(1+exp(-2*1*x));
	end heaviside;
end Functions;