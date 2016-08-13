function Eight_Node_Brick_Inverse = Eight_Node_Brick_Inverse()

	% shape function coordinates  coefficients
	Brick_Coordinates(1 ,:)=[ 1, 1, 1]; 
	Brick_Coordinates(2 ,:)=[-1, 1, 1]; 
	Brick_Coordinates(3 ,:)=[-1,-1, 1]; 
	Brick_Coordinates(4 ,:)=[ 1,-1, 1]; 
	Brick_Coordinates(5 ,:)=[ 1, 1,-1]; 
	Brick_Coordinates(6 ,:)=[-1, 1,-1]; 
	Brick_Coordinates(7 ,:)=[-1,-1,-1]; 
	Brick_Coordinates(8 ,:)=[ 1,-1,-1]; 

	% gauss coordinates coefficients 
	Brick_8_Gauss_Coordinates(1,:)=[-1,-1,-1];
	Brick_8_Gauss_Coordinates(2,:)=[-1,-1, 1];
	Brick_8_Gauss_Coordinates(3,:)=[-1, 1,-1];
	Brick_8_Gauss_Coordinates(4,:)=[-1, 1, 1];
	Brick_8_Gauss_Coordinates(5,:)=[ 1,-1,-1];
	Brick_8_Gauss_Coordinates(6,:)=[ 1,-1, 1];
	Brick_8_Gauss_Coordinates(7,:)=[ 1, 1,-1];
	Brick_8_Gauss_Coordinates(8,:)=[ 1, 1, 1];

	Eight_Node_Brick_Inverse = zeros(8,8);

	SQRT_3 = sqrt(1.0/3.0);

	for j=1:8	% Gauss Point Coordinates
		for i=1:8 % Shape Functions 
			Eight_Node_Brick_Inverse(j,i) = 1/8*(1+Brick_Coordinates(i,1)*Brick_8_Gauss_Coordinates(j,1)*SQRT_3)* ...
												(1+Brick_Coordinates(i,2)*Brick_8_Gauss_Coordinates(j,2)*SQRT_3)* ...
												(1+Brick_Coordinates(i,3)*Brick_8_Gauss_Coordinates(j,3)*SQRT_3);
		end
	end

	Eight_Node_Brick_Inverse

	Eight_Node_Brick_Inverse = inv(Eight_Node_Brick_Inverse);

end