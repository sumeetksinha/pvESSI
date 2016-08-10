function Twenty_Node_Brick_Inverse = Twenty_Node_Brick_Inverse()

	Twenty_Brick_Coordinates = Build_Brick_Coordinates();
	Twenty_Brick_Coordinates = Twenty_Brick_Coordinates(1:20,:);

	Twenty_Node_Brick_Inverse = zeros(20,20);

	SQRT_3_5 = sqrt(3/5);

	for j=1:20
		for i=1:8
			Twenty_Node_Brick_Inverse(j,i) = 1/8*(1+Twenty_Brick_Coordinates(j,1)*Twenty_Brick_Coordinates(i,1)*SQRT_3_5)* ...
												(1+Twenty_Brick_Coordinates(j,2)*Twenty_Brick_Coordinates(i,2)*SQRT_3_5)* ...
												(1+Twenty_Brick_Coordinates(j,3)*Twenty_Brick_Coordinates(i,3)*SQRT_3_5)* ...
												(-2+ (Twenty_Brick_Coordinates(j,1)*Twenty_Brick_Coordinates(i,1)*SQRT_3_5)+ ...
												     (Twenty_Brick_Coordinates(j,2)*Twenty_Brick_Coordinates(i,2)*SQRT_3_5)+ ...
												     (Twenty_Brick_Coordinates(j,3)*Twenty_Brick_Coordinates(i,3)*SQRT_3_5));
		end
		for i=[9,11,13,15]
			Twenty_Node_Brick_Inverse(j,i) = 1/4*(1-(Twenty_Brick_Coordinates(j,1)*SQRT_3_5)^2)* ...
												(1+Twenty_Brick_Coordinates(j,2)*Twenty_Brick_Coordinates(i,2)*SQRT_3_5)* ...
												(1+Twenty_Brick_Coordinates(j,3)*Twenty_Brick_Coordinates(i,3)*SQRT_3_5);
		end
		for i=[10,12,14,16]
			Twenty_Node_Brick_Inverse(j,i) = 1/4*(1-(Twenty_Brick_Coordinates(j,2)*SQRT_3_5)^2)* ...
												(1+Twenty_Brick_Coordinates(j,1)*Twenty_Brick_Coordinates(i,1)*SQRT_3_5)* ...
												(1+Twenty_Brick_Coordinates(j,3)*Twenty_Brick_Coordinates(i,3)*SQRT_3_5);
		end
		for i=[17,18,19,20]
			Twenty_Node_Brick_Inverse(j,i) = 1/4*(1-(Twenty_Brick_Coordinates(j,3)*SQRT_3_5)^2)* ...
												(1+Twenty_Brick_Coordinates(j,1)*Twenty_Brick_Coordinates(i,1)*SQRT_3_5)* ...
												(1+Twenty_Brick_Coordinates(j,2)*Twenty_Brick_Coordinates(i,2)*SQRT_3_5);
		end

	end

	Twenty_Node_Brick_Inverse

	Twenty_Node_Brick_Inverse = inv(Twenty_Node_Brick_Inverse);

end