function Eight_Node_Brick_Inverse = Eight_Node_Brick_Inverse()

	Eight_Brick_Coordinates = Build_Brick_Coordinates();
	Eight_Brick_Coordinates = Eight_Brick_Coordinates(1:8,:);

	Eight_Node_Brick_Inverse = zeros(8,8);

	SQRT_3 = sqrt(1/3);


	for j=1:8
		for i=1:8
			Eight_Node_Brick_Inverse(j,i) = 1/8*(1+Eight_Brick_Coordinates(j,1)*Eight_Brick_Coordinates(i,1)*SQRT_3)* ...
												(1+Eight_Brick_Coordinates(j,2)*Eight_Brick_Coordinates(i,2)*SQRT_3)* ...
												(1+Eight_Brick_Coordinates(j,3)*Eight_Brick_Coordinates(i,3)*SQRT_3);
		end
	end

	Eight_Node_Brick_Inverse

	Eight_Node_Brick_Inverse = inv(Eight_Node_Brick_Inverse);

end