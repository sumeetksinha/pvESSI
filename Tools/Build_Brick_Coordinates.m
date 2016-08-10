
function  Brick_Coordinates = Build_Brick_Coordinates()

	Brick_Coordinates = zeros(27,3);

	Brick_Coordinates(1,:)=[1,1,1]; 
	Brick_Coordinates(2,:)=[-1,1,1]; 
	Brick_Coordinates(3,:)=[-1,-1,1]; 
	Brick_Coordinates(4,:)=[1,-1,1]; 
	Brick_Coordinates(5,:)=[1,1,-1]; 
	Brick_Coordinates(6,:)=[-1,1,-1]; 
	Brick_Coordinates(7,:)=[-1,-1,-1]; 
	Brick_Coordinates(8,:)=[1,-1,-1]; 
	Brick_Coordinates(9,:)=[0,1,1]; 
	Brick_Coordinates(10,:)=[-1,0,1];
	Brick_Coordinates(11,:)=[0,-1,1];
	Brick_Coordinates(12,:)=[1,0,1];
	Brick_Coordinates(13,:)=[0,1,-1];
	Brick_Coordinates(14,:)=[-1,0,-1];
	Brick_Coordinates(15,:)=[0,-1,-1];
	Brick_Coordinates(16,:)=[1,0,-1];
	Brick_Coordinates(17,:)=[1,1,0];
	Brick_Coordinates(18,:)=[-1,1,0];
	Brick_Coordinates(19,:)=[-1,-1,0];
	Brick_Coordinates(20,:)=[1,-1,0];
	Brick_Coordinates(21,:)=[0,0,0];
	Brick_Coordinates(22,:)=[0,1,0];
	Brick_Coordinates(23,:)=[-1,0,0];
	Brick_Coordinates(24,:)=[0,-1,0];
	Brick_Coordinates(25,:)=[1,0,0];
	Brick_Coordinates(26,:)=[0,0,1];
	Brick_Coordinates(27,:)=[0,0,-1];

end