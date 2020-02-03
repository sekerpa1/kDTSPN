

# finds roots of quadratic equation a*x^2 + b*x + c = 0
function solve_quadratic_equation(a, b, c)
	discriminant = b^2 - 4*a*c;
	if(discriminant < 0)
		error("No solution which lies within set of real numbers found");
	end
	return ((-b + sqrt(discriminant)) / (2*a), (-b - sqrt(discriminant)) / (2*a));
end


#=
 returns two points (x1,y1) and (x2,y2) that lie on the intersection of circle with line
 a,b,c determine the equation of the line and (s1,s2) - center of the circle, r determines
 the radius of the circle
=#
function find_intersection_of_line_and_circle(a1, b1, c1, s1, s2, r)
	if(a1 == 0)
		y = -c1/b1;
		b = -2*s1;
		c = s1^2 + y^2 - 2*s2*y + s2^2 - r^2;
		(x1, x2) = solve_quadratic_equation(1, b, c);
		return [(x1,y),(x2,y)];
	end
	a = b1^2/a1^2 + 1;
	b = 2*b1*c1/a1^2 + 2*s1*b1/a1 - 2*s2;
	c = c1^2/a1^2+2*s1*c1/a1 + s1^2 + s2^2 - r^2;
	
	(y1, y2) = solve_quadratic_equation(a, b, c);
	
	if y1 == y2
		b = -2*s1;
		c = s1^2 + y1^2 - 2*s2*y1 + s2^2 - r^2;
		(x1, x2) = solve_quadratic_equation(1, b, c);
		return [(x1,y1),(x2,y1)];
	else
		x1 = (-b1*y1 - c1) / a1;
		x2 = (-b1*y2 - c1) / a1;
		return [(x1,y1),(x2,y2)];
	end
end


function move_point_closer_to_circle(x, y, s1, s2,r)
	while !point_lies_within_circle(x,y,s1,s2,r)
		v = (s1 - x, s2 - y);
		x = x + 0.001*v[1];
		y = y + 0.001*v[2];
	end
	
	return (x,y);
end

#=
  Find projection for point point_c. Where point_a, point_b are neighbourhood points of point point_c.
  point_c is the projected point and region_c is the region of point point_c.
=#
#=function projected_point(point_a::Point{T}, point_b::Point{T}, point_c::Point{T}, region_c::Region{T}) where T <: AbstractFloat
	
	#println("Finding projected point for : [ $(point_a.x) , $(point_a.y)] , [ $(point_b.x), $(point_b.y)] and point c : [ $(point_c.x), $(point_c.y)],
	# and region: [$(region_c.center.x) , $(region_c.center.y) , $(region_c.radius)]");
	
	
	# for points point_a, point_b find general equation of line
	normal_vector_ab = (point_b.y - point_a.y, point_a.x - point_b.x);
	general_line_equation_ab = (normal_vector_ab[1], normal_vector_ab[2], -normal_vector_ab[1]*point_a.x -normal_vector_ab[2]*point_a.y);
	#println(normal_vector_ab);
	#println(general_line_equation_ab);
	
	
	# for line between points point_a, point_b we find general equation for perpendicular line
	# to this line intersecting point_c
	normal_vector_c_ab = (-normal_vector_ab[2], normal_vector_ab[1]);
	general_line_equation_c_ab = (normal_vector_c_ab[1], normal_vector_c_ab[2], -normal_vector_c_ab[1]*point_c.x 
	-normal_vector_c_ab[2]*point_c.y);
	#println(normal_vector_c_ab);
	#println(general_line_equation_c_ab);
	
	#find intersection of line between points point_a and point_b and perpendicular line
	# to this line containing point_c
	(intersection_x, intersection_y) = find_intersection_of_line_and_line(general_line_equation_ab[1], general_line_equation_ab[2],general_line_equation_ab[3],
		general_line_equation_c_ab[1],general_line_equation_c_ab[2],general_line_equation_c_ab[3]);
	#println("Intersection",intersection_x, " ", intersection_y);
	
	
	# find out if intersection lies on the line between points point_a and point_b
	is_between = max(distance(intersection_x, intersection_y, point_a.x, point_a.y),
		distance(intersection_x, intersection_y, point_b.x, point_b.y)) <=  distance(point_a.x, point_a.y, point_b.x, point_b.y);
	if is_between
		# find intersection of perpendicular line with region
		#println(is_between);
		points = find_intersection_of_line_and_circle(general_line_equation_c_ab[1],general_line_equation_c_ab[2],general_line_equation_c_ab[3],
			region_c.center.x, region_c.center.y, region_c.radius);
		#println(points);
		# return point closer to intersection point of two lines
		
		move_point_closer_to_circle!(points[1][1], points[1][2], region_c.center.x, region_c.center.y, region_c.radius);
		move_point_closer_to_circle!(points[2][1], points[2][2], region_c.center.x, region_c.center.y, region_c.radius);
		
		if(!point_lies_within_circle(points[1][1], points[1][2], region_c.center.x, region_c.center.y, region_c.radius))
			println(points, region_c);
			error("point outside circle");
		end
		
		if(!point_lies_within_circle(points[2][1], points[2][2], region_c.center.x, region_c.center.y, region_c.radius))
			println(points, region_c);
			error("point outside circle");
		end
		
		if distance(points[1][1], points[1][2],intersection_x, intersection_y) < distance(points[2][1], points[2][2],intersection_x, intersection_y)
			return Point(points[1]);
		else
			return Point(points[2]);
		end
	else
		# find center of points point_a and point_b
		(center_x, center_y) = ((point_b.x + point_a.x)/2, (point_b.y + point_a.y)/2);
		
		# calculate line equation for line connecting point_c and center
		normal_vector = (point_c.y - center_y, center_x - point_c.x);
		general_line_equation = (normal_vector[1], normal_vector[2], -normal_vector[1]*point_c.x -normal_vector[2]*point_c.y);
		
		# find intersection of line and region
		points = find_intersection_of_line_and_circle(general_line_equation[1],general_line_equation[2],general_line_equation[3],
			region_c.center.x, region_c.center.y, region_c.radius);
		
		move_point_closer_to_circle!(points[1][1], points[1][2], region_c.center.x, region_c.center.y, region_c.radius);
		move_point_closer_to_circle!(points[2][1], points[2][2], region_c.center.x, region_c.center.y, region_c.radius);
		
		if(!point_lies_within_circle(points[1][1], points[1][2], region_c.center.x, region_c.center.y, region_c.radius))
			println(points, region_c);
			error("point outside circle");
		end
		
		if(!point_lies_within_circle(points[2][1], points[2][2], region_c.center.x, region_c.center.y, region_c.radius))
			println(points, region_c);
			error("point outside circle");
		end
		
		
		# return point closer to center of points point_a and point_b
		if distance(points[1][1], points[1][2],center_x, center_y) < distance(points[2][1], points[2][2],center_x, center_y)
			return Point(points[1]);
		else
			return Point(points[2]);
		end
	end
end=#

#=
  takes adjacency matrix of regions as input and finds all the overlapping
  regions from largest number of overlapping elements to smallest number
  of overlapping elements
=#
function squash_regions(A)
	# find cliques in the given graph
	C = populate_cliques(A);
	
	# order cliques from largest to smallest
	ordered = [];
	for i in 1:length(C)
		push!(ordered, (C[i], length(C[i])));
	end
	sort!(ordered, by = x -> x[2]);
	
	# separate cliques into arrays according to sizes
	separated = [filter(x -> x[2] == i, ordered) for i in 1:size(A,1)];
	println(separated);
	
	
end