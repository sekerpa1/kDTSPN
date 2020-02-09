include("./linked_stack.jl");
include("./finding_cliques.jl");
using .Stack	
using Test
using PyPlot

Circle = matplotlib.patches.Circle;

struct Point{T<:AbstractFloat}
	x::T
	y::T
end

struct Region{T}
	center::Point{T}
	radius::T
end

struct IntersectionPoint{T<:AbstractFloat}
	x::T
	y::T
	parents::Vector{Region{T}}
end

function isequal(a::Region, b::Region)
	return a.center.x == b.center.x && a.center.y == b.center.y && a.radius == b.radius;
end

function mutual_parents(parents_1::Vector{Region}, parents_2::Vector{Region})
	
	parents = [];
	for i in 1:length(parents_1)
		for j in 1:length(parents_2)
			if isequal(parents_1[i], parents_2[j])
				push!(parents);
			end
		end
	end
	return parents;
end

#= takes as input list of element arrays
 and groups them according to length
 for [[1,2,3],[2,3],[3,4],[1]] -> [[[1,2,3]],[[2,3],[3,4]],[[1]]]
=#
function group_list(list)
	# remove empty elements
	
	list = filter(x -> !isempty(x), list);
	# sort elements by length
	list_sorted = sort(list, by = x -> length(x), rev=true);
	
	# group them by length
	grouped = [[] for i in 1:length(list_sorted[1])];
	println(length(grouped));
	map(x -> push!(grouped[length(x)], x) ,list_sorted);
	flattened = [];
	for i in 1:length(grouped)
		if length(grouped) != 0 || isequal(grouped[i][1],[]) 
			push!(flattened, grouped[i]);
		end
	end
	
	flattened = filter(x -> !isempty(x), flattened);
	sort!(flattened, by = x -> length(x[1]), rev=true);
	return flattened;
end

# remove all elements of list that contain used elements
function remove_used(list, used_elements)
	new_list = [];
	for i in 1:length(list)
		inner_list = [];
		map(x -> !contains_any_element(x, used_elements) ? push!(inner_list, x) : return, list[i]);
		if !isempty(inner_list)
			push!(new_list, inner_list);
		end
	end
	return new_list;
end

function point_lies_within_all_circles(point::Point, regions::Vector{Region})
	return sum(x -> point_lies_within_circle(point, x), regions) == length(regions);
end

function point_lies_within_circle(point::Point, region::Region)
	return distance(point, region.center) < region.radius;
end

function center_of_points(points::Vector{Point})
	
	center = Point(0.0,0.0);
	for i in 1:length(points)
		center.x += points[i].x;
		center.y += points[i].y;
	end
	
	return center/length(points);
end

function calculate_intersection_area(points::Vector{IntersectionPoint}, circles::Vector{Region})
	
	# filter out points that do not lie in all circles
	points_inside = filter(x -> point_lies_within_all_circles(x, circles), points);
	
	# find center of points
	center = center_of_points(points_inside);
	
	#find angle of points and center
	angles = (p -> atan(p.x - center.x, p.y - center.y), points_inside);
	
	
	# calculate radius of circle centered in the center of polygon area
	radius = sqrt((points_inside[1].x - center.x)^2 + (points_inside[2].y - center.y)^2);
	
	# sort polygon points on sizes of angles
	polygon_points = [(angles[i], points_inside[i]) for i in 1:length(angles)];
	sort!(polygon_points);
	
	# calculate polygon area
	area = 0;
	for i in 1:length(polygon_points)
	
		# update polygon area
		alpha = i == 1 ? angles[1] : angles[i] - angles[i-1];
		area += 1/2*r^2*sin(alpha);
	end
	
	p1 = points_inside[length(points_inside)];
	
	# calculate arc area
	for i in 1:length(points_inside)
		p2 = points_inside[i];
		
		# get mutual parents of points
		parents = mutual_parents(p1.parents, p2.parents);
		arc_length = typemax(Float64);
		index = 0;
		arc_angle = 0;
		
		for j in 1:length(parents)
			center = parents[i].center;
			radius = parents[i].radius;
			
			# calculate arc length
			
			#find angle of points and center
			alpha_p1 = atan(p1.x - center.x, p1.y - center.y);
			alpha_p2 = atan(p2.x - center.x, p2.y - center.y);
			angle_between = abs(alpha_p1 - alpha_p2);
			
			length = angle_between/2π * (2*π*radius);
			
			# keep minimum length arc, region
			if(length < arc_length)
				#update parent index and arc length
				index = j;
				arc_length = length;
				arc_angle = angle_between;
			end
		end
		
		# calculate arc area with minimum length
		
		# calculate area of the part of the circle
		circle_part = π*r^2 * (arc_angle/(2*π));
		
		#calculate area in between center and arc
		area_center_arc = 1/2*r^2*sin(arc_angle);
		
		arc_area = circle_part - area_center_arc;
		# add up to the whole area
		area += arc_area;
		
		p1 = p2;
	end
	
	return (center, area);
end




#= find intersection points of two regions. If two circles intersect,
 the algorithm first projects circles, so that the left circle is placed at 
 point (0,0) and both centers intersect x axis. It calculates the intersection
 points in this space and then projects them back according to original circles positions
=#
function intersection(a::Region, b::Region)
	
	center_distance = distance(a.center.x, a.center.y, b.center.x, b.center.y);
	if center_distance > a.radius + b.radius
		return IntersectionPoint();
	elseif center_distance == a.radius + b.radius
		# if circles are touching each other
		mag = a.radius/center_distance;
		return IntersectionPoint(a.center.x + (a.center.x - b.center.x)*mag, a.center.y + (a.center.y - b.center.y)*mag,
			[a,b]);
	elseif center_distance + a.radius == b.radius || center_distance  + b.radius == a.radius
		# if one circle is inside another circle
		(smaller, bigger) = a.radius > b.radius ? (b,a) : (a,b);
		
		# calculate vector between centers
		v = (smaller.center.x - bigger.center.x, smaller.center.y - bigger.center.y);
		
		# normalize vector
		v = v/sqrt(v[1]^2, v[2]^2);
		
		# calculate intersection point using radius of bigger circle as magnitude 
		v = v .* bigger.radius;
		return IntersectionPoint(v[1], v[2], [a,b]);
	end
	
	# alias the left most circle
	if a.center.x < b.center.x
		left = a;
	else
		left = b;
	end
	right = isequal(left,a) ? b : a;
	
	# move centers by x axis
	a_center_x = 0;
	b_center_x = right.center.x - left.center.x;
	
	# move centers by y axis
	displacement = min(abs(right.center.y), abs(left.center.y));
	a_center_y = a.center.y - displacement;
	b_center_y = b.center.y - displacement;
	
	# get angle between right point and x axis
	cos_alpha = b_center_x / sqrt(b_center_x^2 + b_center_y^2);
	alpha = acos(cos_alpha);
	
	# get angle direction of projection of centers a and b
	direction_positive = b_center_y > 0;
	
	# rotate right center to intersect x axis
	b_center_x = sqrt(b_center_x^2+b_center_y^2);
	b_center_y = 0;
	
	# find intersections
	x = (left.radius^2 - right.radius^2 + b_center_x^2) / (2 * b_center_x);
	y1 = sqrt(left.radius^2 - x^2);
	y2 = -y1;
	
	# project points back
	m = sqrt(x^2 + y1^2);
	
	# get angle between upper point and x axis
	cos_beta = x/m;
	beta = acos(cos_beta);
	
	# project point x back
	
	# project upper intersection point
	if direction_positive
		gamma_1 = alpha + beta;
	else
		gamma_1 = beta - alpha;
	end
	
	x1, y1 = cos(gamma_1)*m, sin(gamma_1)*m;
	
	# project lower intersection point
	if direction_positive
		gamma_2 = alpha - beta;
	else
		gamma_2 = -beta - alpha;
	end
	
	x2, y2 = cos(gamma_2)*m, sin(gamma_2)*m;
	
	# move points back by x and y axis
	x1,x2 = x1 + left.center.x, x2 + left.center.x;
	y1,y2 = y1 + displacement, y2 + displacement;
	
	return (IntersectionPoint(x1,y1, [a,b]), IntersectionPoint(x2,y2, [a,b]));
end

function distance(x1, y1, x2, y2)
	return √((x1 - x2) ^ 2
		+ (y1 - y2) ^ 2);
end

function distance(a::Point, b::Point)
	return distance(a.x, a.y, b.x, b.y);
end

function project_points_on_center(x1,y1,x2,y2)
	minimum = min(x1,x2);
	x1 -= minimum;
	x2 -= minimum;
end

#=
	Returns true if list contains any of the used elements
=#
function contains_any_element(list, used_elements)
	return sum(map(y -> !isempty(filter(x -> x == y, list)), used_elements)) > 0;
end

function testIntersectionArea()
	a = Region(Point(0.02,0.77), 0.75);
	b = Region(Point(0.89,0.08), 0.78);
	c = Region(Point(0.74,0.91), 0.95);
	
	intersection_points = Vector{IntersectionPoint}();
	intersection_p  = intersection(a,b);
	push!(intersection_points, intersection_p[1]);
	push!(intersection_points, intersection_p[2]);
	
	intersection_p  = intersection(b,c);
	push!(intersection_points, intersection_p[1]);
	push!(intersection_points, intersection_p[2]);
	
	intersection_p  = intersection(a,c);
	push!(intersection_points, intersection_p[1]);
	push!(intersection_points, intersection_p[2]);
	println(typeof(intersection_points);
	println(typeof([a,b,c]));
	
	(center, area) = calculate_intersection_area(intersection_points, [a,b,c]);
	println(area);
	#Area:	0.28245
	#Polygon Area:	0.09531
	#Arc Area:	0.18714
end

testIntersectionArea();

function testRegionIntersections()
	#a = Region(Point(3.0,4.0),2.0);
	#b = Region(Point(5.0,6.5),3.0);
	a = Region(Point(1.0,2.0), 3.5);
	b = Region(Point(4.0,4.0), 3.0);
	
	pygui(true);
	plt.figure(1);
    plt.clf();
	ax = plt.gca();
	
	circle = Circle([a.center.x,a.center.y], a.radius, edgecolor= "red", 
		facecolor="red", alpha=0.2);
	ax.add_patch(circle);
	circle = Circle([b.center.x, b.center.y], b.radius, edgecolor= "red", 
		facecolor="red", alpha=0.2);
	ax.add_patch(circle);
	
	points = intersection(a,b);
	
	plt.scatter(points[1][1],points[1][2], color="black");
	plt.scatter(points[2][1],points[2][2], color="black");
	
	plt.title("Length = ");
    plt.axis("equal");
    plt.tight_layout();
    plt.pause(0.1);
end

#testRegionIntersections();
#println(contains_any_element([1,2,3,4,5] , [2,3]));

#input = [[1,2,3],[2,3],[3,4],[1],[3,5,6,7], [], [8,7,2,5,7,8],[4,5]];
#grouped = group_list(input);
#println(grouped);
#println("");
#unused =remove_used(grouped, [1,2]);
#println(unused);

