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

struct Region{T<:AbstractFloat}
	center::Point{T}
	radius::T
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

#= find intersection points of two regions. If two circles intersect,
 the algorithm first projects circles, so that the left circle is placed at 
 point (0,0) and both centers intersect x axis. It calculates the intersection
 points in this space and then projects them back according to original circles positions
=#
function intersection(a::Region, b::Region)
	
	center_distance = distance(a.center.x, a.center.y, b.center.x, b.center.y);
	if center_distance > a.radius + b.radius
		return ();
	elseif center_distance == a.radius + b.radius
		# if circles are touching each other
		mag = a.radius/center_distance;
		return (a.center.x + (a.center.x - b.center.x)*mag, a.center.y + (a.center.y - b.center.y)*mag);
	elseif center_distance + a.radius == b.radius || center_distance  + b.radius == a.radius
		# if one circle is inside another circle
		(smaller, bigger) = a.radius > b.radius ? (b,a) : (a,b);
		
		# calculate vector between centers
		v = (smaller.center.x - bigger.center.x, smaller.center.y - bigger.center.y);
		
		# normalize vector
		v = v/sqrt(v[1]^2, v[2]^2);
		
		# calculate intersection point using radius of bigger circle as magnitude 
		return v .* bigger.radius;
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
	
	return ((x1,y1),(x2,y2));
end

function distance(x1, y1, x2, y2)
	return √((x1 - x2) ^ 2
		+ (y1 - y2) ^ 2);
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

