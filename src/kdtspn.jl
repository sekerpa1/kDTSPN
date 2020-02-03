#!/usr/bin/env julia

module kDTSPN

import Base: show

using Dubins
using Random
using Test
using PyPlot

include("./linked_stack.jl")
include("./finding_cliques.jl")
using .Stack

Circle = matplotlib.patches.Circle

input_file_path = "./../input/problem.txt";
delimiter = " ";
population_size = 200;
crossover_ratio = 0.2;

mutation_probability = 0.5;
number_of_iterations = 200;
turning_radius = 0.7;

DEFAULT_FITTNESS = typemax(Float64);

# local search stops when optimized points converge to a value that is
# less than this number away from previously optimized value
local_search_stop_criterion = 0.001;


struct Point{T<:AbstractFloat}
	x::T
	y::T
end

Point(p::Tuple{T,T}) where T <: AbstractFloat = Point{T}(p[1], p[2]);

struct Region{T<:AbstractFloat}
	center::Point{T}
	radius::T
end

function show(r::Region, i)
	println("Region $i: Center [$(r.center.x) $(r.center.y)] Radius $(r.radius)");
end

mutable struct Individual
	permutation::Vector{Int32}
	starting_points::Vector{Int32}
	starting_points_positions::Vector{Point{Float64}}
	starting_points_directions::Vector{Float64}
	points_positions::Vector{Point{Float64}}
	points_directions::Vector{Float64}
	fittness::Float64
end

isless(a::Main.kDTSPN.Individual, b::Main.kDTSPN.Individual) = isless(a.fittness, b.fittness);

#=
  Takes input and produces set of regions to visit
=#
function parseInput()
	regions = Vector{Region}();
	number_vehicles = 0;
	starting_region = 0;
	i = 1;
	for l in eachline(input_file_path)
		item = split(l, delimiter);
		if i == 1
			if(length(item) > 1)
				error("Wrong input, provide number of vehicles in the first line.");
			end
			number_vehicles = parse(Int32, l);
		elseif i == 2
			# parse starting region
			if(length(item) != 4)
				error("Wrong input, error at line $i");
			end
			starting_region = Region(Point(parse(Float64, convert(String, item[2])),parse(Float64, convert(String, item[3]))), parse(Float64, convert(String, item[4])));
		else
			if(length(item) != 4)
				error("Wrong input, error at line $i");
			end
			r = Region(Point(parse(Float64, convert(String, item[2])),parse(Float64, convert(String, item[3]))), parse(Float64, convert(String, item[4])));
			push!(regions, r);
		end
		i = i + 1;
    end
	
	return (number_vehicles, starting_region, regions);
end

#=
  for given array of regions, it returns adjacency matrix in form 
  A[i,j] == 1 - if region i overlaps with region j
  A[i,j] == 0 otherwise
=#
function overlaps(a)
	A = [(i == j) ? 1 : 0 for i in 1:length(a), j in 1:length(a)];
	for i in 1:length(a), j in i+1:length(a)
		A[i,j] = overlaps(a[i],a[j]);
	end
	
	return A;
end

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

#=
  Calculate Dubins shortest path for individual using Dubins maneuvers
=#
function evaluate_individual!(ind::Individual, regions::Vector{Region})

	max_path_length = 0;
	for i in 1:length(ind.starting_points)
	
		last = i == length(ind.starting_points) ? length(ind.permutation) : ind.starting_points[i+1]-1;
		current_points = vcat(ind.starting_points_positions[i], ind.points_positions[ind.permutation[ind.starting_points[i]:last]]);
		directions = vcat(ind.starting_points_directions[i], ind.points_directions[ind.permutation[ind.starting_points[i]:last]]);
		
		d_path_vehicle = [];
		for j in 1:length(current_points)-1
			d_path = dubins_shortest_path([current_points[j].x, current_points[j].y, directions[j]], 
				[current_points[j+1].x,current_points[j+1].y, directions[j+1]], turning_radius);
			if d_path[1] != 0
				error("Error when calculating dubins path, error code $(d_path[1])");
			end
			push!(d_path_vehicle, d_path[2]);
		end
		
		d_path = dubins_shortest_path([current_points[length(current_points)].x, current_points[length(current_points)].y, directions[length(current_points)]], 
				[current_points[1].x,current_points[1].y, directions[1]], turning_radius);
		if(d_path[1] != 0)
			error("Error when calculating dubins path, error code $(d_path[1])");
		end
		
		max_path_length = max(max_path_length, sum(broadcast(x -> dubins_path_length(x), d_path_vehicle)));
	end
	
	ind.fittness = max_path_length;
end

#=
  calculate direction for dubins path for point c and its neighbors point a and point b
  point a lies immediately before point c in sequence and point b lies immediately after point c
=#
function get_direction(a::Point{T}, b::Point{T}, c::Point{T}) where T <: AbstractFloat
	v = (b.x - a.x, b.y - a.y);
	n = sqrt(v[1]^2 + v[2]^2);
	v = v ./ n;
	cos_l1 = v[1]*1/(sqrt(v[1]^2+v[2]^2)*sqrt(1));
	cos_l2 = v[2]*1/(sqrt(v[1]^2+v[2]^2)*sqrt(1));
	
	if acos(cos_l2) > pi/2
		return 2*pi - acos(cos_l1);
	else
		return acos(cos_l1);
	end
end

#=
  calculates new positions for points of individual, accepts individual definition
  and vector of regions, which contains union of starting region with rest of the regions
  regions = starting_region ∪ regions. Updates the supplied individual with new set of points
=#
function local_search!(ind::Individual, regions::Vector{Region})
	
	points = [ind.permutation[i] + 1 for i in 1:length(ind.permutation)];
	
	#add starting points to sequence of points
	all_points = [];
	for i in 1:length(ind.starting_points)-1
		all_points = vcat(all_points, 1, points[ind.starting_points[i]:ind.starting_points[i+1]-1]);
	end
	
	all_points = vcat(all_points, 1, points[ind.starting_points[length(ind.starting_points)]:end]);
	
	# update starting points indexes 
	starting_points = [ind.starting_points[i]+i-1 for i in 1:length(ind.starting_points)];
	
	#reinitialize point positions to center of circles
	positions = [regions[i].center for i in 2:length(regions)];
	ind.points_positions[1:end] = positions[1:end];
	
	for i in 1:length(ind.starting_points_positions)
		ind.starting_points_positions[i] = regions[1].center;
	end
	
	len = length(starting_points);
	for i in 1:len
		diff = typemax(Float64);
		current_sequence = (i==len ? all_points[starting_points[i]:end] : all_points[starting_points[i]:starting_points[i+1]-1]);
		
		# add starting point specific to sequence and rest of the points to current points
		current_sequence_to_ind_points_positions = [x - 1 for x in current_sequence[2:end]];
		current_points = vcat(ind.starting_points_positions[i], ind.points_positions[current_sequence_to_ind_points_positions]);
		
		# while not converged optimize current subsequence
		while diff > local_search_stop_criterion
			
			# optimize points
			optimized = [j == length(current_sequence) ? 
			projected_point(current_points[j-1], current_points[1], 
				current_points[j], regions[current_sequence[j]]) :
			j == 1 ? projected_point(current_points[length(current_sequence)], current_points[j+1], 
				current_points[j], regions[current_sequence[j]]) :
			projected_point(current_points[j-1],current_points[j+1], 
				current_points[j], regions[current_sequence[j]])
			for j in 1:length(current_sequence)];
			
			diff = distance(optimized, current_points);
			current_points = optimized;
		end
		
		# copy optimized points back to individual
		ind.starting_points_positions[i] = current_points[1];
		
		# start from second point as first is starting point
		for j in 2:length(current_sequence)
			ind.points_positions[current_sequence[j]-1] = current_points[j];
		end
		
		#=
		# optimize point directions
		current_directions = [j == 1 ? get_direction(current_points[end],current_points[j+1],current_points[j]) :
			j == length(current_sequence) ? get_direction(current_points[j-1],current_points[1],current_points[j]) :
			get_direction(current_points[j-1], current_points[j+1], current_points[j])
			for j in 1:length(current_sequence)];
		
		# copy directions back to individual
		ind.starting_points_directions[i] = current_directions[1];
		
		# start from second point as first is starting point
		for j in 2:length(current_sequence)
			ind.points_directions[current_sequence[j]-1] = current_directions[j];
		end =#
	end
end


function projected_point(point_a::Point{T}, point_b::Point{T}, point_c::Point{T}, region_c::Region{T}) where T <: AbstractFloat
	
	# for points point_a, point_b find general equation of line
	normal_vector_ab = (point_b.y - point_a.y, point_a.x - point_b.x);
	general_line_equation_ab = (normal_vector_ab[1], normal_vector_ab[2], -normal_vector_ab[1]*point_a.x -normal_vector_ab[2]*point_a.y);
	
	# for line between points point_a, point_b we find general equation for perpendicular line
	# to this line intersecting point_c
	normal_vector_c_ab = (-normal_vector_ab[2], normal_vector_ab[1]);
	general_line_equation_c_ab = (normal_vector_c_ab[1], normal_vector_c_ab[2], -normal_vector_c_ab[1]*point_c.x 
	-normal_vector_c_ab[2]*point_c.y);

	#find intersection of line between points point_a and point_b and perpendicular line
	# to this line containing point_c
	(intersection_x, intersection_y) = find_intersection_of_line_and_line(general_line_equation_ab[1], general_line_equation_ab[2],general_line_equation_ab[3],
		general_line_equation_c_ab[1],general_line_equation_c_ab[2],general_line_equation_c_ab[3]);
	
	# find out if intersection lies on the line between points point_a and point_b
	is_between = max(distance(intersection_x, intersection_y, point_a.x, point_a.y),
		distance(intersection_x, intersection_y, point_b.x, point_b.y)) <=  distance(point_a.x, point_a.y, point_b.x, point_b.y);
	
	projection_vector = normal_vector_c_ab;
	if !is_between
		#= if intersection point does not lie between the points point_a and point_b, calculate projection vector
		 using the middle point of point_a and point_b =#
		(center_x, center_y) = ((point_b.x + point_a.x)/2, (point_b.y + point_a.y)/2);
		projection_vector = (center_x - point_c.x, center_y - point_c.y);
	end
	
	#normalize the projection vector
	val = sqrt(projection_vector[1]^2+projection_vector[2]^2);
	projection_vector = projection_vector ./ val;
	
	# optimize point position using projection vector
	stopping_diff = 0.1;
	n = 1;
	initial_point = (point_c.x, point_c.y);
	halved_times = 0;
	while(true)
		projected_point = initial_point .+ projection_vector .* n;
		if !point_lies_within_circle(projected_point[1], projected_point[2], region_c.center.x, region_c.center.y, region_c.radius)
			# half the projection magnitude maximum five times
			if halved_times > 5
				return Point(initial_point[1], initial_point[2]);
			end
			n = n/2;
			halved_times = halved_times + 1;
		else
			if distance(initial_point[1], initial_point[2], projected_point[1], projected_point[2]) < stopping_diff
				return Point(projected_point[1], projected_point[2]);
			else
				initial_point = projected_point;
				halved_times = 0;
				n = n*2;
			end
		end
	end
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

function move_point_closer_to_circle(x, y, s1, s2,r)
	while !point_lies_within_circle(x,y,s1,s2,r)
		v = (s1 - x, s2 - y);
		x = x + 0.001*v[1];
		y = y + 0.001*v[2];
	end
	
	return (x,y);
end

function point_lies_within_circle(x, y, s1, s2, r)
	return (x-s1)^2 + (y-s2)^2 <= r^2;
end

#=
  takes general equation of line 1 represented by parameters a1, b1, c1
  and general equation of line 2 represented by parameters a2, b2, c2
  and returns intersection point 
=#
function find_intersection_of_line_and_line(a1,b1,c1,a2,b2,c2)
  if a1 == 0
	if a2 == 0
	  return (a1,b1,c1);
	end
    y = c1/b1;
	x = -c2/a2 - (b2*c1)/(a2*b1);
	return (x,y);
  end
  
  y = (-c2 + c1*a2/a1)/(-b1*a2/a1 +b2);
  x = (-b1*y - c1)/a1;
  
  return (x,y);
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

function distance(a::Vector{Point{Float64}},b::Vector{Point{Float64}})
	return sum([distance(a[i].x,a[i].y,b[i].x,b[i].y) for i in 1:length(a)])
end

function distance(x1, y1, x2, y2)
	return √((x1 - x2) ^ 2
		+ (y1 - y2) ^ 2);
end

function distance(a::Point{Float64}, b::Point{Float64})
	return √((a.x - b.x) ^ 2
		+ (a.y - b.y) ^ 2);
end

# finds roots of quadratic equation a*x^2 + b*x + c = 0
function solve_quadratic_equation(a, b, c)
	discriminant = b^2 - 4*a*c;
	if(discriminant < 0)
		error("No solution which lies within set of real numbers found");
	end
	return ((-b + sqrt(discriminant)) / (2*a), (-b - sqrt(discriminant)) / (2*a));
end


function overlaps(a::Region, b::Region)

	center_distance = √((a.center.x - b.center.x) ^ 2
		+ (a.center.y - b.center.y) ^ 2);

	return center_distance - a.radius - b.radius < 0
end

#=
  Generates random individual based on the input. Accepts total_points - total points count without starting point.
  starting_points_count - number of starting points used, is equal to number of vehicles used
  starting_point - coordinates of starting point
  points - coordinates of individual points without starting point
=#
function generate_individual(total_points_count::Integer, starting_points_count::Integer, starting_point::Point{Float64}, points::Vector{Point{Float64}})

	if starting_points_count < 2
		throw(ArgumentError("individual should have at least two starting points"))
	end

	if total_points_count < 2
		throw(ArgumentError("Individual should have at least two total points"))
	end


	permutation = randperm(total_points_count);

	taken = [];
	
	# automatically add 1 as starting point
	push!(taken, 1);
	
	# add other starting points randomly, starting from 2...
	free = LinkedStack{Integer}();
	free_size = total_points_count-1;
	for i in 2:total_points_count
		push!(free, i);
	end
	
	for i in 2:starting_points_count
		pos = rand(1:free_size);
		push!(taken, popAt!(free, pos));
		free_size = free_size - 1;
	end

	sort!(taken);

	starting_points = [starting_point for i in 1:starting_points_count];
	
	points_directions = [pi for i in 1:total_points_count];
	
	starting_points_directions = [pi for i in 1:starting_points_count];
	
	return Individual(permutation, taken, starting_points, starting_points_directions, points[permutation], points_directions, DEFAULT_FITTNESS);
end

function crossover(left_ind::Individual, right_ind::Individual, point::Integer, keep_number_of_starting_points)

	if point >= length(left_ind.permutation)-1 || point <= 1
		error("Crossover point outside possible range");
	end

	# copy the left part of the permutation from left individual
	child_perm = copy(left_ind.permutation[1:point]);

	# copy the right leftover part of the permutation from left individual
	leftover_perm = LinkedStack{Integer}();
	for i in point+1:length(left_ind.permutation)
		push!(leftover_perm, left_ind.permutation[i]);
	end

	# copy the right part of the permutation from right individual
	duplicates = [];
	for i in point+1:length(left_ind.permutation)
		isavailable = pop!(leftover_perm, right_ind.permutation[i]);
		push!(child_perm, right_ind.permutation[i]);

		# if element has been added already, store the index of the element
		if !isavailable
			push!(duplicates, i);
		end
	end

	# assign the rematining points that are not part of the permutation
	for i in 1:length(duplicates)
		child_perm[duplicates[i]] = pop!(leftover_perm);
	end

	# copy the left part of the starting points from left individual
	left_starting_points = filter((x) -> x <= point, left_ind.starting_points);

	# copy the right part of the starting points from right individual
	right_starting_points = filter((x) -> x > point, right_ind.starting_points);

	total_points = length(left_starting_points) + length(right_starting_points);
	#= calculate the difference between actual number of points(vehicles) and initial
	  number of points
	=#
	diff = total_points - length(left_ind.starting_points);

	if diff == 0
		# if actual number of points would be the same then simply concatenate arrays
		left_starting_points = vcat(left_starting_points, right_starting_points);
	elseif diff > 0
		#= if actual number of points would be greater than initial number of points copy only
		  part of the right array
		=#
		left_starting_points = vcat(left_starting_points,
			right_starting_points[1:length(right_starting_points) - diff]);
	else
		if keep_number_of_starting_points

			#= if actual number of points would be fewer, collect all unused points,
			  that are not part of the right starting points of the right individual
			=#
			unused_points = LinkedStack{Integer}();
			unused_points_length = 0;
			j = 1;
			for i in point+1:length(left_ind.permutation)

				if j > length(right_starting_points) || right_starting_points[j] != i
					push!(unused_points, i);
					unused_points_length = unused_points_length + 1;
				else
					j = j + 1;
				end
			end

			#= randomly select n points from this unused points, where n is equal
			  to the number of points that need to be added to sum up to initial Number
			  of points
			=#
			add_points = Vector{Int32}();
			for i in 1:abs(diff)
				next = rand(1:unused_points_length);
				p = popAt!(unused_points, next);
				push!(add_points, p);
				unused_points_length = unused_points_length - 1;
			end
			sort!(add_points);

			#= combine the randomly selected starting points and starting points
			  from the right part of the right individual
			=#
			new_right_starting_points = combine_two_ordered_arrays(right_starting_points, add_points);

			left_starting_points = vcat(left_starting_points, new_right_starting_points);
		else
			left_starting_points = vcat(left_starting_points, right_starting_points);
		end
	end

	return Individual(child_perm, left_starting_points, copy(left_ind.starting_points_positions), copy(left_ind.starting_points_directions),
		copy(left_ind.points_positions), copy(left_ind.points_directions), DEFAULT_FITTNESS);
end

function mutate!(ind::Individual)
	pos = rand(1:length(ind.permutation));
	if pos > length(ind.permutation) / 2
		pos_other = rand(1:pos-1);
	else
		pos_other = rand(pos+1:length(ind.permutation));
	end

	swap!(ind.permutation, pos, pos_other);
end

function swap!(arr::Array{T}, i::Integer, j::Integer) where {T<:Number}
	swap = arr[i];
	arr[i] = arr[j];
	arr[j] = swap;
end

function combine_two_ordered_arrays(arr1::Array{T}, arr2::Array{T}) where {T<:Number}
	i,j = 1,1;
	new_arr = [];

	while i <= length(arr1) && j <= length(arr2)
		if arr1[i] > arr2[j]
			push!(new_arr, arr2[j]);
			j = j + 1;
		else
			push!(new_arr, arr1[i]);
			i = i + 1;
		end
	end

	if i <= length(arr1)
		new_arr = vcat(new_arr,
		 arr1[i:length(arr1)]);
	elseif j <= length(arr2)
		new_arr = vcat(new_arr,
		 arr2[j:length(arr2)]);
	end

	return new_arr;
end

function show(ind::Individual, regions::Vector{Region})
	
	plt.figure(1)
    plt.clf()
	ax = plt.gca();
	
	# plot trajectories
	palette = ["red", "blue", "cyan", "green", "purple"];
	
	for i in 1:length(ind.starting_points)
		last = (i == length(ind.starting_points)) ? length(ind.permutation) : ind.starting_points[i+1]-1;
		
		# add points to trajectory
		trajectory = [];
		push!(trajectory, [ind.starting_points_positions[i].x, ind.starting_points_positions[i].y, ind.starting_points_directions[i]]);
		for j in ind.permutation[ind.starting_points[i]:last]
			push!(trajectory, [ind.points_positions[j].x, ind.points_positions[j].y, ind.points_directions[j]]);
		end
		
		# add starting point again to the end of trajectory
		push!(trajectory, [ind.starting_points_positions[i].x, ind.starting_points_positions[i].y, ind.starting_points_directions[i]]);
		
		c = palette[rand(1:length(palette))];
		for j in 1:length(trajectory)-1
			dubins = dubins_shortest_path(trajectory[j], trajectory[j+1], turning_radius)[2];
			samples = dubins_path_sample_many(dubins, turning_radius * 0.1)[2];
			
			x_val = [x[1] for x in samples];
			y_val = [x[2] for x in samples];
			
			plt.plot(x_val, y_val, color=c);
		end
	end
	
	# plot points
	for i in 1:length(ind.starting_points_positions)
		plt.scatter(ind.starting_points_positions[i].x,ind.starting_points_positions[i].y, color="black");
	end
	
	for i in 1:length(ind.points_positions)
		plt.scatter(ind.points_positions[i].x,ind.points_positions[i].y, color="black");
	end
	
	# plot regions
	for i in 1:length(regions)
		center = [regions[i].center.x, regions[i].center.y];
		radius = regions[i].radius;
		circle = Circle(center, radius, facecolor="red", edgecolor="red", 
        linewidth=1, alpha=0.2);
		ax.add_patch(circle)
	end
	
	plt.title("Length = ");
    plt.axis("equal");
    plt.tight_layout();
    plt.pause(0.1);
	plt.savefig("result.pdf");
	
end

function main()
	(number_vehicles, starting_region, regions) = parseInput();
	centers = [regions[i].center for i in 1:length(regions)];
	radiuses = [regions[i].radius for i in 1:length(regions)];
	region_count = length(regions);
	
	population = [generate_individual(region_count, number_vehicles, starting_region.center, centers) for i in 1:population_size];
	
	crossover_count = round(crossover_ratio*population_size);
	crossover_count = convert(Integer, crossover_count);
	if isodd(crossover_count)
			crossover_count = crossover_count - 1;
	end
	
	crossover_point = 0.5;
	crossover_point = round(crossover_point*region_count);
	crossover_point = convert(Integer, crossover_point);
	
	empty_individual = Individual([],[],[],[],[],[], typemax(Float64));
	best_found = empty_individual;
	
	local_regions = Vector{Region}();
	push!(local_regions, starting_region);
	for i in 1:length(regions)
		push!(local_regions, regions[i]);
	end
	regions = local_regions;
	
	for i in 1:number_of_iterations
		# local search
		
		
		for i in 1:length(population)
			local_search!(population[i], regions);
		end
		
		#evaluation
		for i in 1:length(population)
			evaluate_individual!(population[i], regions);
		end
	
		#selection
		sort!(population, by = p -> p.fittness);
		#update best found individual
		if best_found.fittness > population[1].fittness
			best_found = population[1];
			println(best_found.fittness);
		end
		population = population[1:population_size];
		
		pairs = randperm(crossover_count);
		half = round(crossover_count/2);
		half = convert(Integer, half);
		
		lefts = population[pairs[1:half]];
		rights = population[half+1:end];
		
		#crossover
		offspring = [];
		for i in 1:length(lefts)
			child_1 = crossover(lefts[i], rights[i], crossover_point, true);
			if rand() > 0.5
				mutate!(child_1);
			end
			child_2 = crossover(rights[i], lefts[i], crossover_point, true);
			if rand() > 0.5
				mutate!(child_2);
			end
			push!(offspring, child_1);
			push!(offspring, child_2);
		end
		
		population = vcat(population, offspring);
	end 
	show(best_found, regions);
end

main();

export Point, Region

end
