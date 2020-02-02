#!/usr/bin/env julia

module kDTSPN

import Base: show

using Dubins
using Random
using Test
using PyPlot
include("./linked_stack.jl")
include("./FindingCliques.jl")
using .Stack

Circle = matplotlib.patches.Circle

input_file_path = "./problem.txt";
delimiter = " ";
population_size = 200;
number_of_iterations = 1000;

# local search stops when optimized points converge to a value that is
# less than this number away from previously optimized value
local_search_stop_criterion = 0.001;


struct Point{T<:AbstractFloat}
	x::T
	y::T
end

struct Region{T<:AbstractFloat}
	center::Point{T}
	radius::T
end

function show(r::Region, i)
	println("Region $i: Center [$(r.center.x) $(r.center.y)] Radius $(r.radius)");
end

struct Configuration
	points::Dict{Int32, Vector{Region}}
	number_of_regions::Int32
	number_of_starting_points::Int32
end

struct Individual
	permutation::Vector{Int32}
	starting_points::Vector{Int32}
	starting_points_positions::Vector{Point{Float32}}
	points_positions::Vector{Point{Float32}}
end

#=
  Takes input and produces set of regions to visit
=#
function parseInput()
	regions = [];
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
			starting_region = Region(Point(parse(Float32, convert(String, item[2])),parse(Float32, convert(String, item[3]))), parse(Float32, convert(String, item[4])));
		else
			if(length(item) != 4)
				error("Wrong input, error at line $i");
			end
			r = Region(Point(parse(Float32, convert(String, item[2])),parse(Float32, convert(String, item[3]))), parse(Float32, convert(String, item[4])));
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
function evaluate_individual(ind::Individual)

	for i in length(ind.starting_points)
		# calculate path for selected vehicle
		part_path = []
	end
end

#=
  calculates new positions for points of individual, accepts indivudal definition
  and vector of regions, which contains union of starting region with rest of the regions
  regions = starting_region ∪ regions. Updates the supplied individual with new set of points
=#
function local_search!(ind::Individual, regions::Vector{Region})
	
	points = [i + 1 for i in 1:length(ind.permutation)];
	
	#add starting points to sequence of points
	all_points = [];
	for i in 1:length(ind.starting_points)-1
		all_points = vcat(all_points, 1, points[ind.starting_points[i]:ind.starting_points[i+1]]);
	end
	
	all_points = vcat(all_points, 1, points[ind.starting_points[length(ind.starting_points) - 1]:end]);
	println(all_points);
	
	# update starting points indexes 
	starting_points = [ind.starting_points[i]+i-1 for i in 1:ind.starting_points];
	
	len = length(starting_points)
	for i in 1:len
		diff = typemax(Float32);
		current_sequence = (i==len ? all_points[starting_points[i]:end] : all_points[starting_points[i]:starting_points[i+1]]);
		
		# add starting point specific to sequence and rest of the points to current points
		current_points = vcat(ind.starting_points[i], ind.points_positions);
		
		while diff > local_search_stop_criterion
			# while not converged optimized current subsequence
			
			# optimize subsequence
			optimized = [j == length(current_sequence) ? 
			projected_point(ind.points_positions[current_sequence[j-1]], ind.points_positions[current_sequence[1]], 
				ind.points_positions[current_sequence[j]], regions[current_sequence[j]]) :
			j == 1 ? projected_point(ind.points_positions[current_sequence[end]], ind.points_positions[current_sequence[j+1]], 
				ind.points_positions[current_sequence[j]], regions[current_sequence[j]]) :
			projected_point(ind.points_positions[current_sequence[j-1]],ind.points_positions[current_sequence[j+1]], 
				ind.points_positions[current_sequence[j], regions[current_sequence[j]]])
			for j in 1 : length(current_sequence)];
			
			diff = distance(optimized, current_points[current_sequence]);
			current_points[current_sequence] = optimized;
		end
		
		# copy optimized points back to individual
		ind.starting_points[i] = current_points[1];
		
		# start from second point as first i starting point
		for j in 2:length(current_sequence)
			ind.points_positions[current_sequence[j]-1] = current_points[j];
		end
	end
end

#=
  Find projection for point point_c. Where point_a, point_b are neighbourhood points of point point_c.
  point_c is the projected point and region_c is the region of point point_c.
=#
function projected_point(point_a::Point{T}, point_b::Point{T}, point_c::Point{T}, region_c::Region{T}) where T <: AbstractFloat
	
	# for points point_a, point_b find general equation of line
	normal_vector_ab = (point_b.y - point_a.y, point_a.x - point_b.x);
	general_line_equation_ab = (normal_vector_ab[1], normal_vector_ab[2], -normal_vector_ab[1]*point_a.x -normal_vector_ab[2]*point_a.y);
	
	# for line between points point_a, point_b we find general equation for perpendicular line
	# to this line intersecting point_c
	normal_vector_c_ab = (point_b.x - point_a.x, point_a.y - point_b.y);
	general_line_equation_c_ab = (normal_vector_c_ab[1], normal_vector_c_ab[2], -normal_vector_c_ab[1]*point_c.x 
	-normal_vector_c_ab[2]*point_c.y);
	
	#find intersection of line between points point_a and point_b and perpendicular line
	# to this line containing point_c
	intersection_y = (general_line_equation_ab[1]*general_line_equation_c_ab[3] - general_line_equation_ab[3]*general_line_equation_c_ab[1]) /
		(general_line_equation_c_ab[1]*general_line_equation_ab[2] - general_line_equation_ab[1]*general_line_equation_c_ab[2]);
	
	
	
	intersection_x = general_line_equation_ab[1] != 0 ? 
		((-general_line_equation_ab[2]*intersection_y - general_line_equation_ab[3]) / general_line_equation_ab[1]) :
		((-general_line_equation_c_ab[2]*intersection_y - general_line_equation_c_ab[3]) / general_line_equation_c_ab[1]);
	
	# find out if intersection lies on the line between points point_a and point_b
	is_between = max(distance(intersection_x, intersection_y, point_a.x, point_a.y),
		distance(intersection_x, intersection_y, point_b.x, point_b.y)) <=  distance(point_a.x, point_a.y, point_b.x, point_b.y);
	if is_between
		# find intersection of perpendicular line with region
		points = find_intersection_of_line_and_circle(general_line_equation_c_ab[1],general_line_equation_c_ab[2],general_line_equation_c_ab[3],
			region_c.center.x, region_c.center.y, region_c.radius);
		# return point closer to intersection point of two lines
		if distance(points[1][1], points[1][2],intersection_x, intersection_y) < distance(points[2][1], points[2][2],intersection_x, intersection_y)
			return points[1];
		else
			return points[2];
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
		
		# return point closer to center of points point_a and point_b
		if distance(points[1][1], points[1][2],center_x, center_y) < distance(points[2][1], points[2][2],center_x, center_y)
			return points[1];
		else
			return points[2];
		end
	end
end

function testProjectPoint()
	r = Region(Point(5.0,5.0), 2.0);
	a = Point(3.0,0.0);
	b = Point(6.0,0.0);
	c = Point(5.0, 5.0);
	
	projected = projected_point(a, b, c, r);
	println(projected);
	
	
	
	pyplot();
	xy = Vector{Float32}();
	xy.push(r.center.x);
	xy.push(r.center.y);
    circle = Circle(xy, r.radius, facecolor="yellow", edgecolor="orange",
       linewidth=1, alpha=0.2);
	
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
		return [(x1,y1),(x2,y1)];
	end
	a = b1^2/a1^2 + 1;
	b = -2*b1*c1/a1^2 - 2*s1*b1/a1 - 2*s2;
	c = c1^2/a1^2+2*s1*c1/a1 + s1^2 + s2^2 - r^2;
	(y1, y2) = solve_quadratic_equation(a, b, c);
	if y1 == y2
		b = -2*s1;
		c = s1^2 + y1^2 - 2*s2*y1 + s2^2 - r^2;
		(x1, x2) = solve_quadratic_equation(1, b, c);
		return [(x1,y1),(x2,y1)];
	else
		x1 = (b1*y1 - c1) / a1;
		x2 = (b1*y2 - c1) / a1;
		return [(x1,y1),(x2,y2)];
	end
end

function distance(a::Vector{Point{Float32}},b::Vector{Point{Float32}})
	return sum([distance(a[i].x,a[i].y,b[i].x,b[i].y) for i in 1:length(a)])
end

function distance(x1, y1, x2, y2)
	return √((x1 - x2) ^ 2
		+ (y1 - y2) ^ 2);
end

function distance(a::Point{Float32}, b::Point{Float32})
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

function show(ind::Individual)
	str = "Permutation [";
	for i in 1:length(ind.permutation)-1
		str = string(str, ind.permutation[i], ", ")
	end
	str = string(str, ind.permutation[length(ind.permutation)], "] StartingPoints: [");
	for i in 1:length(ind.starting_points)-1
		str = string(str, ind.starting_points[i], ", ")
	end
	str = string(str, ind.starting_points[length(ind.starting_points)], "]");
	println(str);
end

function generate_individual(total_points_count::Integer, starting_points_count::Integer, starting_point::Point{Float32}, points::Vector{Point{Float32}})

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
	
	return Individual(permutation, taken, starting_points, points[permutation]);
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

	return Individual(child_perm, left_starting_points);
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

function testcrossoverproducesvalidoffspring()
	perm1 = [1,2,3,4,5,6];
	starting_points1 = [1,3,4];
	ind1 = Individual(perm1, starting_points1);

	perm2 = [2,3,4,5,6,1];
	starting_points2 = [1,5,6];
	ind2 = Individual(perm2, starting_points2);

	st1 = LinkedStack{Integer}();
	for i in 1:length(perm1)
		push!(st1, perm1[i]);
	end

	st2 = LinkedStack{Integer}();
	for i in 1:length(perm2)
		push!(st2, perm2[i]);
	end

	child1 = crossover(ind1, ind2, 3, true);
	child2 = crossover(ind2, ind1, 3, true);

	println(child1.permutation);
	println(child1.starting_points);
	println(child2.permutation);
	println(child2.starting_points);

	@assert length(child1.permutation) == length(perm1)
	for i in 1:length(child1.permutation)
		@assert pop!(st1, child1.permutation[i]);
	end
	@assert isempty(st1);

	@assert length(child2.permutation) == length(perm2)
	for i in 1:length(perm2)
		@assert pop!(st2, child2.permutation[i]);
	end
	@assert isempty(st2);

	st_startingpoints1 = LinkedStack{Integer}();
	@assert length(child1.starting_points) == 3;
	for i in 1:length(child1.starting_points)
		@assert !pop!(st_startingpoints1, child1.starting_points[i]);
		push!(st_startingpoints1, child1.starting_points[i]);
	end

	st_startingpoints2 = LinkedStack{Integer}();
	@assert length(child2.starting_points) == 3;
	for i in 1:length(child2.starting_points)
		@assert !pop!(st_startingpoints2, child2.starting_points[i]);
		push!(st_startingpoints2, child2.starting_points[i]);
	end
end

function test_crossover_bulk()
	for i in 1:10000
		test_crossover();
	end
end

function test_crossover()

	total_points = 10;
	starting_points = 3;

	ind1 = generate_individual(total_points, starting_points);
	ind2 = generate_individual(total_points, starting_points);

	keep = true;
	point = 5;

	child = crossover(ind1, ind2, point, keep);

	# verify that left part of the child remains the same as in left individual
	@assert isequal(ind1.permutation[1:point], child.permutation[1:point]);

	left_start_points = filter(x -> x <= point, ind1.starting_points);
	child_left_start_points = filter(x -> x <= point, child.starting_points);

	@assert isequal(left_start_points, child_left_start_points);


	# check that child contains all points
	points = [false for x in 1:total_points];

	@assert isequal(length(child.permutation), total_points);
	for i in 1:length(child.permutation)
		points[child.permutation[i]] = true;
	end

	for i in 1:length(points)
		@assert points[i]
	end

	# check starting points of the child
	stack = LinkedStack{Integer}();
	@assert isequal(length(child.starting_points), starting_points)
	for i in 1:length(child.starting_points)
		@assert pop!(stack, child.starting_points[i]) == false
		push!(stack, child.starting_points[i]);
	end
end

function test_mutate()
	total_points = 10;
	starting_points = 3;

	ind1 = generate_individual(total_points, starting_points);

	mutate!(ind1);
	stack = LinkedStack{Integer}();
	for i in 1:length(ind1.permutation)
		push!(stack, i)
	end

	# test individual contains all points excatly once
	for i in 1:length(ind1.permutation)
		@assert pop!(stack, ind1.permutation[i])
	end

	# test starting points do not have duplicates
	stack = LinkedStack{Integer}();
	for i in 1:length(ind1.starting_points)
		@assert pop!(stack, ind1.starting_points[i]) == false;
		push!(stack, ind1.starting_points[i])
	end

end

function test_generates_individual()
	total_points = 10;
	starting_points = 3;
	ind = generate_individual(total_points, starting_points);

	stack = LinkedStack{Integer}();
	for i in 1:total_points
		push!(stack, i);
	end

	for i in 1:total_points
		@assert pop!(stack, ind.permutation[i]);
	end

	stack = LinkedStack{Integer}();
	for i in 1:length(ind.starting_points)
		@assert pop!(stack, ind.starting_points[i]) == false
		push!(stack, ind.starting_points[i]);
	end
end

function test_combine_arrays()
	arr1 = [1,4,5,6];
	arr2 = [2,7,8];

	arr = combine_two_ordered_arrays(arr1, arr2);
	@assert arr[1] === 1;
	@assert arr[2] === 2;
	@assert arr[3] === 4;
	@assert arr[4] === 5;
	@assert arr[5] === 6;
	@assert arr[6] === 7;
	@assert arr[7] === 8;
end

function test_overlapping_regions()

	a = Region(Point(1.0,1.0), 3.0);
	b = Region(Point(2.0,2.0), 3.0);
	@test overlaps(a, b);

	a = Region(Point(1.0,1.0), 3.0);
	b = Region(Point(10.0,10.0), 2.0);

	@test !overlaps(a, b);
end

function main()
	#=(number_vehicles, starting_region, regions) = parseInput();
	centers = [regions[i].center for i in 1:length(regions)];
	radiuses = [regions[i].radius for i in 1:length(regions)];

	a = Vector{Point{Float32}}();
	push!(a, Point(3.0,4.0);
	push!(a, Point(3.0,4.0);
	push!(a, Point(3.0,4.0);
	a[1] = Point(3.0,5.0);
	println(a[1]);
	
	#= population = [generate_individual(length(regions), number_vehicles, starting_region.center, centers) for i in 1:population_size];
	
	for i in 1:number_of_iterations
		# local search
		local_regions = vcat(starting_region, regions);
		for i in 1:length(population)
			local_search!(population[i], local_regions);
		end
		#evaluation
		#selection
		#crossover
		#mutace
	end =#
	=#
	
end
	
testProjectPoint();

#main();

#test_overlapping_regions();
#test_mutate();
#@time test_crossover_bulk();
#test_generates_individual();
#testcrossoverproducesvalidoffspring();
#test_combine_arrays();

#s = parseInput();
#for i in 1:length(s)
#	show(s[i],i);
#end

#A = overlaps(s);
#C = populate_cliques(A);
#println(C);
#squash_regions(A);


end
