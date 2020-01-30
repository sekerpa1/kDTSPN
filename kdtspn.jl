#!/usr/bin/env julia

module kDTSPN

import Base: show

#using Dubins
using Random
using Test
include("./linked_stack.jl")
using .Stack

struct Point{T<:AbstractFloat}
	x::T
	y::T
end

struct Region{T<:AbstractFloat}
	center::Point{T}
	radius::T
end

struct Configuration
	points::Dict{Int32, Vector{Region}}
	number_of_regions::Int32
	number_of_starting_points::Int32
end

struct Individual
	permutation::Vector{Int32}
	starting_points::Vector{Int32}
end

function overlaps(a::Region, b::Region)

	center_distance = √(abs(a.center.x - b.center.x) ^ 2
		+ abs(a.center.y - b.center.y) ^ 2);

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

function generate_individual(total_points::Integer, starting_points::Integer)

	if starting_points < 2
		throw(ArgumentError("individual should have at least two starting points"))
	end

	if total_points < 2
		throw(ArgumentError("Individual should have at least two total points"))
	end


	permutation = randperm(total_points);

	taken = [];
	free = LinkedStack{Integer}();
	free_size = total_points;
	for i in 1:total_points
		push!(free, i);
	end

	for i in 1:starting_points
		pos = rand(1:free_size);
		push!(taken, popAt!(free, pos));
		free_size = free_size - 1;
	end

	sort!(taken);

	return Individual(permutation, taken);
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

function crossover(ind1::Vector{Int32}, ind2::Vector{Int32}, point::Int32)
	# check that individual length does not exceed the point of crossover
	if	point >= length(ind1) ||
		point == length(ind1) - 1 && ind2[length(ind1)] == 1 ||
		point == length(ind1) - 2 && ind1[length(ind1)-2] == 1 && (ind2[length(ind1) - 1] == 1 || ind2[length(ind1)] == 1) ||

		error("Crossover point outside possible range");
	end

	# if crossover would cause empty inner sequence at start of ind2 of first child
	if ind1[point] == 1 && ind2[point+1] == 1
		#swap two beginning points in individual2
		swap = ind2[point+1];
		ind2[point+1] = ind2[point+2];
		ind2[point+2] = swap;
	end

	# if crossover would cause empty inner sequence at start of ind1 of second child
	if ind2[point] == 1 && ind1[point+1] == 1
		#swap two beginning points in individual1
		swap = ind1[point+1];
		ind1[point+1] = ind1[point+2];
		ind1[point+2] = swap;
	end


	left_ind1 = copy(ind1[1:point]);

	right_ind1 = copy(ind1[point+1:end]);
	duplicates_ind1 = [];

	for i in point+1:length(ind1)
		not_found = 1;
		for j in 1:length(right_ind1)
			if ind2[i] == right_ind1[j]
				right_ind1[j] = 0;
				not_found = 0;
				break;
			end
		end
		if not_found == 1
			push!(duplicates, i);
		end
		push!(left_ind1, ind2[i]);
	end

	# replace duplicate in child with available points
	for i in 1:length(duplicates)
		for j in 1:length(right_ind1)
			if (right_ind1[j] != 0 && (right_ind1[j] != 1 || left_ind1[duplicates[i]-1] != 1))
				left_ind1[duplicates[i]] = right_ind1[j];
				break
			end
		end
	end

	left_ind2 = copy(ind2[1:point]);
	right_ind2 = copy(ind2[point+1:end]);
	duplicates = [];
	for i in point+1:length(ind1)
		not_found = 1;
		for i in 1:length(right_ind2)
			if(ind1[i] == right_ind2[j])
				right_ind2[j] = 0;
				not_found = 0;
				break;
			end
		end
		if not_found == 1
			push!(duplicates, i);
		end
		push!(left_ind2, ind1[i]);
	end

	# replace duplicate in child with available points
	for i in 1:length(duplicates)
		for j in 1:length(right_ind2)
			if (right_ind2[j] != 0 && (right_ind2[j] != 1 || left_ind2[duplicates[i]-1] != 1))
				left_ind2[duplicates[i]] = right_ind2[j];
				break;
			end
		end
	end

	return [left_ind1, left_ind2];
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


#test_overlapping_regions();
test_mutate();
test_crossover_bulk();
test_generates_individual();
#testcrossoverproducesvalidoffspring();
#test_combine_arrays();
end
