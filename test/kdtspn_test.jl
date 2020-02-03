#!/usr/bin/env julia

using Test


#include("./../src/linked_stack.jl")
include("./../src/kdtspn.jl")
using .kDTSPN

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


function testProjectPoint()
	r = Region(Point(8.0,5.0), 2.0);
	a = Point(3.0,2.0);
	b = Point(6.0,2.0);
	c = Point(8.0,5.0);
	
	projected = projected_point(a, b, c, r);
	println(projected);
	
	xy = [r.center.x, r.center.y];
	
	plt.figure(1)
    plt.clf()
	ax = plt.gca();
    circle = Circle(xy, r.radius, facecolor="red", edgecolor="red", 
        linewidth=1, alpha=0.2);
    ax.add_patch(circle);
	plt.scatter([a.x, b.x, c.x],[a.y, b.y, c.y]);
	plt.scatter(projected[1],projected[2]);
	#plt.scatter([a.x, b.x, c.x],[a.y, b.y, c.y],markersize=3,alpha=.8,legend=false);
	#plt.plot(c.x, c.y, specs = "b");
	#plt.plot(projected[1], projected[2], specs = "b");
	
	plt.title("Length = ")
    plt.axis("equal")
    plt.tight_layout()
    plt.pause(0.1)

    # save into a file
    #plt.savefig("something.pdf")
end

test_overlapping_regions();
test_mutate();
@time test_crossover_bulk();
test_generates_individual();
testcrossoverproducesvalidoffspring();
test_combine_arrays();

