#!/usr/bin/env julia

module kDTSPN

using Dubins

struct Individual
	Sequence::Vector{Int32}
	StartingPoints::Vector{Int32}
end

function crossover(ind1::Individual, ind2::Individual, point::Int32)
	# check that individual length does not exceed the point of crossover
	if	point >= length(ind1.Sequence) || (point == length(ind1.Sequence) - 1 && ind2[length] == 1)
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
	for i in point+1:length(right_ind1)
		for j in point+1:length(right_ind1)
			if(ind2[i] == right_ind1[j])
				right_ind1[j] = 0;
			else
				push!(duplicates, i);
			end
		end
		push!(left_ind1, ind2[i]);
	end
	
	
end

function contains()
	
end

function doSomething()
	a = Individual([3,4], [3,5]);
	#a.Sequence = [3,4];
	#a.StartingPoints = [3,5];
	println(a.Sequence, a.StartingPoints)
end

doSomething();

end