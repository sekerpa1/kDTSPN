#!/usr/bin/env julia

module kDTSPN

using Dubins

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

function doSomething()
end

doSomething();

end