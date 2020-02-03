module Stack

using Test

import Base: isempty, empty!, length, first, last, iterate, isdone
import Base: eltype, pushfirst!, popfirst!, deleteat!
import Base: show, println, indexin, lastindex, push!, pop!, insert!, splice!, keys
import Base: in , append!, prepend!
import Base: getindex, setindex!

abstract type AbstractNode{T} end
abstract type AbstractStack{T} end

mutable struct StackNode{T} <: AbstractNode{T}
	prev::StackNode{T}
	data::T
	StackNode{T}() where T =(x = new(); x.prev = x;)
	StackNode{T}(p,d) where T = new(p,d)
end

StackNode(prev::StackNode{Integer}, data::Integer) = StackNode{Integer}(prev, data);
StackNode(prev::StackNode{Int64}, data::Int64) = StackNode{Int64}(prev, data);
StackNode(prev::StackNode{Array{Integer}}, data::Array{Int64}) = StackNode{Array{Integer}}(prev, convert(Array{Integer}, data));
StackNode(prev::StackNode{Array{Integer}}, data::Array{Integer}) = StackNode{Array{Integer}}(prev, data);
StackNode(prev::StackNode{AbstractFloat}, data::AbstractFloat) = StackNode{AbstractFloat}(prev, data);

mutable struct LinkedStack{T} <: AbstractStack{T}
	head::StackNode{T}
	LinkedStack{T}() where T = new(StackNode{T}())
end

isempty(s::AbstractStack) = (s.head.prev === s.head);

function first(s::AbstractStack)
	return s.head.data;
end

function last(s::AbstractStack)
	head = s.head;
	element = head.data;
	while head.prev != head.prev.prev
		element = head.prev.data;
		head = head.prev;
	end
	
	return element;
end

function empty!(s::AbstractStack)
	while s.head != s.head.prev
		s.head = s.head.prev;
	end
end

function push!(s::AbstractStack, item)
	a = StackNode(s.head, item);
	s.head = a;
end

function pop!(s::AbstractStack)
	if(s.head === s.head.prev)
		throw(UndefRefError());
	else
		item = s.head.data;
		s.head = s.head.prev;
		return item;
	end
end

function pop!(s::AbstractStack, item)

	if(s.head === s.head.prev)
		return false;
	end

	if(s.head.data == item)
		s.head = s.head.prev;
		return true;
	end

	head = s.head;
	while(true)
		if(head.prev === head.prev.prev)
			return false;
		end

		if(head.prev.data == item)
			head.prev = head.prev.prev;
			return true;
		end

		head = head.prev;
	end
end

function popAt!(s::AbstractStack, pos)
    if pos < 1
        throw(BoundsError());
    end

    if pos == 1
        if s.head == s.head.prev
            throw(BoundsError());
        end
        item = s.head.data;
        s.head = s.head.prev;
        return item;
    end

    head = s.head;
    while pos != 2
        if head == head.prev
            throw(BoundsError());
        end
        head = head.prev;
        pos = pos - 1;
    end

    if head == head.prev
        throw(BoundsError())
    end

    item = head.prev.data;
    head.prev = head.prev.prev;

    return item;
end

function first(s::AbstractStack)
    if s.head == s.head.prev
        throw(UndefRefError());
	end
    return s.head.data;
end

function length(s::AbstractStack)
	length = 0;
	head = s.head;
	while(head != head.prev)
		length = length + 1;
	end

	return length;
end

#=
    inserts the item into the stack, after specified element.
    If stack does not contain element, do not add item
=#
function insert!(s::AbstractStack, item, element)

    cur = s.head;
    while cur != cur.prev
        if cur.data == item
            prev = cur.prev;
            new_item = StackNode(prev, item);
            cur.prev = new_item;
            return true;
        end
    end

    return false;
end

#=
	insert item into the stack, keep order of elements
=#
function insertAndKeepOrder!(s::AbstractStack, item)

	if s.head === s.head.prev || s.head.data >= item
		prev = s.head;
		head = StackNode(prev, item);
		s.head = head;
		return;
	end

	head = s.head;
	while true
		if head.prev === head.prev.prev
			break;
		end

		if head.prev.data >= item
			break;
		end

		head = head.prev;
	end

	element = StackNode(head.prev, item);
	head.prev = element;
	return;
end

function show(s::AbstractStack)

	if s.head === s.head.prev
		return "empty"
	end
	str = string("[", s.head.data);
	head = s.head.prev;
	while head.prev != head
		str = string(str, ",", head.data);
		head = head.prev;
	end
	str = string(str, "]");
	println(str);
end

function length(s::AbstractStack)
	length = 0;
	head = s.head;
	while head != head.prev
		length = length + 1;
		head = head.prev;
	end

	return length;
end

function tests()
	s = StackNode{Int64}()
	@test s === s.prev
	s2 = StackNode{Int64}(s, 10)
	@test (s2.data === 10 && s2.prev === s2.prev.prev)
	st = LinkedStack{Int64}();
	@test st.head === st.head.prev
	@test isempty(st);

	# test pop for empty stack returns error
	@test_throws UndefRefError pop!(st);

	for i in 1:10
		push!(st, i);
	end
	@test pop!(st, 5);
	@test length(st) == 9;
	@test pop!(st, 4);
	@test pop!(st, 3);
	@test length(st) == 7;
	@test pop!(st, 20) == false;
	@test length(st) == 7;

    sta = LinkedStack{Int64}();
    for i in 1:10
		push!(sta, 10 - i + 1);
	end

    @test popAt!(sta, 2) == 2;
    @test popAt!(sta, 2) == 3;
    @test popAt!(sta, 2) == 4;

	stack = LinkedStack{Int64}();
	for i in 1:5
		push!(stack, i);
	end

	@test length(stack) == 5;
	pop!(stack);
	@test length(stack) == 4;
	push!(stack, 6);
	push!(stack, 7);
	@test length(stack) == 6;
	testInsertAndKeepOrder();

	s = StackNode{Float64}();
	st = LinkedStack{AbstractFloat}();
	push!(st, 4.5);
	@test pop!(st, 4.5);
	
	s = StackNode{Array{Integer}}();
	st = LinkedStack{Array{Integer}}();
	push!(st, [1]);
	push!(st, [1,2,3]);
	push!(st, [1,2,3,4]);
	push!(st, [1,2,3,7]);
	@test pop!(st, [1,2,3]);
	@test pop!(st, [1,2,3,4]);
	@test pop!(st, [1,2,3,7]);
	
	st = LinkedStack{Integer}();
	push!(st,3);
	empty!(st);
	@test isempty(st);
	push!(st,3);
	push!(st,5);
	@test last(st) == 3;
	@test first(st) == 5;
end

function testInsertAndKeepOrder()
	stack = LinkedStack{Integer}();
	for i in 1:1000
		insertAndKeepOrder!(stack, rand(1:1000));
	end

	arr = [];
	for i in 1:1000
		push!(arr, pop!(stack));
	end

	for i in 1:999
		@assert arr[i] <= arr[i+1]
	end
end

#@time @testset "Stack" begin tests() end

export StackNode, LinkedStack

export popAt!

end
