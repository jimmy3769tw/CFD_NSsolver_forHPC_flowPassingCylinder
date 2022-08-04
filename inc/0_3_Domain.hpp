#pragma once
#include<vector>
#include<iostream>

class divideLocal{

	public:
	int word_size;
	int i_begin;
	int j_begin;
	int k_begin;

	int i_endof;
	int j_endof;
	int k_endof;

	std::vector<int> i_begin_table;
	std::vector<int> j_begin_table;
	std::vector<int> k_begin_table;

	std::vector<int> i_endof_table;
	std::vector<int> j_endof_table;
	std::vector<int> k_endof_table;

	std::vector<int> i_length_table;
	std::vector<int> j_length_table;
	std::vector<int> k_length_table;


	void show_begin(){
		cout << "show begin";

		for (auto x :i_begin_table){
			std::cout << x << ", ";
		}

		for (auto x :j_begin_table){
			std::cout << x << ", ";
		}

		for (auto x :k_begin_table){
			std::cout << x << ", ";
		}
	}

	void initTable(	
		std::vector<int> gridASize, 
		std::vector<int> dims
	){
		this->word_size = dims.at(0) * dims.at(1) * dims.at(2);
		this->divideDomain(dims.at(0), gridASize.at(0)-4, i_begin_table, i_endof_table, i_length_table);
		this->divideDomain(dims.at(1), gridASize.at(1)-4, j_begin_table, j_endof_table, j_length_table);
		this->divideDomain(dims.at(2), gridASize.at(2)-4, k_begin_table, k_endof_table, k_length_table);
	}

	void initLocal(int i, int j, int k){

		this->i_begin = this->i_begin_table.at(i);
		this->j_begin = this->j_begin_table.at(j);
		this->k_begin = this->k_begin_table.at(k);

		this->i_endof = this->i_endof_table.at(i);
		this->j_endof = this->j_endof_table.at(j);
		this->k_endof = this->k_endof_table.at(k);
	}

	private:

	void divideDomain(
		int dims,
		int n,
		std::vector<int> &begin ,
		std::vector<int> &endof ,
		std::vector<int> &length
	)
	{
		if (dims == 0){
			throw std::invalid_argument( "divideDomain !!!!\n");
		}

		begin.resize(dims);

		endof.resize(dims);
		
		length.resize(dims);

		for (int i = 0; i < dims; ++i)
		{
			length[i] = (n - n % dims) / (dims);
		}

		for (int i = 0; i < n % dims; ++i)
		{
			length[i] += 1;
		}

		int gC = 2;

		begin[0] = gC;

		for (int i = 0; i < dims - 1; i++)
		{
			begin[i + 1] = begin[i] + length[i];
		}

		for (int i = 0; i < dims; i++)
		{
			endof[i] = begin[i] + length[i];
		}

	}

};


