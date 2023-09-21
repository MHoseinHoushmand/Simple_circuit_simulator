// Simple_circuit_simulator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <string>
#include <iomanip>
using namespace std;
template <class type = double> class Array {
public:

	Array() : Cell(new type *[1]), row(1), col(1) { this->Cell[0] = new type[1]; }
	Array(int argrow, int argcol) : Cell(new type *[row]), row(argrow), col(argcol) {
		for (int i = 0; i < this->row; i++)
			this->Cell[i] = new type[this->col];
	}
	Array(const Array < type > &other) {
		this->row = other.row;
		this->col = other.col;
		this->Cell = new type *[this->row];
		for (int i = 0; i < this->row; i++) {
			this->Cell[i] = new type[this->col];
			for (int j = 0; j < this->col; j++)
				this->Cell[i][j] = other.Cell[i][j];
		}
	}


	~Array() {
		delete[] this->Cell;
	}


	inline type* &operator [ ](int argrow) { return this->Cell[argrow]; }

	Array < type > operator *(Array < type > other) {
		Array < type > _result(this->row, other.col);

		if (this->col == other.row) {
			for (int i = 0; i < _result.row; i++)
				for (int j = 0; j < _result.col; j++) {
					_result[i][j] = 0.0;
					for (int k = 0; k < this->col; k++)
						_result[i][j] += this->Cell[i][k] * other.Cell[k][j];
				}
			return _result;
		}
	}

	/* Algebra Methods */

	// return minor of this
	Array <type> minor(int x, int y) {
		int x1 = 0, y1 = 0;
		Array <type> _result(this->row - 1, this->col - 1);

		for (int i = 0; i < this->row; i++) {
			if (i == x) {
				x1 = 1;
				continue;
			}

			for (int j = 0; j < this->col; j++) {
				if (j == y) {
					y1 = 1;
					continue;
				}
				_result[i - x1][j - y1] = this->Cell[i][j];
			}
			y1 = 0;
		}
		return _result;
	}

private:
	class POINT {
	public:
		char option;
		int  coordinates;
		POINT(char op, int co) :option(op), coordinates(co) {}
	};

public:
	// return determinant of this
	inline long double determinant() {
		POINT coor = bestChoice();
		long double sum = 0.0;
		if (coor.option == 'r')
			if (this->row == 1)
				return this->Cell[0][0];
			else {
				for (int c = 0; c < this->col; c++)
					if (this->Cell[coor.coordinates][c])
						sum += ((coor.coordinates + c) % 2 ? -1 : 1) * this->Cell[coor.coordinates][c] * this->minor(coor.coordinates, c).determinant();
			}
		else
			if (this->col == 1)
				return this->Cell[0][0];
			else {
				for (int r = 0; r < this->row; r++)
					if (this->Cell[r][coor.coordinates])
						sum += ((coor.coordinates + r) % 2 ? -1 : 1) * this->Cell[r][coor.coordinates] * this->minor(r, coor.coordinates).determinant();
			}
			return sum;
	}

	/*
	search for the line or column that
	has most ZERO element
	*/
	POINT bestChoice() {
		int RowZeroCount = 0, _row = 0;
		for (int i = 0; i < this->row; i++) {
			int temp = 0;
			for (int j = 0; j < this->col; j++)
				temp += this->Cell[i][j] ? 0 : 1;
			RowZeroCount = temp > RowZeroCount ? temp : RowZeroCount;
			_row = temp > RowZeroCount ? i : _row;
		}

		int ColZeroCount = 0, _col = 0;
		for (int j = 0; j < this->col; j++) {
			int temp = 0;
			for (int i = 0; i < this->row; i++)
				temp += this->Cell[i][j] ? 0 : 1;
			ColZeroCount = temp > ColZeroCount ? temp : ColZeroCount;
			_col = temp > ColZeroCount ? j : _col;
		}

		return RowZeroCount > ColZeroCount ? POINT('r', _row) : POINT('c', _row);
	}

	// return transported of this
	inline Array < type > transport() {
		Array < type > temp(this->col, this->row);
		for (int i = 0; i < this->row; i++)
			for (int j = 0; j < this->col; j++)
				temp[j][i] = this->Cell[i][j];
		return temp;
	}

	// return adjacent of this
	inline Array < type > adj() {
		Array < type > temp(this->row, this->col);

		for (int i = 0; i < this->row; i++)
			for (int j = 0; j < this->col; j++)
				temp.Cell[i][j] = (((i + j + 2) % 2) ? 1 : -1) * this->minor(i, j).determinant();

		return temp.transport();
	}

	// return invert of this = (this ^ -1)
	inline Array < type > invers() {
		Array < type > temp(this->adj());

		for (int i = 0; i < this->row; i++)
			for (int j = 0; j < this->col; j++)
				temp.Cell[i][j] *= -1.0 / this->determinant();

		return temp;
	}

	/*
	display this
	x	y	z
	a	b	c
	d	e	f
	*/
	bool display() {
		if (this->row < 1 || 1 > this->col)
			return 0;
		else {
			for (int i = 0; i < this->row; i++) {
				for (int j = 0; j < this->col; j++)
					cout << this->Cell[i][j] << '\t';
				cout << '\n';
			}
			cout << '\n';
		}
		return 1;
	}

	// read this
	bool read() {
		for (int i = 0; i < this->row; i++)
			for (int j = 0; j < this->col; j++)
				cin >> this->Cell[i][j];
		return 1;
	}

	void readEquation(float** A, float** B) {
		char constant = 'a', var = 't';
	//	cout << "input each equation like : ";
	//	for (int i = 0; i < (this->col - 1); i++) cout << constant++ << "." << var++ << "   ";
	//	cout << "= Y" << endl;

		for (int i = 0; i < this->row; i++)
			for (int j = 0; j < this->col; j++) {
				if (j == this->col - 1)
					Cell[i][j] = B[i][0];
				else
				    Cell[i][j] = A[i][j];
			}
	}
/*	void writeEquation() {
		char constant = 't';

		cout << "Equations System" << endl;
		for (int i = 0; i < this->row; i++) {
			for (int j = 0; j < (this->col - 1); j++) cout << Cell[i][j] << constant++ << "   ";
			cout << "=\t" << Cell[i][this->col - 1] << endl;
			constant = 't';
		}
		cout << endl;
	}*/

	/*   Methods.Get.Set
	*   ...  @protected fields
	*/

	void SetRow(int argrow) {
		delete[] this->Cell;
		this->Cell = new type *[this->row = argrow];
	}

	void SetCol(int argcol) {
		for (int i = 0; i < this->row; i++) this->Cell[i] = new type[this->col = argcol];
	}

	int Row() { return this->row; }
	int Col() { return this->col; }

private:
	int row, col;
	type **Cell;

};

//Main Code

class Resistor
{
private:
	float Value;
	float Voltage;
	float current;
	bool Is_or_not;
public:
	void Set_Voltage(float v) {
		this->Voltage = v;
	}
	void Set_value(float v) {
		this->Value = v;
	}
	void Set_current(float v) {
		this->current = v;
	}
	float Get_Voltage() {
		return Voltage;
	}
	float Get_current() {
		return current;
	}
	float Get_value() {
		return Value;
	}
	void set_Is_or_not(bool n) {
		Is_or_not = n;
	}
	bool get_Is_or_not() {
		return Is_or_not;
	}
};
class Voltage_source
{
private:
	float value;
public:
	void Set_value(float v) {
		this->value = v;
	}
	float Get_value() {
		return this->value;
	}


};
class Current_source
{
private:
	float value;
public:
	void Set_value(float v) {
		this->value = v;
	}
	float Get_value() {
		return this->value;
	}



};
class Node
{
private:
	int Code;
	int size;
	bool reference;
	float Voltage;
public:
	Node() {
		reference = false;
		size = 0;
		Code = 0;
		Voltage = 0;
	}
	int Get_Code() {
		return Code;
	}
	void Set_Code(int c) {
	    Code=c;
	}
	bool is_reference() {
		return reference;
	}
	void set_reference(bool r) {
		reference = r;
	}
	int get_size() {
		return size;
	}
	void size_plusplus () {
		size++;
	}
	float get_Voltage() {
		return Voltage;
	}
	void set_Voltage(float V) {
		 Voltage=V;
	}
};



class Branch
{
private:
	Node *first_Node;
	Node *final_Node;
	Resistor* resistor;
	Current_source* curent_source;
	Voltage_source* voltage_source;
	float voltage;

public:
	Node* Get_first_Node() {
		return first_Node;
	}
	Node* Get_final_Node() {
		return final_Node;
	}
	void Set_first_Node(Node* fn) {
	    first_Node=fn;
	}
	void Set_final_Node(Node* fn) {
	    final_Node=fn;
	}
	Resistor* get_resistor() {
		return this->resistor;
	}
	void set_resistor(Resistor* r) {
		resistor = r;
	}
	Current_source* get_curent_source() {
		return curent_source;
	}
	void set_curent_source(Current_source* c) {
		curent_source = c;
	}
	Voltage_source* get_voltage_source() {
		return voltage_source;
	}
	void set_voltage_source(Voltage_source* v) {
		voltage_source = v;
	}
	void set_voltage(float v) {
		voltage = v;
	}
	float get_voltage() {
		return voltage;
	}



};
template <class type>  inline Array < type > InvertSolutions(Array < type >);
class Circuit
{
private:
	int Number_of_Nodes;
	int Number_of_Branchs;
	Node * reference;
	Branch * branchlist;
	Node * Nodelist;
	int** Incidence_matrix;
	float** Conductunce_matrix;
	int** Incidence_transpose_matrix;
	float* Voltage_Source_of_Branches;
	float* Current_source_of_Branches;
	float **IC_matrix;
	float **ICIt_matrix;
	float **ICSv_matrix;
	float **ISc_matrix;
	float **Right_Equation_matrix;
	float **solved_matrix;


public:
	Circuit() {
		int i,first,final;
		Resistor*res;
		Voltage_source* Vsource;
		Current_source*Csource;

		float value;
		string answer;
		cout << "How many Nodes do we have?";
		cin >> Number_of_Nodes;
		cout << "How many Branchs do we have?";
		cin >> Number_of_Branchs;
		Nodelist = new Node[Number_of_Nodes];
		branchlist = new Branch[Number_of_Branchs];
		for (i = 0; i < Number_of_Nodes; i++)
			Nodelist[i].Set_Code(i);
		///////////////////////////////////////////////
//		Nodelist[0].set_reference(true);
		reference = &Nodelist[0];
		for (i = 0; i < Number_of_Branchs; i++) {
			cout << "/////////////////////////////////////////" << endl;
			cout << "Enter first & final Node for Branch " << i <<":";
			cin >> first;
			Nodelist[first].size_plusplus();
			if (Nodelist[first].get_size() >= reference->get_size()) {
			//	Nodelist[first].set_reference(true);
				reference = &Nodelist[first];
			}
			cin >> final;
			Nodelist[final].size_plusplus();
			if (Nodelist[final].get_size() >= reference->get_size()) {
			//	Nodelist[final].set_reference(true);
				reference = &Nodelist[final];
			}


			branchlist[i].Set_first_Node(&Nodelist[first]);
			branchlist[i].Set_final_Node(&Nodelist[final]);
			cout << "Has Resistor (" << i << ")? Yes or No :";
			res = new Resistor();
			cin >> answer;
			if (answer == "Yes") {
				cout << "Value:";
				cin >> value;
				res->Set_value(value);
				res->set_Is_or_not(true);
				branchlist[i].set_resistor(res);
			}
			else {
				res->Set_value(float(0.0005));
				res->set_Is_or_not(false);
				branchlist[i].set_resistor(res);
			}
			/////////////////////////////////////////////////////////
			cout << "Has voltage_source (" << i << ")? Yes or No :";
			Vsource = new Voltage_source();
			cin >> answer;
			if (answer == "Yes") {
				cout << "Value:";
				cin >> value;
				Vsource->Set_value(value);
				branchlist[i].set_voltage_source(Vsource);
			}
			else {
				Vsource->Set_value(0);
				branchlist[i].set_voltage_source(Vsource);
			}
			//////////////////////////////////////////////////////////
			cout << "Has curent source (" << i << ")? Yes or No :";
			cin >> answer;
			Csource = new Current_source();
			if (answer == "Yes") {
				cout << "Value :";
				cin >> value;
				Csource->Set_value(value);
				branchlist[i].set_curent_source(Csource);
			}
			else {
				Csource->Set_value(0);
				branchlist[i].set_curent_source(Csource);
			}
		}
		reference->set_reference(true);
		//cout<<reference->get_size();
		

	}
	void set_Incidence_matrix() {
		int i ,j,Node_index=0;
		Incidence_matrix = new int*[this->Number_of_Nodes-1];
		for (i = 0; i < this->Number_of_Nodes - 1; i++)
			Incidence_matrix[i] = new int[this->Number_of_Branchs];
		//////////////////////////////////////////////////////////////
       	for (i = 0; i < this->Number_of_Nodes ; i++) {
			if (Nodelist[i].is_reference() == true) {
			  //Node_index++;
				continue;
			}
			for (j = 0; j < this->Number_of_Branchs; j++) {
				if (branchlist[j].Get_first_Node()->Get_Code() == Nodelist[i].Get_Code())
					Incidence_matrix[Node_index][j] = 1;
				else 
					if (branchlist[j].Get_final_Node()->Get_Code() == Nodelist[i].Get_Code())
						Incidence_matrix[Node_index][j] = -1;
					else
						Incidence_matrix[Node_index][j] = 0;
				
			}
			Node_index++;
		}
		/*for (i = 0; i < this->Number_of_Nodes - 1; i++) {
			for (j = 0; j < this->Number_of_Branchs; j++) {
				cout << Incidence_matrix[i][j]<<" ";
				if (j == this->Number_of_Branchs - 1)
					cout << endl;
			}
		}
		cout << "////////////////////////////A"<<endl;*/

	}
	void set_conductunce_matrix() {
		int i,j;
		Conductunce_matrix = new float*[Number_of_Branchs];
		for (i = 0; i < Number_of_Branchs; i++)
			Conductunce_matrix[i] = new float[Number_of_Branchs];
		for (i = 0; i < Number_of_Branchs; i++) {
			for (j = 0; j < Number_of_Branchs; j++) {
				if (i == j)
					Conductunce_matrix[i][j] = 1 / branchlist[i].get_resistor()->Get_value();
				else
					Conductunce_matrix[i][j] = 0;
			}

		}
	   /*for (i = 0; i < Number_of_Branchs; i++) 
			for (j = 0; j < Number_of_Branchs; j++) {
				cout << Conductunce_matrix[i][j] << " ";
				if (j == this->Number_of_Branchs - 1)
					cout << endl;
			}
	   cout << "//////////////////////////////////G"<<endl;*/
	}
	void set_Incidence_transpose_matrix() {
		int i, j;
		Incidence_transpose_matrix = new int*[Number_of_Branchs];
		for (i = 0; i < Number_of_Branchs; i++)
			Incidence_transpose_matrix[i] = new int[Number_of_Nodes - 1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			for (j = 0; j < Number_of_Branchs; j++)
				Incidence_transpose_matrix[j][i] = Incidence_matrix[i][j];
		/*for (i = 0; i < this->Number_of_Branchs; i++) {
			for (j = 0; j < this->Number_of_Nodes - 1; j++) {
				cout << Incidence_transpose_matrix[i][j] << " ";
				if (j == this->Number_of_Nodes - 2)
					cout << endl;
			}

		}
		cout << "/////////////////////////////At"<<endl;*/
	}
	void set_Voltage_Source_of_Branches(){
		int i;
		Voltage_Source_of_Branches = new float[Number_of_Branchs];
		for (i = 0; i < Number_of_Branchs; i++)
			Voltage_Source_of_Branches[i] = branchlist[i].get_voltage_source()->Get_value();
		/*for (i = 0; i < Number_of_Branchs; i++)
			cout << Voltage_Source_of_Branches[i] << endl;
		cout << "///////////////////////////////////////Vs"<<endl;*/
	
	}
	void set_Current_Source_of_Branches() {
		int i;
	    Current_source_of_Branches = new float[Number_of_Branchs];
		for (i = 0; i < Number_of_Branchs; i++)
			Current_source_of_Branches[i] = branchlist[i].get_curent_source()->Get_value();
		/*	for (i = 0; i < Number_of_Branchs; i++)
		cout << Source_Voltage_of_Branches[i] << endl;*/

	}
	void calculate_IC_matrix() {
		int i, j, k;
		IC_matrix = new float*[Number_of_Nodes - 1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			IC_matrix[i] = new float[Number_of_Branchs];

		for (i = 0; i < Number_of_Nodes - 1; i++)
			for (j = 0; j < Number_of_Branchs; j++)
				IC_matrix[i][j] = 0;
/////////////////////////////////////////////////////////////////
		for (i = 0; i < Number_of_Nodes - 1; i++) {
			for (j = 0; j < Number_of_Branchs; j++) {
				for (k = 0; k < Number_of_Branchs; k++)
					IC_matrix[i][j] += Incidence_matrix[i][k] * Conductunce_matrix[k][j];
			}
		}
		/*for (i = 0; i < Number_of_Nodes - 1; i++) {
			for (j = 0; j < Number_of_Branchs; j++) {
				cout << IC_matrix[i][j]<<" "       ;
				if (j == Number_of_Branchs - 1)
					cout << endl;
			}
		}
		cout << "///////////////////////////AG" << endl;*/
	}
	void calculate_ICIt_matrix() {
		int i, j, k;
		ICIt_matrix = new float*[Number_of_Nodes - 1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			ICIt_matrix[i] = new float[Number_of_Nodes - 1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			for (j = 0; j < Number_of_Nodes - 1; j++)
				ICIt_matrix[i][j] = 0;
		for (i = 0; i < Number_of_Nodes - 1; i++) {
			for (j = 0; j < Number_of_Nodes - 1; j++) {
				for (k = 0; k < Number_of_Branchs; k++)
					ICIt_matrix[i][j] += IC_matrix[i][k] * Incidence_transpose_matrix[k][j];
			}
		}
		/*for (i = 0; i < Number_of_Nodes - 1; i++) {
		    for (j = 0; j < Number_of_Nodes - 1; j++) {
		         cout << ICIt_matrix[i][j]<<" ";
		    if (j == Number_of_Nodes - 2)
		         cout << endl;
		    }
		}
		cout << "//////////////////////////////////////AGAt" << endl;*/
	}
	void calculate_ICSv_matrix() {
		int i, k;
		ICSv_matrix = new float*[Number_of_Nodes - 1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			ICSv_matrix[i] = new float[1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
		    ICSv_matrix[i][0] = 0;
		for (i = 0; i < Number_of_Nodes - 1; i++)
			for (k = 0; k < Number_of_Branchs; k++)
				ICSv_matrix[i][0] += IC_matrix[i][k] * Voltage_Source_of_Branches[k];

	/*	for (i = 0; i < Number_of_Nodes - 1; i++){
			cout << ICSv_matrix[i][0] << " ";
		        cout << endl;

	    }
		cout << "/////////////////////////////AGVs"<<endl;*/

	/*	for (i = 0; i < Number_of_Branchs; i++) {
			cout << Voltage_Source_of_Branches[i] << " ";
			cout << endl;
		}*/
	}
	void calculate_ISc_matrix() {
		int i, k;
		ISc_matrix = new float*[Number_of_Nodes - 1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			ISc_matrix[i] = new float[1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			ISc_matrix[i][0] = 0;
		for (i = 0; i < Number_of_Nodes - 1; i++)
			for (k = 0; k < Number_of_Branchs; k++)
				ISc_matrix[i][0] += Incidence_matrix[i][k] * Current_source_of_Branches[k];

	  /*for (i = 0; i < Number_of_Nodes - 1; i++) {
			cout << ISc_matrix[i][0] << " ";
			cout << endl;
		}
		cout << "////////////////////////////////AJs"<<endl;*/
		/*for (i = 0; i < Number_of_Branchs; i++) {
			cout << Current_source_of_Branches[i] << " ";
			cout << endl;
		}*/


		}

	void calculate_Right_Equation_matrix() {
		int i;
		Right_Equation_matrix = new float*[Number_of_Nodes - 1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			Right_Equation_matrix[i] = new float[1];
		for (i = 0; i < Number_of_Nodes - 1; i++)
			Right_Equation_matrix[i][0] = ICSv_matrix[i][0] - ISc_matrix[i][0];
	/*	for (i = 0; i < Number_of_Nodes - 1; i++)
			cout << Right_Equation_matrix[i][0]<<endl;*/
	}

	void Solve() {
		unsigned char option(0);
		int row = 1, col = 1, index=0;
		row = Number_of_Nodes - 1;
		col = Number_of_Nodes - 1;
		register Array < double > *system = new Array < double >(row, col + 1);
		system->readEquation(ICIt_matrix, Right_Equation_matrix);
//		system->writeEquation();

		register Array < double > solve(InvertSolutions(*system));//solve(col, 1);
		char constant = 't';
//		cout << "Solutions:" << endl;
		for (int i = 0; i < system->Col() - 1; i++) {
//			cout << constant++ << "  =  " << solve[i][0] << endl;
			if (Nodelist[index].is_reference() == false){
				if(col!=1 && row!=1)
				    Nodelist[index].set_Voltage(float(solve[i][0]));
				else
					Nodelist[index].set_Voltage(Right_Equation_matrix[i][0]/ICIt_matrix[i][0]);

			}
			else {
				Nodelist[index].set_Voltage(0);
				i--;
	   		}
	 		index++;
		}
		

	//	for (int i = 0; i < Number_of_Nodes; i++)
		//	cout << Nodelist[i].get_Voltage()<<endl;
	}
	void calculate_current_and_Voltage() {
		int i;
		float value,v,r,vs;
		for (i = 0; i < Number_of_Branchs; i++) {
			branchlist[i].set_voltage(branchlist[i].Get_first_Node()->get_Voltage() - branchlist[i].Get_final_Node()->get_Voltage());
		}
		for (i = 0; i < Number_of_Branchs; i++) {
			if (branchlist[i].get_resistor()->get_Is_or_not() == true) {
				v = branchlist[i].get_voltage();
				r = branchlist[i].get_resistor()->Get_value();
		//		js = branchlist[i].get_curent_source()->Get_value();
				vs = branchlist[i].get_voltage_source()->Get_value();
				value = (v / r) - (vs / r);
				branchlist[i].get_resistor()->Set_current(value);
			}
		}
		for (i = 0; i < Number_of_Branchs; i++) {
			if (branchlist[i].get_resistor()->get_Is_or_not() == true) {
				value = branchlist[i].get_resistor()->Get_current();
				branchlist[i].get_resistor()->Set_Voltage(value*branchlist[i].get_resistor()->Get_value());
			}
		}
	}
	void print_result() {
		int i;
		cout << endl << "Result: "<<endl<<endl;
		for (i = 0; i < Number_of_Branchs; i++) {
			if (branchlist[i].get_resistor()->get_Is_or_not() == true) {
				cout  << "electrode(+): " << branchlist[i].Get_first_Node()->Get_Code() <<endl << "electrode(-): " << branchlist[i].Get_final_Node()->Get_Code() << endl;
				cout << "voltage(" << branchlist[i].get_resistor()->Get_value() << " Ohm) :" << branchlist[i].get_resistor()->Get_Voltage() << " Volt " << endl;
				cout << "current(" << branchlist[i].get_resistor()->Get_value() << " Ohm) :" << branchlist[i].get_resistor()->Get_current() << " Amper" << endl;
				cout << "////////////////////////////////////////////////"<<endl;
			}

		}


	}



	void Set_number_of_Node(int n) {
		this->Number_of_Nodes = n;
	}
	void Set_number_of_Branch(int n) {
		this->Number_of_Branchs = n;
	}
	int Get_number_of_Node() {
		return this->Number_of_Nodes;
	}
	int Get_number_of_Branch() {
		return this->Number_of_Branchs;
	}
	~Circuit() {

		delete[] IC_matrix;
		delete[] ICIt_matrix;
		delete[] ICSv_matrix;
		delete[] ISc_matrix;
		delete[] Right_Equation_matrix;
		delete[] solved_matrix;
		delete[] Conductunce_matrix;
		delete[] Incidence_transpose_matrix;
		delete[] Voltage_Source_of_Branches;
		delete[] Current_source_of_Branches;
		delete[] Incidence_matrix;
	//	delete[] reference;
		delete[] Nodelist;
		delete[] branchlist;
	}

};



int main()
{
	Circuit circuit;
	circuit.set_Incidence_matrix();
	circuit.set_conductunce_matrix();
	circuit.set_Incidence_transpose_matrix();
	circuit.set_Voltage_Source_of_Branches();
	circuit.set_Current_Source_of_Branches();
	circuit.calculate_IC_matrix();
	circuit.calculate_ICIt_matrix();
	circuit.calculate_ICSv_matrix();
	circuit.calculate_ISc_matrix();
	circuit.calculate_Right_Equation_matrix();
	circuit.Solve();
	circuit.calculate_current_and_Voltage();
	circuit.print_result();




    return 0;
}

/////////////////////////////////////
template <class type>  inline Array < type > InvertSolutions(Array < type > MainEquations) {

	Array < type > Oprands(MainEquations.Row(), MainEquations.Col() - 1);

	for (int i = 0; i < Oprands.Row(); i++)
		for (int j = 0; j < Oprands.Col(); j++)
			Oprands[i][j] = MainEquations[i][j];

	double det = Oprands.determinant();

	if (det) {
		Array < type > Left(MainEquations.Row(), 1);
		for (int i = 0; i < MainEquations.Row(); i++)
			Left[i][0] = MainEquations[i][MainEquations.Col() - 1];

		return Oprands.invers() * Left;
	}
}