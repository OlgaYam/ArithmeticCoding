// ArithmeticCoding.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <ctime>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cstdio>
#include<map>

using namespace std;

const int size = 10000; //длина выборки
const int sizeBlock = 30;

string alphabet;
vector<double> probability; //вероятность
vector<double> Q_probability; //кумулятивная вероятность
//vector<bool> Code; //закодированное число
map<int, pair<double, double>> boardSym;
vector<int> lenBinaryBloks;
vector<int> lenBlocks;
map<char,vector<double>> matrix;



void createQ_Prob(vector<double> pr)
{
	Q_probability.push_back(0.0);
	for (int i = 1; i < alphabet.length(); ++i)
	{
		Q_probability.push_back(Q_probability[i-1]+pr[i-1]);
	}
}

//генерирование вероятностей
void GenerateProbab()
{
	cout << "Probabilities:" << endl;
	srand(time(NULL));
	double sum = 0.0;
	for(int j = 0; j != alphabet.length()-1; ++j)
	{
		double r = (rand()%100)/(100 * 1.0)*(1-sum);
		sum += r;
		probability.push_back(r);
		cout << r << " ";
	}
	double r = 1- sum;
	probability.push_back(r);
	cout << r << endl;

	createQ_Prob(probability);
}

//генерируем дискретный источник
int generateDiscrete(string source)
{
	vector<double> thresh;
	int len = 0;
	ofstream fileSource(source, ios::out | ios::binary);
	srand(time(NULL));

	for (int i = 0; i < size; ++i)
    {
		// генерируем рандомное число
		float r = (rand()%100)/(100 * 1.0);
		for (int j = probability.size()-1; j >= 0; j--) 
		{
			if (Q_probability[j] <= r) {
				fileSource << alphabet[j]; //запись в файл SourceCode
				len++;
				break;
			}
		}
	}
	fileSource.close();
	return len;
}


//генерируем текст с помощью марковского источника
void GenerateMarkov(map<char, double> init)
{
	//находим первый символ по вектору инициализации
	char sym;
	srand(time(NULL));
	double r = (rand()%100)/(100 * 1.0);
	double sum = 0.0;
	for (auto it = init.begin(); it != init.end(); ++it)
	{
		sum += it->second;
		if (r < sum)
			sym = it ->first;
	}
	ofstream fileSource("SourceCode.txt", ios::out | ios::binary);
	fileSource << sym; //записываем первый символ в файл
	for (int i = 1; i < size; ++i)
    {
		// генерируем рандомное число
		float r = (rand()%100)/(100 * 1.0);
		//определяем в каком состоянии находимся сейчас
		for (auto it = matrix.begin(); it != matrix.end(); ++it) 
			if (sym == it->first) {//попали сюда -> значит, нашли нужное состояние 
				float m = 0;
				vector<double> pr = it->second;
				//определяем, куда будем переходить дальше
				for (int j = 0; j < pr.size(); ++j) {
					m += pr[j];
					if (r < m) {
						sym = alphabet[j]; 
						break;
					}
				}
				break;
			}
			fileSource << sym; //запись в файл SourceCode
	}
	fileSource.close();
}

//генерация матрицы переходных вероятностей рандомно
void GenerationRandom()
{	
	cout << "\nGenerate random" << endl;
	srand(time(NULL));
	vector<double> v;
	for(int i = 0; i != alphabet.length(); ++i)
	{
		double sum = 0.0;
		for(int j = 0; j != alphabet.length()-1; ++j)
		{
			double r = (rand()%100)/(100 * 1.0)*(1-sum);
			sum += r;
			v.push_back(r);
			cout << r << " ";
		}
		double r = 1- sum;
		v.push_back(r);
		cout << r << endl;
		matrix.insert(pair<char, vector<double> > (alphabet[i], v));
		v.clear();
	}
}


//считываем матрицу переходных вероятностей из текстового файла
void ReadFromFile()
{
	float str = 0;
	string matrfile;
	cout << " Enter filename: ";
	cin >> matrfile;
	ifstream matr(matrfile, ios::out | ios::binary);
	if (matr.fail()) 
	{
		cout<< "\n Ошибка открытия файла";
		exit(1);
	}
	int i = 0; int k = 0;
	vector<double> v;
	while (!matr.eof())
	{ 
		matr >> str;
		v.push_back(str);
		i++;
		if (i == alphabet.length()) {
			matrix.insert(pair<char, vector<double> > (alphabet[k], v));
			v.clear();
			i = 0;
			k++;
		}
	}
}


//энтропия практическая
void PracticalEntropy(int countBytes, int len)
{
	double entropy2 = double(countBytes*8)/double(len);
    cout << "Practical entropy = " << entropy2 << endl;
}

//Энтропия теоритическая (дискретный источник)
void TeoreticalEntropy(vector<double> pr)
{
	double entropy1 = 0;
	for (int k = 0; k < pr.size(); ++k )
	{
		if (pr[k] != 0.0) {
			entropy1 += pr[k]*(log(pr[k])/log(2));
		}
	}
	entropy1 *= -1;
	cout << "Teoretical entropy = " << entropy1 << endl;
}


//поиск символа в алфавите
int findFromAlphabet(char c)
{
	for (int i = 0; i != alphabet.length(); ++i)
		if (alphabet[i] == c)
			return i;
	return -1;
}


//создаем вектор инициализации
map<char, double> createVecInit()
{
	vector<double> pr;
	map<char, double> init;
	srand(time(NULL));
	double sum = 0.0;
	for(int i = 0; i != alphabet.length()-1; ++i)
	{
		double r = (rand()%100)/(100 * 1.0)*(1-sum);
		sum += r;
		init.insert(pair<char,double>(alphabet[i], r) );
	}
	double r = 1- sum;
	init.insert(pair<char,double>(alphabet[alphabet.length()-1], r) );
	return init;
}


//перевод в двоичное число
vector<bool> ToBits(double word, int n)
{
	vector<bool> bits;
	for (int k = 0; k != n; ++k)
	{
		word *= 2;
		if (word >= 1) { 
			bits.push_back(1);
			word -= 1;
		}
		else  {
			bits.push_back(0);
		}
	}
	return bits;
}


void Coder(string source, string binary, vector<double> pr)
{
	//cout << "\nCoding:" << endl;
	Q_probability.clear();
	createQ_Prob(pr);
	ifstream fileIn(source, ios::in | ios::binary);
	ofstream fileBinary(binary, ios::out | ios::binary);
	string str; int len = 0;
	int countBytes = 0;
	while (!fileIn.eof())
    {
		fileIn >> str; //считываем из файла SourceCode
		for (int i = 0; i < str.length(); i+=sizeBlock)
		{
			double G = 1.0;
			double F = 0.0;
			string data = str.substr(i, sizeBlock);
			lenBlocks.push_back(data.size());
			vector<int> masdata;
			for (int j = 0; j != data.size(); ++j) 
				masdata.push_back(findFromAlphabet(data[j]));
			for (int j = 0; j < masdata.size(); ++j) {
				F += Q_probability[masdata[j]]*G;
				G *= pr[masdata[j]];
			}
			double value = F + G / 2; 
			//int n = 50;
			int n = ceil( (-1) * log(G)/log(2) + 1) + 5;
			lenBinaryBloks.push_back(n);
			vector<bool> bits = ToBits(value, n);
			/*cout << "Data = " << data << "\nVal = " << value << "\nBits = ";
			int k = 0;
			for (auto it = bits.begin() ; it!=bits.end() ; ++it) {
				if (k == 4) {
					cout << " "; 
					k = 0;
				}
                cout << *it;
				k++;
			}
			cout << endl;*/
			masdata.clear();
			
			//записываем код побитово
			int count=0; char buf=0;
			for(int k = 0; k < bits.size(); k++)
			{
				buf = buf | bits[k]<<(7-count);   
				count++;   
				//если образовался байт, то записываем его в файл
				if (count==8) { 
					count=0;   
					fileBinary<<buf; 
					buf=0; 
					countBytes++;
				} 
			}
			//дбиваем нулями
			while (count < 8)
			{
				buf = buf | 0 << (7-count);
				count++;
			}
			fileBinary<<buf; 
			countBytes++;
		}
		len += str.length();
	}

	fileIn.close();
	fileBinary.close();
	//энтропия практическая
	double entropy2 = double(countBytes*8)/double(len);
    cout << "\nPractical entropy = " << entropy2 << endl;
}


//перевод в десятичное число
double fromBits(string bits)
{
    double value = 0;
	for(int i = 0; i != bits.length(); i++)
    {
        value += bits[i] * pow(2,(-(i+1)));
    }
    return value;
}


void Decoder(string binary, string fOut, vector<double> pr)
{
	//cout << "Decoding:" << endl;
	Q_probability.push_back(1.0);

	// считывание из файла output.txt и преобразование обратно в биты
	string bits;
	ifstream binaryfile(binary, ios::in | ios::binary);
	ofstream fileOut(fOut, ios::out | ios::binary);

	int count=0;
	char byte = binaryfile.get();
	while(!binaryfile.eof())
	{   
		bool b = byte & (1 << (7-count) ) ; 
		bits += b;
		count++;
		if (count==8) {
			count=0; 
			byte = binaryfile.get();
		}
	}
	int len = 0;
	for(int k = 0; k < lenBinaryBloks.size(); ++k)
	{
		double S = 0.0;
		double G = 1.0;
		double F = 0.0;
		int l = len;
		string bit = bits.substr(len, lenBinaryBloks[k]);
		double val = fromBits(bit);
		/*cout << "\nBits = ";
		int g = 0;
		for (int it = 0; it < bit.length() ; ++it)
		{
			if (g == 4) {cout << " "; g = 0;}
			if (bit[it] == 0)
                cout << 0;
			else cout << 1;
			g++;
		}
		cout << "\nVal = " << val;
		cout << "\nData = "; */
		for (int i = 0; i < lenBlocks[k]; ++i)
		{
			int j = 0;
			double ss = S+Q_probability[j+1]*G;
			while (S + Q_probability[j+1]*G < val) 
			{
				ss = S+Q_probability[j+1]*G;
				if (j < Q_probability.size()-2) j++;
				else 
					break;
				if (j == 2) 
					int h = 0;
			}
			S = S + Q_probability[j]*G;
			G = pr[j]*G;
			fileOut << alphabet[j];
			//cout << alphabet[j];
		}
		len += lenBinaryBloks[k] + (8 - (lenBinaryBloks[k] % 8));
	}
	binaryfile.close();
	fileOut.close();
	lenBinaryBloks.clear();
	lenBlocks.clear();
}


//нахождение распредения символов по тексту
vector<double> FindDistr(string SourceFile)
{
	//распределение встречающихся символов
	vector<double> m; 
	//заполняем распределение нулями, так пока символы не встречались
	for(int i = 0; i != alphabet.length(); ++i)
		m.push_back(0.0);

	int lenText = 0;
	ifstream fileSource(SourceFile, ios::in | ios::binary);
	string str;
	while (!fileSource.eof())
	{ 
		fileSource >> str; //считали симовл
		for (int i = 0; i != str.length(); ++i)
		{
			int k = findFromAlphabet(str[i]);
			++m[k];
		}
		lenText += str.length();
	}
	for (int i = 0; i != m.size(); ++i) 
	{
		m[i] = m[i]/lenText;
	}
	fileSource.close();
	return m;
}


//для нахождения стационарного распределения вероятностей запишем матрицу в файл и вычислим на матлабе
void WriteMartInFile(map<char,vector<double>> matr1, string filename)
{
	double str = 0;
	ofstream matr(filename, ios::out | ios::binary);
	for (auto it = matr1.begin(); it != matr1.end(); ++it) 
    {
		for(int j = 0; j < it->second.size(); ++j)
        {
			matr << it->second[j] << " ";
		}
		matr << endl;
	}
	matr.close();
}


//чтение вероятностей из файла
vector<double> pReadFromFile()
{
	vector<double> m;
	double str = 0;
	ifstream matr("Probabilities.txt", ios::out | ios::binary);
	while (!matr.eof())
	{ 
		matr >> str;
		m.push_back(str);
	}
	return m;
}


//пересчет матрицы вероятностей по тексту
vector<double> CreateMatrOnText()
{
	//распределение встречающихся символов
	map<char,vector<double> > m;
	vector<double> v;
	for(int i = 0; i != alphabet.length(); ++i) v.push_back(0.0);
	//заполняем распределение нулями, так пока символы не встречались
	for(int i = 0; i != alphabet.length(); ++i)
		m.insert ( pair<char,vector<double>>(alphabet[i], v) );
	int lenText = 0;
	ifstream fileSource("SourceCode.txt", ios::out | ios::binary);
	string str; int i = 1;
	char c1;
	char c2;
	bool flag = true;
	while (!fileSource.eof())
	{ 
		fileSource >> str; //считали симовл
		if (flag) { c1 = str[0]; flag = false;}
		while (i != str.length())
		{
			c2 = str[i];
			int k = findFromAlphabet(c2);
			if (k != -1) ++m.find(c1)->second[k];
			else cout << "error in Coder4" << endl;
			++i;
			c1 = c2;
		}
		i = 0;
	}
	fileSource.close();

	for (auto it = m.begin(); it != m.end(); ++it)
	{
		vector<double> pr = it->second;
		double sum = 0;
		for (int j = 0; j < pr.size(); ++j) {
			sum += pr[j];
		}
		if (sum != 0.0) {
			for (int j = 0; j < pr.size(); ++j) {
				pr[j] = pr[j]/sum;
			}
		}
		it->second = pr;
	}
	WriteMartInFile(m, "matrixDistr.dat");
	cout << "Press 1 after the execution of code on Matlab" << endl;
	int g;
	cin >> g;
	return pReadFromFile();
}


//вывод вероятностей на экран
void writePr(vector<double> pr)
{
	cout << "Probabilities:" << endl;
	for (int i = 0; i < pr.size(); ++i)
		cout << pr[i] << " ";
	cout << endl;
}


int main (int argc, char *argv[])
{
		
	cout << "Enter alphabet (without spaces, one line)" << endl;
	cin >> alphabet;
	
	cout << "Discrete Source" << endl;
	GenerateProbab();
	int len = generateDiscrete("Source.txt");

	//Сжатие выхода дискретного стационарного источника без памяти, используя истинное распределение вероятностей
	cout << "\nCoding, using the original probability distribution" << endl;
	Coder("Source.txt", "Binary1.txt", probability);
	Decoder("Binary1.txt", "Decode1.txt", probability);
	TeoreticalEntropy(probability);


	//Сжатие выхода дискретного стационарного источника без памяти, оцененного по потоку распределение вероятностей
	cout << "\nCoding for the evaluation of the probability distribution" << endl;
	vector<double> p = FindDistr("Source.txt");
	writePr(p);
	Coder("Source.txt", "Binary2.txt", p);
	Decoder("Binary2.txt", "Decode2.txt", p);
	TeoreticalEntropy(p);

	//Сжатие выхода дискретного стационарного источника без памяти, при использовании равных вероятностей	
	cout << "\nCoding, using the wrong probability distribution" << endl;
	vector<double> p1;
	for(int j = 0; j != alphabet.length(); ++j)
	{
		double r = 1.0/alphabet.length();
		p1.push_back(r);
	}
	writePr(p1);
	Coder("Source.txt", "Binary3.txt", p1);
	Decoder("Binary3.txt", "Decode3.txt", p1);
	TeoreticalEntropy(p1);


	
	//alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 ";
	
	cout << "\nMarkov's source" << endl;

	cout << " Select the method of specifying the matrix:" << endl;
	cout << " 1 - generate random" << endl;
	cout << " 2 - from file *.txt" << endl;
	int k;
	cin >> k;
	switch ( k )
	{
	case 1: 
		GenerationRandom(); break;
	case 2:
		ReadFromFile(); break;
	default: 
		cout << " Wrong operation" << endl; 
		system("pause");
		return 0;
	}

	map<char, double> init = createVecInit();
	string Sourcefile = "SourceCode.txt";
	GenerateMarkov(init); 	
	
	
	cout << "\nCoding for the evaluation of the probability distribution" << endl;
	vector<double> prMark1 = FindDistr("SourceCode.txt");
	writePr(prMark1);
	Coder(Sourcefile, "Binary4.txt", prMark1);
	Decoder("Binary4.txt", "Decode4.txt", prMark1);
	TeoreticalEntropy(prMark1);
	

	cout << "\n Coding, assessing the transition probability matrix for power output" << endl;
	vector<double> prMark2 = CreateMatrOnText();
	writePr(prMark2);
	Coder(Sourcefile, "Binary5.txt", prMark2);
	Decoder("Binary5.txt", "Decode5.txt", prMark2);


	cout << "\n Coding, using the original matrix" << endl;
	WriteMartInFile(matrix, "matrixDistr.dat");
	cout << "Press 1 after the execution of code on Matlab" << endl;
	int m;
	cin >> m;
	vector<double> prMark3 = pReadFromFile();
	writePr(prMark3);
	Coder(Sourcefile, "Binary6.txt", prMark3);
	Decoder("Binary6.txt", "Decode6.txt", prMark3);

	cout << "Decoding is complete" << endl;
	system("pause");
	return 0;
}