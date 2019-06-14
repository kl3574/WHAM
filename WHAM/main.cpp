#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cctype>

const double unit = 120.2723625;
const double R = 8.314462145;

using namespace std;
bool isnum(string s);

int main(int argc, char** argv) {
	if (argc == 1) {
		cout << "No input command!\n";
		exit(-1);
	}
	string command = argv[1];
	// ./out t-wham Emin Emax Estep T0 Nlipid *.txt (T.txt)
	if (command == "t-wham") {
		string test;
		stringstream ss;
		test = argv[2];
		if (!isnum(test)) {
			cout << "Invalid input Emin: " << test << " !\n";
			exit(-1);
		}
		test = argv[3];
		if (!isnum(test)) {
			cout << "Invalid input Emax: " << test << " !\n";
			exit(-1);
		}
		test = argv[4];
		if (!isnum(test)) {
			cout << "Invalid input Estep: " << test << " !\n";
			exit(-1);
		}
		test = argv[5];
		if (!isnum(test)) {
			cout << "Invalid input reference Temperature: " << test << " !\n";
			exit(-1);
		}
		test = argv[6];
		if (!isnum(test)) {
			cout << "Invalid input Nlipid: " << test << " !\n";
			exit(-1);
		}
		ss << argv[2];
		double temp;
		ss >> temp;
		const double Emin = temp;
		ss.clear();
		ss << argv[3];
		ss >> temp;
		const double Emax = temp;
		ss.clear();
		ss << argv[4];
		ss >> temp;
		const double Estep = temp;
		ss.clear();
		ss << argv[5];
		ss >> temp;
		const double T0 = temp;
		ss.clear();
		ss << argv[6];
		ss >> temp;
		const int Nlipid = int(temp);
		ss.clear();
		if (Emin > Emax) {
			cout << "Emin should be lower than Emax!\n";
			exit(-1);
		}
		if (Emin + Estep > Emax) {
			cout << "Invalid Estep!\n";
			exit(-1);
		}
		const int Nbin = int((Emax - Emin) / Estep);
		const int Nrepl = argc - 7;
		//read T.txt
		ifstream fi2;
		fi2.open("T.txt");
		if (!fi2) {
			cout << "No Replica Temperature file: T.txt!\n";
			exit(-1);
		}
		double* T;
		T = new double[Nrepl] {};
		for (int i = 0; i < Nrepl; ++i) {
			fi2 >> T[i];
			cout << "Replica " << i << " : " << T[i] << "K\n";
		}
		fi2.close();
		//initialize hist
		int** hist;
		hist = new int*[Nbin];
		for (int i = 0; i < Nbin; ++i) {
			hist[i] = new int [Nrepl] {};
		}
		//write to twham.txt
		string suffix = argv[5];
		ofstream fo;
		fo.open("twham_" + suffix + ".txt");
		fo.setf(ios::fixed);
		fo.setf(ios::left);
		fo.precision(8);
		fo << setw(16) << "Potential+pV" << " Probability Density Function" << endl;
		//calculate Cp[J/mol/K] & Binder Cumulant write to E_cal.txt
		ofstream fo2;
		fo2.open("E_cal.txt", ios::app);
		fo2.setf(ios::left);
		//read enthalpy_*.txt
		//Total counts in simulation i
		int* N;
		N = new int[Nrepl] {};
		for (int j = 7; j < argc; ++j) {
			ifstream fi;
			string input = argv[j];
			fi.open(input);
			string temp = input.substr(input.find_last_of('_') + 1, input.find_first_of('.') - input.find_last_of('_') - 1);
			stringstream ss;
			int i;
			ss << temp;
			ss >> i;
			ss.clear();
			cout << "input file = " << input << ", index = " << i << endl;
			double data;
			while (fi >> data) {
				++N[i];
				int k = int((data - Emin) / Estep);
				++hist[k][i];
			}
			fi.close();
		}
		cout << endl;
		//wham iteration
		const double epsilon = 1e-14; //convergence criterion
		//calculate lnN[i]
		double* ln_N;
		ln_N = new double[Nrepl] {};
		for (int i = 0; i < Nrepl; ++i) {
			ln_N[i] = log(N[i]);
		}
		//initialize ln(c[i][k]) = -(bi-b0)*Ek
		double** ln_c;
		ln_c = new double*[Nrepl];
		for (int i = 0; i < Nrepl; ++i) {
			ln_c[i] = new double[Nbin] {};
			for (int k = 0; k < Nbin; ++k) {
				ln_c[i][k] = (1 / T0 - 1 / T[i]) * unit * (Emin + k * Estep);
			}
		}
		//initialize n[k]
		double* n;
		n = new double[Nbin] {};
		for (int k = 0; k < Nbin; ++k) {
			for (int i = 0; i < Nrepl; ++i) {
				n[k] += hist[k][i];
			}
		}
		//initialize ln_f[i]=0
		double* ln_f;
		ln_f = new double[Nrepl] {};
		for (int i = 0; i < Nrepl; ++i) {
			ln_f[i] = 0;
		}
		//initialize ln_f_old
		double* ln_f_old;
		ln_f_old = new double[Nrepl] {};
		//initialize ln_p
		double* ln_p;
		ln_p = new double[Nbin] {};
		double sum3;
		do {
			for (int i = 0; i < Nrepl; ++i) {
				ln_f_old[i] = ln_f[i];
			}
			for (int k = 0; k < Nbin; ++k) {
				double ln_sum1 = ln_N[0] + ln_f[0] + ln_c[0][k];
				for (int i = 1; i < Nrepl; ++i) {
					double temp = ln_N[i] + ln_f[i] + ln_c[i][k];
					if (ln_sum1 > temp)
						ln_sum1 = ln_sum1 + log(1 + exp(temp - ln_sum1));
					else
						ln_sum1 = temp + log(1 + exp(ln_sum1 - temp));
				}
				ln_p[k] = log(n[k]) - ln_sum1;
			}
			//norm p[k]
			double ln_sump = ln_p[0];
			for (int k = 1; k < Nbin; ++k) {
				if (ln_sump > ln_p[k])
					ln_sump = ln_sump + log(1 + exp(ln_p[k] - ln_sump));
				else
					ln_sump = ln_p[k] + log(1 + exp(ln_sump - ln_p[k]));
			}
			for (int k = 0; k < Nbin; ++k) {
				ln_p[k] -= ln_sump;
			}
			//update f[i]
			for (int i = 0; i < Nrepl; ++i) {
				double ln_sum2 = ln_c[i][0] + ln_p[0];
				for (int k = 1; k < Nbin; ++k) {
					double temp = ln_c[i][k] + ln_p[k];
					if (ln_sum2 > temp)
						ln_sum2 = ln_sum2 + log(1 + exp(temp - ln_sum2));
					else
						ln_sum2 = temp + log(1 + exp(ln_sum2 - temp));
				}
				ln_f[i] = -ln_sum2;
			}
			/*cout << "f_old is as follows,\n";
			for (int i = 0; i < Nrepl; ++i)
				cout << f_old[i] << " ";
			cout << endl << "f_new is as follows,\n";
			for (int i = 0; i < Nrepl; ++i)
				cout << f[i] << " ";
			cout << endl;*/
			sum3 = 0;
			for (int i = 0; i < Nrepl; ++i) {
				sum3 += (1 - ln_f[i] / ln_f_old[i]) * (1 - ln_f[i] / ln_f_old[i]);
				/*if (i == Nrepl - 1)
				cout << "difference = " << sum3 << endl;*/
			}
		} while (sum3 > epsilon);
		//norm p[k]
		double ln_sump = ln_p[0];
		for (int k = 1; k < Nbin; ++k) {
			if (ln_sump > ln_p[k])
				ln_sump = ln_sump + log(1 + exp(ln_p[k] - ln_sump));
			else
				ln_sump = ln_p[k] + log(1 + exp(ln_sump - ln_p[k]));
		}
		for (int k = 0; k < Nbin; ++k) {
			ln_p[k] -= ln_sump;
		}
		cout << "log of f_old is as follows, ";
		for (int i = 0; i < Nrepl; ++i)
			cout << ln_f_old[i] << " ";
		cout << endl << "log of f_new is as follows, ";
		for (int i = 0; i < Nrepl; ++i)
			cout << ln_f[i] << " ";
		cout << endl;
		cout << "difference = " << sum3 << endl << endl;
		for (int k = 0; k < Nbin; ++k) {
			fo << setw(16) << Emin + Estep * k << setw(16) << exp(ln_p[k]) / Estep << endl;
		}
		fo.close();
		//unit of Cp = kJ/mol/K
		double E1 = 0;
		double E2 = 0;
		double E4 = 0;
		for (int k = 0; k < Nbin; ++k) {
			double Ek = ((Emin + Estep * k) * 1000 + 1.5*Nlipid*R*T0);// unit J/mol
			E1 += exp(ln_p[k]) * Ek;
			E2 += exp(ln_p[k]) * Ek*Ek;
			E4 += exp(ln_p[k]) * Ek*Ek*Ek*Ek;
		}
		fo2 << setw(16) << T0 << setw(16) << E1 / 1e3 << setw(16) << E2 / 1e6 << setw(16) << E4 / 1e12 << setw(16) << (E2 - E1 * E1)*unit / T0 / T0 / 1e6 << setw(16) << 1 - E4 / 3 / E2 / E2 << endl;
		fo2.close();
		//delete hist
		for (int i = 0; i < Nbin; ++i) {
			delete[] hist[i];
		}
		//delete T, n, f, f_old, N
		delete[] T;
		delete[] n;
		delete[] ln_f;
		delete[] ln_f_old;
		delete[] N;
		//delete c
		for (int i = 0; i < Nrepl; ++i) {
			delete[] ln_c[i];
		}
	}
	// ./out wham Emin Emax Estep *.txt (T.txt)
	else if (command == "wham") {
	string test;
	stringstream ss;
	test = argv[2];
	if (!isnum(test)) {
		cout << "Invalid input Emin: " << test << " !\n";
		exit(-1);
	}
	test = argv[3];
	if (!isnum(test)) {
		cout << "Invalid input Emax: " << test << " !\n";
		exit(-1);
	}
	test = argv[4];
	if (!isnum(test)) {
		cout << "Invalid input Estep: " << test << " !\n";
		exit(-1);
	}
	ss << argv[2];
	double temp;
	ss >> temp;
	const double Emin = temp;
	ss.clear();
	ss << argv[3];
	ss >> temp;
	const double Emax = temp;
	ss.clear();
	ss << argv[4];
	ss >> temp;
	const double Estep = temp;
	ss.clear();
	if (Emin > Emax) {
		cout << "Emin should be lower than Emax!\n";
		exit(-1);
	}
	if (Emin + Estep > Emax) {
		cout << "Invalid Estep!\n";
		exit(-1);
	}
	const int Nbin = int((Emax - Emin) / Estep);
	const int Nrepl = argc - 5;
	//read T.txt
	ifstream fi2;
	fi2.open("T.txt");
	if (!fi2) {
		cout << "No Replica Temperature file: T.txt!\n";
		exit(-1);
	}
	double* T;
	T = new double[Nrepl] {};
	for (int i = 0; i < Nrepl; ++i) {
		fi2 >> T[i];
		cout << "Replica " << i << " : " << T[i] << "K\n";
	}
	fi2.close();
	int** hist;
	hist = new int*[Nbin];
	for (int i = 0; i < Nbin; ++i) {
		hist[i] = new int [Nrepl] {};
	}
	//write to state_density.txt
	ofstream fo;
	fo.open("state_density.txt");
	fo.setf(ios::scientific, ios::uppercase);
	fo.setf(ios::left);
	fo.precision(8);
	fo << setw(16) << "Potential+pV" << "log of Probability Density Function" << endl;
	//read enthalpy_*.txt
	
	//Total counts N in simulation i
	int* N;
	N = new int[Nrepl] {};
	for (int j = 5; j < argc; ++j) {
		ifstream fi;
		string input = argv[j];
		fi.open(input);
		string temp = input.substr(input.find_last_of('_') + 1, input.find_first_of('.') - input.find_last_of('_') - 1);
		stringstream ss;
		int i;
		ss << temp;
		ss >> i;
		ss.clear();
		cout << "input file = " << input << ", index = " << i << endl;
		double data;
		while (fi >> data) {
			++N[i];
			int k = int((data - Emin) / Estep);
			++hist[k][i];
		}
		fi.close();
	}
	cout << endl;
	//wham iteration
	const double epsilon = 0; //convergence criterion
	//calculate lnN[i]
	double* ln_N;
	ln_N = new double[Nrepl] {};
	for (int i = 0; i < Nrepl; ++i) {
		ln_N[i] = log(N[i]);
	}
	//initialize ln(c[i][k]) = -(bi-b0)*Ek
	double** ln_c;
	ln_c = new double*[Nrepl];
	for (int i = 0; i < Nrepl; ++i) {
		ln_c[i] = new double[Nbin] {};
		for (int k = 0; k < Nbin; ++k) {
			ln_c[i][k] = -1 / T[i] * unit * (Emin + k * Estep);
		}
	}
	//initialize n[k] = sum_i (hist[k][i])
	double* n;
	n = new double[Nbin] {};
	for (int k = 0; k < Nbin; ++k) {
		for (int i = 0; i < Nrepl; ++i) {
			n[k] += hist[k][i];
		}
	}
	//initialize f[i]=1
	double* f;
	f = new double[Nrepl] {};
	for (int i = 0; i < Nrepl; ++i) {
		f[i] = 1;
	}
	//initialize f_old
	double* f_old;
	f_old = new double[Nrepl] {};
	//initialize ln_p
	double* ln_p;
	ln_p = new double[Nbin] {};
	double sum3;
	do {
		for (int i = 0; i < Nrepl; ++i) {
			f_old[i] = f[i];
		}
		for (int k = 0; k < Nbin; ++k) {
			double ln_sum1 = ln_N[0] + f[0] + ln_c[0][k];
			for (int i = 1; i < Nrepl; ++i) {
				double temp = ln_N[i] + f[i] + ln_c[i][k];
				if (ln_sum1 > temp)
					ln_sum1 = ln_sum1 + log(1 + exp(temp - ln_sum1));
				else
					ln_sum1 = temp + log(1 + exp(ln_sum1 - temp));
			}
			ln_p[k] = log(n[k]) - ln_sum1;
		}
		//norm p[k]
		double ln_sump = ln_p[0];
		for (int k = 1; k < Nbin; ++k) {
			if (ln_sump > ln_p[k])
				ln_sump = ln_sump + log(1 + exp(ln_p[k] - ln_sump));
			else
				ln_sump = ln_p[k] + log(1 + exp(ln_sump - ln_p[k]));
		}
		for (int k = 0; k < Nbin; ++k) {
			ln_p[k] -= ln_sump;
		}
		//update f[i]
		for (int i = 0; i < Nrepl; ++i) {
			double ln_sum2 = ln_c[i][0] + ln_p[0];
			for (int k = 1; k < Nbin; ++k) {
				double temp = ln_c[i][k] + ln_p[k];
				if (ln_sum2 > temp)
					ln_sum2 = ln_sum2 + log(1 + exp(temp - ln_sum2));
				else
					ln_sum2 = temp + log(1 + exp(ln_sum2 - temp));
			}
			f[i] = -ln_sum2;
		}
		/*cout << "f_old is as follows,\n";
		for (int i = 0; i < Nrepl; ++i)
			cout << f_old[i] << " ";
		cout << endl << "f_new is as follows,\n";
		for (int i = 0; i < Nrepl; ++i)
			cout << f[i] << " ";
		cout << endl;*/
		sum3 = 0;
		for (int i = 0; i < Nrepl; ++i) {
			sum3 += (1 - f[i] / f_old[i]) * (1 - f[i] / f_old[i]);
			/*if (i == Nrepl - 1)
				cout << "difference = " << sum3 << endl;*/
		}
	} while (sum3 > epsilon);
	//norm p[k]
	double ln_sump = ln_p[0];
	for (int k = 1; k < Nbin; ++k) {
		if (ln_sump > ln_p[k])
			ln_sump = ln_sump + log(1 + exp(ln_p[k] - ln_sump));
		else
			ln_sump = ln_p[k] + log(1 + exp(ln_sump - ln_p[k]));
	}
	for (int k = 0; k < Nbin; ++k) {
		ln_p[k] -= ln_sump;
	}
	cout << "f_old is as follows, ";
	for (int i = 0; i < Nrepl; ++i)
		cout << f_old[i] << " ";
	cout << endl << "f_new is as follows, ";
	for (int i = 0; i < Nrepl; ++i)
		cout << f[i] << " ";
	cout << endl;
	cout << "difference = " << sum3 << endl << endl;

	for (int k = 0; k < Nbin; ++k) {
		fo << setw(16) << Emin + Estep * k << setw(16) << (ln_p[k]) - log(Estep) << endl;
	}
	fo.close();
}
	//invalid command!
	else {
		cout << "Invalid command!\n";
		cout << "WHAM.out t-wham Emin Emax Estep T0 Nlipid *.txt (T.txt)\n";
		cout << "WHAM.out wham Emin Emax Estep *.txt (T.txt)\n";
		exit(-1);
	}
	return 0;
}

bool isnum(string s) {
	stringstream ss(s);
	double n;
	char c;
	if (!(ss >> n))
		return false;
	if (ss >> c)
		return false;
	else
		return true;
}