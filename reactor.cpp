#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#define GROUPS 26
#define Pi 3.1415926536

struct RESULT
{
	int Day;
	double Concentration[100];
	std::string Name_Concentration[100];
	double K_effective, betta_neutron, lifetime_neutron, breeding_ratio, sum_flux;
	double flux[26];
	double SIGMA_F[26];
	double SIGMA_C[26];
	double SIGMA_A[26];
};

class table_of_element //Класс описывающий таблицу групповых констант для конкретного элемента
{
	public: // Здесь название элемента и константы
		std::string name_of_element;
		double sigma_tot[GROUPS];
		double sigma_f[GROUPS];
		double nu[GROUPS];
		double sigma_c[GROUPS];
		double sigma_in[GROUPS];
		double sigma_e[GROUPS];
		double mu_e[GROUPS];
		double epsilon[GROUPS];
		double sigma_3e[GROUPS];
		double mu_3e[GROUPS];
		double sigma_a[GROUPS];
		double sigma_tr[GROUPS];
		double sigma_in26[100][100];
		double sigma_in26_ik[100][100];
		double sigma_e26_ik[100][100];
		double sigma_e_plus_i_26_ik[100][100];//временная переменная
	
		double table_f[7][100];
		double table_c[7][100];
	
		double sigma_0[26];
		
		double const_sigma_c[26];
		double const_sigma_f[26];
		
		double betta[26];
		
		double betta_i_betta[6];
		double re_labda[6];
		double lambda[6];
		double L_delay=0;
		
		table_of_element(char *path_file, char *name)//конструктор принимает путь к табице элемента и название самого элемента
		{
			for(int i = 0; i < 7;i++)
			{
				for(int j = 0; j < 100; j++)
				{
					table_f[i][j]=1;
					table_c[i][j]=1;
				}
			}
			
			for(int i = 0; i < 26; i++)
			{
				sigma_0[i]=0;
			}
			
			char pf[100]="elements\\";
			strcat(pf,path_file);
			double sigma_in26_temp[27][25];//переменная для таблицы сечений переходов i i+k ein
			name_of_element = name;
			std::ifstream file(pf);
			std::string buff;
			if(!file)
			std::cout << "file: " << path_file << " is not opened" << "\n";
			file >> buff;//пропускаем лишнюю строку;
			buff.clear();
			for(int p = 0; p <GROUPS; p++)
			{//в цикле переносятся константы из файла в переменные из private
				file >> sigma_tot[p] >> sigma_f[p]>> nu[p]>> 
				sigma_c[p] >> sigma_in[p] >> sigma_e[p]>> mu_e[p]>>
				epsilon[p] >> sigma_3e[p] >> mu_3e[p];
			}
			
			for(int r = 0; r < GROUPS; r++)//сечения поглащения и транспортное считаются отдельно
			{
				sigma_a[r] = sigma_f[r] + sigma_c[r];
				sigma_tr[r] = sigma_a[r] + sigma_in[r] + sigma_e[r]*(1-mu_e[r]);
			}
			//пропускаем строку
			file >> buff;
			buff.clear();
			//далее заполняем таблицу неупругих сечений вида i,i+k в массив sigma_in26_temp
			for(int i=0;i<25;i++)
			{
				for(int p=0;p<27;p++)
				{
					file >> sigma_in26[p][i];
				}
			}
			
			//Заполняем table_c и table_f 1ами(на всякий случай)
			for(int y = 0; y < 100; y++)
			{
				for(int x = 0; x < 7; x++)
				{
					table_c[x][y]=1;
					table_f[x][y]=1;
				}
			}
			
			
			//Заполняем table_c
			//пропускаем строку
			file >> buff;
			buff.clear();
			for(int y = 0; y < 100; y++)
			{
				for(int x = 0; x < 7; x++)
				{
					file >> table_c[x][y];
				}
			}
			//Заполняем table_f
			//пропускаем строку
			file >> buff;
			buff.clear();
			for(int y = 0; y < 100; y++)
			{
				for(int x = 0; x < 7; x++)
				{
					file >> table_f[x][y];
				}
			}
			
			//заполняем долю нейтронов зап. betta
			file >> buff;
			for(int y = 0; y < 26; y++)
			{
				file >> betta[y];
				//std::cout << get_name() << " " << betta[y] << "\n";
			}			
			
			//заполняем lambda
			file >> buff;
			for(int y = 0; y < 6; y++)
			{
				file >> lambda[y];
				//std::cout << get_name() << " " << betta[y] << "\n";
			}	
			
			//Получаем 1/lamda
			for(int y = 0; y < 6; y++)
			{
				if(lambda[y]!=0)
				re_labda[y]= (1/lambda[y]);
				else
				re_labda[y]=0;
				//std::cout << get_name() << " " << re_labda[y] << "\n";
			}
			//заполняем betta_i_betta
			file >> buff;
			for(int y = 0; y < 6; y++)
			{
				file >> betta_i_betta[y];
				//std::cout << get_name() << " " << betta[y] << "\n";
			}	
			
			//Ищем L_delay
			for(int y=0;y<6;y++)
			{
				L_delay+=betta_i_betta[y]*re_labda[y];
				//std::cout << get_name() << " " << L_delay << "\n";
			}
			
			//заполняем массив нулями
			for(int i=0;i<100;i++)
			{
				for(int m=0;m<100;m++)
				{
					sigma_in26_ik[i][m]=0;
				}
			}
			
			//преобразуем sigma_in26_temp в вид (i,k)
			int t=70;
			for(int u=0;u<t;u++)
			{
				for(int i=0;i<26;i++)
				{
					sigma_in26_ik[i+u][u] = sigma_in26[i][u];
				}
			}
			//массив нулями
			for(int i=0;i<100;i++)
			{
				for(int m=0;m<100;m++)
				{
					sigma_e26_ik[i][m]=0;
				}
			}
			//заполняем таблицу sigma_e;
			for(int u=1;u<90;u++)
			{
				sigma_e26_ik[u][u-1]=sigma_3e[u-1];
			}
			//заполняем таблицу суммарного сечения рассеивания
			for(int i=0;i<26;i++)
			{
				for(int m=0;m<26;m++)
				{
					sigma_e_plus_i_26_ik[i][m]=sigma_e26_ik[i][m]+sigma_in26_ik[i][m];
					//std::cout << sigma_e_plus_i_26_ik[m][i] << "\t";
				}
				//std::cout << "\n";
			}
			
			for(int i = 0; i < 26; i++)
			{
				const_sigma_c[i]=sigma_c[i];
				const_sigma_f[i]=sigma_f[i];
			}
			
		}
		
		std::string get_name()
		{
			return name_of_element;
		}
		
		double get_sigma_tot(int group){return sigma_tot[group];}
		double get_sigma_f(int group){return sigma_f[group];}
		double get_sigma_c(int group){return sigma_c[group];}
		double get_nu(int group){return nu[group];}
		double get_sigma_in(int group){return sigma_in[group];}
		double get_sigma_e(int group){return sigma_e[group];}
		double get_mu_e(int group){return mu_e[group];}
		double get_epsilon(int group){return epsilon[group];}
		double get_sigma_3e(int group){return sigma_3e[group];}
		double get_mu_3e(int group){return mu_3e[group];}
		double get_sigma_a(int group){return sigma_a[group];}
		double get_sigma_tr(int group){return sigma_tr[group];}
		double get_sigma_e_plus_i_26_ik(int x, int y){return sigma_e_plus_i_26_ik[x][y];}
		double get_L_delay(){return L_delay;};

};

double inter_mod(double data[7][100], double xn, double yn, int group)
{
	double x = xn;
	double y = yn;
	double f;
	
	//ищем координаты 4 промежуточных значений
	int x1,x2,y1,y2;
	
	if(data[1][0]==1)
		return 1;
	
	if(y==0)
	{
		return 1;
	}
	
	if(x==0)
	{
		x+=0.00001;
	}
	
	for(int i = 1; i < 7; i++)
	{
		if(x>data[i][0])
		{
			x2 = i;
			x1 = i-1;
			//std::cout << x2 << "  " << x1;
			break;
		}
	}
	
	for(int i = 4*group-3; i < 4*group; i++)
	{
		if(y<data[0][i])
		{
			y2 = i;
			y1 = i-1;
			//std::cout << y2 << "  " << y1;
			break;
		}
	}
	
	//ищем промежуточные 2 значения по x;
	double f1, f2;
	f1 = data[x2][y1]+(data[x1][y1]-data[x2][y1])/(data[x1][0]-data[x2][0])*(x-data[x2][0]);
	f2 = data[x2][y2]+(data[x1][y2]-data[x2][y2])/(data[x1][0]-data[x2][0])*(x-data[x2][0]);
	
	//интерполируем по y;
	
	double inter;
	
	inter = f1+(f2-f1)/(data[0][y2]-data[0][y1])*(y-data[0][y1]);
	
	if(inter>1||inter<0)
		return 1;
	
	return inter;	
}

double sum_of_string(double b[100][100], int u)
{
	double r=0;
	for(int i = 0; i<26; i++)
	{
		r += b[i][u];	
	}
	return r;
}

double sum_of_multiplying(int x, double n[100][100], double m[26], int k)//функция сумма произведения sigma_R и F0;
{//x-номер столбца в sigma_R; n - sigma_R; m - F0; k - номер группы 
	double b=0;
	for(int i = 0; i < x; i++)
	{
		b+=n[x][i]*m[i];
		//std::cout << "\t+++" << n[x][i] << "*" << m[i] << "\n";
	}
	return b;
}

//Через конструктор загружается кофиг и создается итоговые таблицы
class table_of_calc
{
	public:
		double epsilon[26] = {0.016, 0.088, 0.184, 0.270, 0.202, 0.141, 
						0.061, 0.024, 0.01, 0.003, 0.001, 
						0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
						
		double macro_sigma_tr[100];	
		double D[100];
		double B[100];
		double macro_sigma_a[100];
		double sigma_R[100][100];
		long double concentration[100];
		double b4c;
					
		double power, volume;
		double H;
		double R;
		
		int iter_number=0;
		int step_days=0;
		
		double F0[100];
		double moderator_temperature;//температура среды
		double mod_26;
		
		double temperature[100];
		
		int rods_on=0;
		double conc_B10, conc_B11;
		
		std::vector<table_of_element> el;//вектор из таблиц table_of_element
		int number_elements=0;
		std::string buff;
		int get_config(char conf[100])
		{
			std::string element[100];

			
			std::ifstream file(conf);
			
			file >> buff >> buff;
			file >> power >> volume;
			
			file >> buff >> buff;
			file >> H >> R;
			
			b4c=0;
			file >> buff;
			file >> moderator_temperature;
			file >> buff;
			file >> iter_number;
			file >> buff;
			file >> step_days;
			file >> buff;
			file >> rods_on;
			file >> buff >> buff;
			file >> conc_B10 >> conc_B11;
			file >> buff>> buff>> buff;
			while(!file.eof())
			{
				file >> element[number_elements] >> concentration[number_elements] >> temperature[number_elements];
				//concentration[number_elements]=concentration[number_elements];
				number_elements++;
			}
			
			char buff[100];
			for(int i=0;i<number_elements;i++)
			{
				strcpy(buff, element[i].c_str());
				el.push_back(table_of_element(buff,buff));
			}
			
		}
		
		int calc()
		{
			//помещаем B4C
			concentration[0]=b4c*4*(conc_B10/100);
			concentration[1]=b4c*4*(conc_B11/100);
			concentration[2]=b4c;
			
			//очищаем прежнии значения
			for(int i=0;i<26;i++)
			{
				macro_sigma_tr[i]=0;
				macro_sigma_a[i]=0;
			}
			for(int x = 0; x < 100; x++)
			{
				for(int y=0;y<100;y++)
				{
					sigma_R[x][y]=0;
				}
			}
			//вводим поправки на самоэкранировку
			
			double sum_all_el[26];
			
			for(int k = 0; k < el.size(); k++)
			{
				for(int i = 0; i < 26; i++)
				{
					sum_all_el[i]+=el[k].sigma_tot[i]*concentration[k];
				}
			}
			
			for(int i = 0; i < 26; i++)
			{
				for(int k = 0; k < el.size(); k++)
				{
					if(concentration[k]==0)
					{
						el[k].sigma_0[i]=0;
						continue;
					}
					el[k].sigma_0[i]=(sum_all_el[i]-el[k].sigma_tot[i]*concentration[k])/concentration[k];
				}
			}
			
			for(int i = 0; i < el.size(); i++)
			{
				if(concentration[i]==0)
					continue;
				
				//std::cout << "Поправка на " << el[i].get_name() << "\n";
				
				for(int m = 0; m < 25; m++)
				{
				el[i].sigma_f[m]=el[i].const_sigma_f[m]*inter_mod(el[i].table_f, el[i].sigma_0[m], temperature[i], m+1);
				el[i].sigma_c[m]=el[i].const_sigma_c[m]*inter_mod(el[i].table_c, el[i].sigma_0[m], temperature[i], m+1);
				el[i].sigma_a[m]=el[i].sigma_f[m]+el[i].sigma_c[m];
				el[i].sigma_tr[m]=el[i].sigma_a[m] + el[i].sigma_in[m] + el[i].sigma_e[m]*(1-el[i].mu_e[m]);
				}
			}	
			
			//вводим поправка на 26 группу
			//ищем температуру нейтронного газа
			double neutron_gas_temperature = moderator_temperature*((1+1.4*((el[16].get_sigma_a(25)*concentration[16]+el[17].get_sigma_a(25)*concentration[16]))/
			(el[16].get_sigma_e(25)*concentration[16]*el[16].get_epsilon(25)+el[17].get_sigma_in(25)*concentration[16])*el[17].get_epsilon(25)));
			//ищем mod_26
			mod_26 = (sqrt(Pi)/2)*sqrt(293/neutron_gas_temperature);
			
			
			for(int i = 0; i < el.size(); i++)
			{
				el[i].sigma_f[25]=el[i].const_sigma_f[25]*mod_26;
				el[i].sigma_c[25]=el[i].const_sigma_c[25]*mod_26;
				el[i].sigma_a[25]=el[i].sigma_f[25]+el[i].sigma_c[25];
				el[i].sigma_tr[25]=el[i].sigma_a[25] + el[i].sigma_in[25] + el[i].sigma_e[25]*(1-el[i].mu_e[25]);
				//std::cout << el[i].get_name() << "\t" <<  el[i].sigma_c[25] << "\n";
			}	
			
			//заполняем macro_sigma_tr
			for(int i=0;i<26;i++)
			{
				for(int m=0;m<number_elements;m++)
				{
					macro_sigma_tr[i]+=el[m].get_sigma_tr(i)*concentration[m];
				}
			//std::cout << macro_sigma_tr[i] << "\n";
			}
			//заполняем D
			for(int i=0;i<26;i++)
			{
				D[i]=(1/(3*macro_sigma_tr[i]));
				//std::cout << D[i] << "\n";
			}
			//заполняем B
			for(int i=0;i<26;i++)
			{
			B[i]=pow(Pi/(H+(2*0.7104/macro_sigma_tr[i])),2) + pow((2.405/(R+(0.7104/macro_sigma_tr[i]))),2);//для цилиндра
			//std::cout << B[i] << "\n";
			}
			//Заполняем macro_sigma_a
			for(int i=0;i<26;i++)
			{
				for(int m=0;m<number_elements;m++)
				{
					macro_sigma_a[i]+=el[m].get_sigma_a(i)*concentration[m];
				}
			//std::cout << macro_sigma_a[i] << "\n";
			}
			//заполняем sigma_R
			for(int x = 0; x < 100; x++)
			{
				for(int y = 0; y < 100; y++)
				{
					for(int m=0;m<number_elements;m++)
					{
						sigma_R[x][y]+=el[m].get_sigma_e_plus_i_26_ik(x,y)*concentration[m];
					}
				}
			}						
		}
		
};

int zero_group(table_of_calc A, double F0[])
{
	for(int i=0; i < 26; i++)
	{
		F0[i]=(A.epsilon[i]+sum_of_multiplying(i, A.sigma_R, F0, i))/(A.D[i]*A.B[i]+A.macro_sigma_a[i]+sum_of_string(A.sigma_R,i));
		//std::cout << sum_of_string(A.sigma_R,i) << "\n";
	}
}

int next_group(table_of_calc A, double First[], double Second[])
{
			double nuFeF[100];//Ню*сигма_ф*F0
			double SnuFeF[100];//Сигма(Ню*сигма_ф*F0)
			double buff=0;
			for(int i=0; i < 26; i++)
			{
				buff=0;
				for(int x=0; x < A.el.size(); x++)
				{
					buff+=A.el[x].get_nu(i)*A.el[x].get_sigma_f(i)*A.concentration[x];
				}
				nuFeF[i]=First[i]*buff;
			}
	
			//заполняем Сигма(Ню*сигма_ф*F0) SnuFeF
			double sum;
			for(int x=0; x < 26; x++)
			{
				sum+=nuFeF[x];
			}
			
			for(int i=0; i < 26; i++)
			{
				SnuFeF[i] = sum-nuFeF[i];
				//std::cout << SnuFeF[i] << "\n";
			}
			
			//Ищем 1 итерацию
			double sigmaF_nu_eps_com[100];
			
			for(int i=0; i < 26; i++)
			{
				sigmaF_nu_eps_com[i]=0;
			}
			
			for(int i=0; i < 26; i++)
			{
				for(int x=0; x < A.el.size(); x++)
				{
					sigmaF_nu_eps_com[i]+=A.el[x].get_sigma_f(i)*A.concentration[x]*A.el[x].get_nu(i);
				}
				//std::cout << sigmaF_com[i] << "\n";
			}
			
			for(int i=0; i < 26; i++)
			{
				Second[i]=(A.epsilon[i]*SnuFeF[i]+sum_of_multiplying(i, A.sigma_R, Second, i))/(A.D[i]*A.B[i]+A.macro_sigma_a[i]+sum_of_string(A.sigma_R,i)-sigmaF_nu_eps_com[i]*A.epsilon[i]);
			}

}

double calculate_k_effective(table_of_calc A, double spector[])
{
	//Определяем доли спектра
	double perc[100];
	double sum_sp=0;
	double nu_sigma_f, DB, sigma_a;
	for(int i = 0; i < 26; i++)
	{
		sum_sp+=spector[i];
	}
	
	for(int i = 0; i < 26; i++)
	{
		perc[i]=spector[i]/sum_sp;
	}
	
	for(int i = 0; i < 26; i++)
	{
		DB+=A.D[i]*A.B[i]*perc[i];
		sigma_a+=A.macro_sigma_a[i]*perc[i];
		for(int x=0; x < A.el.size(); x++)
		{
			nu_sigma_f+=A.el[x].get_sigma_f(i)*A.concentration[x]*A.el[x].get_nu(i)*perc[i];
		}
	}
	return (nu_sigma_f/(sigma_a+DB));
}

int main(int argc, char *argv[])
{	
	std::vector<RESULT> results;
	RESULT buffer;
	
	setlocale(LC_ALL, "Russian");
	
	int day_count=0;
	buffer.Day=day_count;
	
	//Создаем конечную рассчетную таблицу
	table_of_calc A;
	A.get_config(argv[1]);
	
	
	//Считаем поправку на 26 группу
	
	
	double F0[100];
	double Last_Iter[100];
	
	do{
		
	for(int i = 0; i < 26; i++)
	{
		for(int x=0; x < A.el.size(); x++)
		{
			buffer.SIGMA_F[i]=0;
			buffer.SIGMA_C[i]=0;
			buffer.SIGMA_A[i]=0;
		}
	}	
		
	//std::cout << "Поправка на 26 группу равна " << A.mod_26 << "\n\n";
	A.calc();
	
	buffer.Day=day_count;
	//Концентрации
	std::cout << "День "<< day_count <<" \n";
	for(int i = 0; i < A.el.size(); i++)
	{
		//std::cout << A.el[i].get_name() << "\t" << A.concentration[i] << "\n";
	}

	for(int i=0;i < A.el.size(); i++)
	{
		buffer.Name_Concentration[i]= A.el[i].get_name();
		buffer.Concentration[i]=A.concentration[i];
	}

	//Считаем поток в нулевой итерации
	zero_group(A,F0);
	
	//ищем n следующих итераций		
	for(int i=0; i < A.iter_number; i++)
	{
		next_group(A, F0, Last_Iter);
		for(int i=0; i < 26; i++)
		{
			F0[i]=Last_Iter[i];
		}

	}
	
	///////////////////////////////////////////////////////////////
	/////Зачем не знаю, но по другом не работает//////////////////
	//////////////////////////////////////////////////////////////
	
	A.calc();	
	//Концентрации
	
	//ищем n следующих итераций		
	for(int i=0; i < A.iter_number; i++)
	{
		next_group(A, F0, Last_Iter);
		for(int i=0; i < 26; i++)
		{
			F0[i]=Last_Iter[i];
		}

	}

	///////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

	std::cout << "\nKeff:\t" << calculate_k_effective(A, Last_Iter) << "\n";
	buffer.K_effective=calculate_k_effective(A, Last_Iter);
	//Ищем концентрацию B4C при которой keff=1 если rods_on=1;
	
	if(A.rods_on==1)
	{		
		while(calculate_k_effective(A, Last_Iter)>1)
		{

		A.b4c+=0.0000005;
		A.calc();
		for(int i=0; i < A.iter_number; i++)
		{
			next_group(A, F0, Last_Iter);
			for(int i=0; i < 26; i++)
			{
				F0[i]=Last_Iter[i];
			}
		}

		//std::cout << 1/pow(10,b) << "\t";
	}
	A.calc();
	}
	
	for(int i=0; i < A.iter_number; i++)
	{
		next_group(A, F0, Last_Iter);
		for(int i=0; i < 26; i++)
		{
			F0[i]=Last_Iter[i];
		}

	}
	
		//std::cout << "\nKeff b4c:\t" << calculate_k_effective(A, Last_Iter) << " Концентрация " << A.b4c << "\n";
	
	//Находим доли от птока
	double perc_last_iter[100];
	double sum_sp=0;
	for(int i = 0; i < 26; i++)
	{
		sum_sp+=Last_Iter[i];
	}
	
	for(int i = 0; i < 26; i++)
	{
		perc_last_iter[i]=Last_Iter[i]/sum_sp;
	}
	//Определяем суммарный поток при ном мощности
	double sum_sigma_f=0;
	for(int i = 0; i < 26; i++)
	{
		for(int x=0; x < A.el.size(); x++)
		{
			sum_sigma_f+=A.el[x].get_sigma_f(i)*A.concentration[x]*perc_last_iter[i];
		}
	}

	for(int i = 0; i < 26; i++)
	{
		for(int x=0; x < A.el.size(); x++)
		{
			buffer.SIGMA_F[i]+=A.el[x].get_sigma_f(i)*A.concentration[x];
		}
	}

	for(int i = 0; i < 26; i++)
	{
		for(int x=0; x < A.el.size(); x++)
		{
			buffer.SIGMA_C[i]+=A.el[x].get_sigma_c(i)*A.concentration[x];
		}
	}
	
		for(int i = 0; i < 26; i++)
	{
		for(int x=0; x < A.el.size(); x++)
		{
			buffer.SIGMA_A[i]+=A.el[x].get_sigma_a(i)*A.concentration[x];
		}
	}

	double Fsum = (A.power*pow(10,6))/(200*1.6*pow(10,-13)*A.volume*sum_sigma_f);
	
	//std::cout << "Средняя плотность потока при номинальной мощности: " << Fsum << "\n";
	buffer.sum_flux = Fsum;
	//Доли от суммарного потока при ном мощности
	double calcFsum[100];
	
	for(int i = 0; i < 26; i++)
	{
		calcFsum[i] = Fsum*perc_last_iter[i];
		buffer.flux[i]=calcFsum[i];
	}
	
	
	//Ищем долю запаздывающих нейтронов
	double betta_i_t[26]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double betta_i_l[26]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double betta_i[26]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	for(int i = 0; i < 26; i++)
	{
		for(int x = 0; x < A.el.size(); x++)
		{
			betta_i_t[i] += (A.el[x].betta[i]*A.el[x].get_nu(i)*A.el[x].get_sigma_f(i)*A.concentration[x]);
			betta_i_l[i] += (A.el[x].get_nu(i)*A.el[x].get_sigma_f(i)*A.concentration[x]);
			betta_i[i]=betta_i_t[i]/betta_i_l[i];
		}
	}	
	
	double bvfi[26];
	double vfi[26];
	
	for(int i = 0; i < 26; i++)
	{
		for(int x = 0; x < A.el.size(); x++)
		{
			bvfi[i] += (betta_i[i]*A.el[x].get_nu(i)*A.el[x].get_sigma_f(i)*A.concentration[x]*perc_last_iter[i]);
			vfi[i] += (A.el[x].get_nu(i)*A.el[x].get_sigma_f(i)*A.concentration[x]*perc_last_iter[i]);
		}
	}	
	
	double bvfi_s=0;
	double vfi_s=0;
	
	for(int i = 0; i < 26; i++)
	{
		bvfi_s+=bvfi[i];
		vfi_s+=vfi[i];
	}
	
	double BETTA=bvfi_s/vfi_s;
	
	buffer.betta_neutron=BETTA;
	//std::cout << "\nДоля зап. нейтронов " << BETTA*100 << " %" << "\n";
	
	//Ищем среднее время жизни зап. нейтрона
	//Ищем бетта_j
	double Betta_j[100];
	
	for(int i = 0; i < A.el.size(); i++)
	{
		Betta_j[i]=bvfi[i]/vfi[i];
	}	
	
	//
	
	double nu_sigma_f[100];
	
	for(int i = 0; i < A.el.size(); i++)
	{
		for(int x = 0; x < 26; x++)
		{
			nu_sigma_f[i] += (A.el[i].get_nu(x)*A.el[i].get_sigma_f(x)*A.concentration[i]*perc_last_iter[x]);
		}
	}
	
	//ищем  l_del_sum;
	double l_del_sum_buff1=0;
	double l_del_sum_buff2=0;
	double l_del_sum=0;
	double DELAY_TIME=0;
	for(int i = 0; i < A.el.size(); i++)
	{
		l_del_sum_buff1+=Betta_j[i]*nu_sigma_f[i]*A.el[i].get_L_delay();
		l_del_sum_buff2+=Betta_j[i]*nu_sigma_f[i];
	}
	
	l_del_sum=l_del_sum_buff1/l_del_sum_buff2;
	DELAY_TIME=l_del_sum*BETTA;
	buffer.lifetime_neutron=DELAY_TIME;
	//std::cout << "\nСреднее время жизни нейтронов " << DELAY_TIME << "\n";

	//подготавливаем микросечения для продуктов деления
	//в группах в долях
	double s232a[26];
	double s232c[26];
	double s232f[26];
	double s233a[26];
	double s233c[26];
	double s233f[26];
	double s234a[26];
	double s234c[26];
	double s234f[26];
	double s235a[26];
	double s235c[26];
	double s235f[26];
	double s236a[26];
	double s236c[26];
	double s236f[26];
	double s238a[26];
	double s238c[26];
	double s238f[26];
	double s239a[26];
	double s239c[26];
	double s239f[26];
	double s240a[26];
	double s240c[26];
	double s240f[26];
	double s241a[26];
	double s241c[26];
	double s241f[26];
	double s242a[26];
	double s242c[26];
	double s242f[26];
	//среднее
	double sr_s232a=0;
	double sr_s232c=0;
	double sr_s232f=0;
	double sr_s233a=0;
	double sr_s233c=0;
	double sr_s233f=0;
	double sr_s234a=0;
	double sr_s234c=0;
	double sr_s234f=0;
	double sr_s235a=0;
	double sr_s235c=0;
	double sr_s235f=0;
	double sr_s236a=0;
	double sr_s236c=0;
	double sr_s236f=0;
	double sr_s238a=0;
	double sr_s238c=0;
	double sr_s238f=0;
	double sr_s239a=0;
	double sr_s239c=0;
	double sr_s239f=0;
	double sr_s240a=0;
	double sr_s240c=0;
	double sr_s240f=0;	
	double sr_s241a=0;
	double sr_s241c=0;
	double sr_s241f=0;	
	double sr_s242a=0;
	double sr_s242c=0;
	double sr_s242f=0;

	//рассчет в долях
	for(int i = 0; i < 26; i++)
	{
		s232a[i]=perc_last_iter[i]*A.el[3].get_sigma_a(i)*pow(10,-24);
		s232c[i]=perc_last_iter[i]*A.el[3].get_sigma_c(i)*pow(10,-24);
		s232f[i]=perc_last_iter[i]*A.el[3].get_sigma_f(i)*pow(10,-24);
		s233a[i]=perc_last_iter[i]*A.el[4].get_sigma_a(i)*pow(10,-24);
		s233c[i]=perc_last_iter[i]*A.el[4].get_sigma_c(i)*pow(10,-24);
		s233f[i]=perc_last_iter[i]*A.el[4].get_sigma_f(i)*pow(10,-24);
		s234a[i]=perc_last_iter[i]*A.el[5].get_sigma_a(i)*pow(10,-24);
		s234c[i]=perc_last_iter[i]*A.el[5].get_sigma_c(i)*pow(10,-24);
		s234f[i]=perc_last_iter[i]*A.el[5].get_sigma_f(i)*pow(10,-24);
		s235a[i]=perc_last_iter[i]*A.el[6].get_sigma_a(i)*pow(10,-24);
		s235c[i]=perc_last_iter[i]*A.el[6].get_sigma_c(i)*pow(10,-24);
		s235f[i]=perc_last_iter[i]*A.el[6].get_sigma_f(i)*pow(10,-24);
		s236a[i]=perc_last_iter[i]*A.el[7].get_sigma_a(i)*pow(10,-24);
		s236c[i]=perc_last_iter[i]*A.el[7].get_sigma_c(i)*pow(10,-24);
		s236f[i]=perc_last_iter[i]*A.el[7].get_sigma_f(i)*pow(10,-24);
		s238a[i]=perc_last_iter[i]*A.el[8].get_sigma_a(i)*pow(10,-24);
		s238c[i]=perc_last_iter[i]*A.el[8].get_sigma_c(i)*pow(10,-24);
		s238f[i]=perc_last_iter[i]*A.el[8].get_sigma_f(i)*pow(10,-24);
		s239a[i]=perc_last_iter[i]*A.el[9].get_sigma_a(i)*pow(10,-24);
		s239c[i]=perc_last_iter[i]*A.el[9].get_sigma_c(i)*pow(10,-24);
		s239f[i]=perc_last_iter[i]*A.el[9].get_sigma_f(i)*pow(10,-24);
		s240a[i]=perc_last_iter[i]*A.el[10].get_sigma_a(i)*pow(10,-24);
		s240c[i]=perc_last_iter[i]*A.el[10].get_sigma_c(i)*pow(10,-24);
		s240f[i]=perc_last_iter[i]*A.el[10].get_sigma_f(i)*pow(10,-24);	
		s241a[i]=perc_last_iter[i]*A.el[11].get_sigma_a(i)*pow(10,-24);
		s241c[i]=perc_last_iter[i]*A.el[11].get_sigma_c(i)*pow(10,-24);	
		s241f[i]=perc_last_iter[i]*A.el[11].get_sigma_f(i)*pow(10,-24);
		s242a[i]=perc_last_iter[i]*A.el[12].get_sigma_a(i)*pow(10,-24);
		s242c[i]=perc_last_iter[i]*A.el[12].get_sigma_c(i)*pow(10,-24);
		s242f[i]=perc_last_iter[i]*A.el[12].get_sigma_f(i)*pow(10,-24);		
		
	}
	//ищем среднее
	for(int i = 0; i < 26; i++)
	{
		sr_s232a+=s232a[i];
		sr_s232c+=s232c[i];
		sr_s232f+=s232f[i];
		sr_s233a+=s233a[i];
		sr_s233c+=s233c[i];
		sr_s233f+=s233f[i];
		sr_s234a+=s234a[i];
		sr_s234c+=s234c[i];
		sr_s234f+=s234f[i];
		sr_s235a+=s235a[i];
		sr_s235c+=s235c[i];
		sr_s235f+=s235f[i];
		sr_s236a+=s236a[i];
		sr_s236c+=s236c[i];
		sr_s236f+=s236f[i];
		sr_s238a+=s238a[i];
		sr_s238c+=s238c[i];
		sr_s238f+=s238f[i];
		sr_s239a+=s239a[i];
		sr_s239c+=s239c[i];
		sr_s239f+=s239f[i];
		sr_s240a+=s240a[i];
		sr_s240c+=s240c[i];	
		sr_s240f+=s240f[i];	
		sr_s241a+=s241a[i];
		sr_s241c+=s241c[i];
		sr_s241f+=s241f[i];		
		sr_s242a+=s242a[i];
		sr_s242c+=s242c[i];
		sr_s242f+=s242f[i];			
	}

	//Ищем КВ
	double kw = (sr_s232c*A.concentration[3]+sr_s238c*A.concentration[8])/
	(sr_s235a*A.concentration[6]+sr_s239a*A.concentration[9]+sr_s241a*A.concentration[11]+sr_s233a*A.concentration[4]);
	buffer.breeding_ratio=kw;
	//std::cout << "\nКоэф. воспроизводства " << kw << "\n";
	
	//пересчитываем концентрации

	//OD9
	A.concentration[15]=A.concentration[15]+A.concentration[9]*sr_s239f*Fsum*A.step_days*24*60*60+A.concentration[10]*sr_s240f*Fsum*A.step_days*24*60*60
	+A.concentration[11]*sr_s241f*Fsum*A.step_days*24*60*60+A.concentration[12]*sr_s242f*Fsum*A.step_days*24*60*60;	
	//OD5
	A.concentration[14]=A.concentration[14]+A.concentration[6]*sr_s235f*Fsum*A.step_days*24*60*60+A.concentration[7]*sr_s236f*Fsum*A.step_days*24*60*60
	+A.concentration[8]*sr_s238f*Fsum*A.step_days*24*60*60;	
	//OD3
	A.concentration[13]=A.concentration[13]+A.concentration[3]*sr_s232f*Fsum*A.step_days*24*60*60+A.concentration[4]*sr_s233f*Fsum*A.step_days*24*60*60
	+A.concentration[5]*sr_s234f*Fsum*A.step_days*24*60*60;	
	//Pu242
	A.concentration[12]=A.concentration[12]-A.concentration[12]*sr_s242a*Fsum*A.step_days*24*60*60+(A.concentration[11]*sr_s241c*Fsum*A.step_days*24*60*60);
	//Pu241
	A.concentration[11]=A.concentration[11]-A.concentration[11]*sr_s241a*Fsum*A.step_days*24*60*60+(A.concentration[10]*sr_s240c*Fsum*A.step_days*24*60*60);
	//Pu240
	A.concentration[10]=A.concentration[10]-A.concentration[10]*sr_s240a*Fsum*A.step_days*24*60*60+(A.concentration[9]*sr_s239c*Fsum*A.step_days*24*60*60);
	//Pu239
	A.concentration[9]=A.concentration[9]-A.concentration[9]*sr_s239a*Fsum*A.step_days*24*60*60+(A.concentration[8]*sr_s238c*Fsum*A.step_days*24*60*60);
	//U238
	A.concentration[8]=A.concentration[8]-A.concentration[8]*sr_s238a*Fsum*A.step_days*24*60*60;
	//U236
	A.concentration[7]=A.concentration[7]-A.concentration[7]*sr_s236a*Fsum*A.step_days*24*60*60+(A.concentration[6]*sr_s235c*Fsum*A.step_days*24*60*60);
	//U235
	A.concentration[6]=A.concentration[6]-A.concentration[6]*sr_s235a*Fsum*A.step_days*24*60*60+(A.concentration[5]*sr_s234c*Fsum*A.step_days*24*60*60);
	//U234
	A.concentration[5]=A.concentration[5]-A.concentration[5]*sr_s234a*Fsum*A.step_days*24*60*60+(A.concentration[4]*sr_s233c*Fsum*A.step_days*24*60*60);
	//U233
	A.concentration[4]=A.concentration[4]-(A.concentration[4]*sr_s233a*Fsum*A.step_days*24*60*60)+(A.concentration[3]*sr_s232c*Fsum*A.step_days*24*60*60);
	//Th232
	A.concentration[3]=A.concentration[3]-A.concentration[3]*sr_s232a*Fsum*A.step_days*24*60*60;
	
	day_count+=A.step_days;
	A.b4c=0;
	A.calc();
		for(int i=0; i < A.iter_number; i++)
	{
		next_group(A, F0, Last_Iter);
		for(int i=0; i < 26; i++)
		{
			F0[i]=Last_Iter[i];
		}

	}
	
	results.push_back(buffer);
	
}while(calculate_k_effective(A, Last_Iter)>1);

//Формируем отчет в файл

    std::ofstream out;
    out.open(argv[2]);
	
	out << std::setw(20) << "EFFECTIVE DAYS" << std::setw(25);
	
	for(int u=0; u < results.size(); u++)
	{
		out << results[u].Day << std::setw(25);
	}
	
	out << "\n\n";
	
	for(int i=0;i<A.el.size();i++)
	{ 
	out << std::setw(20) <<results[0].Name_Concentration[i] << std::setw(25);
		for(int u=0; u < results.size(); u++)
		{
			out << results[u].Concentration[i] << std::setw(25) << std::setprecision(10);
		}
		out << "\n";
	} 
    
    out << "\n";
    
	out << std::setw(20) << "K-eff" << std::setw(25);
	for(int u=0; u < results.size(); u++)
	{
		out << results[u].K_effective << std::setw(25) << std::setprecision(10);
	}
	out << "\n";
	
		out << std::setw(20) << "AVARAGE FLUX" << std::setw(25);
	for(int u=0; u < results.size(); u++)
	{
		out << results[u].sum_flux << std::setw(25) << std::setprecision(10);
	}
	out << "\n";

	out << std::setw(20) << "DEL.NEUTRON.FRACT" << std::setw(25);
	for(int u=0; u < results.size(); u++)
	{
		out << results[u].betta_neutron << std::setw(25) << std::setprecision(10);
	}
	out << "\n";
	
	out << std::setw(20) << "BREEDING_RATIO" << std::setw(25);
	for(int u=0; u < results.size(); u++)
	{
		out << results[u].breeding_ratio << std::setw(25) << std::setprecision(10);
	}
	out << "\n";
	
	out << std::setw(20) << "DEL.NEUTRON.LIFETIME" << std::setw(25);
	for(int u=0; u < results.size(); u++)
	{
		out << results[u].lifetime_neutron << std::setw(25) << std::setprecision(10);
	}
	out << "\n\n";
	
	out << std::setw(20) <<"*********************************FLUX*********************************\n";

	for(int g=0; g < GROUPS; g++)
	{
		out << std::setw(20) << g+1;
		for(int u=0; u < results.size(); u++)
		{
			out << std::setw(25) << results[u].flux[g] << std::setw(25);
		}
		out << "\n";
	}	
	out << "\n";
    
    	out << std::setw(20) <<"*********************************Sigma-F*********************************\n";

	for(int g=0; g < GROUPS; g++)
	{
		out << std::setw(20) << g+1;
		for(int u=0; u < results.size(); u++)
		{
			out << std::setw(25) << results[u].SIGMA_F[g] << std::setw(25);
		}
		out << "\n";
	}	
	out << "\n";
	
	    	out << std::setw(20) <<"*********************************Sigma-C*********************************\n";

	for(int g=0; g < GROUPS; g++)
	{
		out << std::setw(20) << g+1;
		for(int u=0; u < results.size(); u++)
		{
			out << std::setw(25) << results[u].SIGMA_C[g] << std::setw(25);
		}
		out << "\n";
	}	
	out << "\n";
    
        	out << std::setw(20) <<"*********************************Sigma-A*********************************\n";

	for(int g=0; g < GROUPS; g++)
	{
		out << std::setw(20) << g+1;
		for(int u=0; u < results.size(); u++)
		{
			out << std::setw(25) << results[u].SIGMA_A[g] << std::setw(25);
		}
		out << "\n";
	}	
	out << "\n";
    
}













