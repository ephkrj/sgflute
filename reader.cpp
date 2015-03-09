#include<iostream>
#include<fstream>
#include"reader.h"
#include<vector>
#include<sstream>
#include<iterator>
#include<iomanip>
#include<algorithm>

int main(){

	std::string name;
	std::cout<<"Please input file name to be converted \n(make sure it's in the correct directory!): ";
	std::cin>>name;
	//zone-subzone relation list
        int zone_subzone[]={1,0,1,1,1,2,1,3,1,4,1,5,1,6,1,7,1,8,2,9,2,10,2,11,2,12,2,13,2,14,2,15,2,16,3,17,3,18,3,19,4,215,4,216,4,217,4,218,
5,20,5,21,5,22,5,23,5,24,5,25,5,26,5,27,5,219,6,28,6,29,6,30,6,31,6,32,6,33,6,34,6,35,6,36,6,37,6,38,6,39,6,40,6,41,
6,42,6,42,6,42,7,43,7,44,7,45,7,46,7,47,7,48,7,49,8,50,8,51,8,52,8,53,8,54,8,55,8,56,8,57,9,220,10,58,10,59,10,59,
11,221,12,60,12,61,12,62,12,63,12,64,12,65,12,66,13,67,13,68,13,69,13,70,13,71,13,72,13,73,13,74,13,74,14,75,14,76,
14,76,14,76,14,76,14,76,14,76,14,76,14,76,14,76,14,76,14,76,15,77,15,78,15,79,15,80,15,81,16,82,16,83,16,84,16,85,16,86,
16,87,16,87,17,88,17,89,17,90,17,91,17,92,17,93,17,94,17,94,17,94,17,94,18,95,18,96,18,97,18,98,18,99,18,100,18,101,
18,101,18,101,19,102,19,103,19,104,19,105,19,106,19,107,19,108,19,109,19,109,20,222,21,110,22,223,23,224,24,111,24,112,
24,113,24,225,24,226,25,227,25,228,25,229,26,114,26,115,26,116,26,117,26,117,26,117,27,230,28,118,28,119,28,120,28,121,
28,122,29,231,29,232,29,233,30,214,31,123,31,124,31,125,31,126,32,127,32,128,32,129,32,130,32,131,32,132,32,132,33,234,
33,235,33,236,33,237,33,238,34,239,34,240,34,241,34,242,34,243,35,133,35,134,35,135,35,136,35,137,35,137,36,138,36,139,
36,140,36,141,36,142,36,143,36,144,36,145,36,146,36,147,36,148,36,148,36,148,36,148,36,244,37,149,37,150,37,151,37,152,
37,153,38,154,38,155,38,156,38,157,38,158,38,159,38,159,38,159,38,159,38,159,39,245,40,160,40,161,40,162,40,163,40,164,
40,165,40,165,40,165,40,165,41,166,41,167,41,168,41,169,41,170,41,171,42,172,42,173,42,174,42,175,43,246,44,176,44,177,
44,177,45,247,45,248,46,249,47,250,47,251,47,252,47,253,47,254,48,178,48,179,48,180,48,181,48,255,49,182,49,183,49,184,
49,185,50,256,51,186,51,187,51,188,51,189,51,190,51,191,51,192,51,193,51,194,51,195,51,196,51,197,52,257,52,258,52,259,
52,260,52,261,52,262,53,263,53,264,54,265,55,198,55,199,55,200,55,201,55,202,55,203,55,204,55,204,55,204,56,205,56,206,
56,207,56,208,56,209,56,210,56,211,56,212,56,213};
	
	std::vector< int > zone_values;
	std::vector<Person> pvec;
	std::vector<Zone> zvec;
	std::ostringstream oss;
	std::string dummy;
	oss.str(name+".csv");//input filename
	std::ifstream iss(oss.str().c_str());
	
	if(!iss){
		std::cout<<"Sorry, "<<name<<" not found. Could it be somewhere else?\n";
		exit(-1);
	}
	getline(iss,dummy);
	while(iss.good()){
		Person p;
		char c;	iss>>p.id>>c>>p.exactAge>>c>>p.gender>>c>>p.race>>c>>p.citizen>>c>>p.marital>>c>>p.number_children>>c>>p.nSchoolType>>c>>p.qualification>>c>>p.work>>c>>
p.income>>c>>p.headhouse>>c>>p.occupation>>c>>p.religion>>c>>p.partnerid>>c>>p.family>>c>>p.dwelling>>c>>p.hasmaid>>c>>p.mobility>>c>>p.nhours>>c>>
p.industry>>c>>p.NS>>c>>p.trans_mode>>c>>p.trans_time>>c>>p.oversea>>c>>p.nFamilySize>>c>>p.f_nuclei>>c>>p.generation>>c>>p.youngest>>c>>p.numwork>>c>>
p.hincome>>c>>p.momid>>c>>p.dadid>>c>>p.unit_id>>c>>p.nHomeComm>>c>>p.postcode_home>>c>>p.x_coor_home>>c>>p.y_coor_home>>c>>p.nWorkplace>>c>>p. subzone_work>>c>>p.postcode_work>>c>>p.x_coor_work>>c>>p.y_coor_work>>c>>p.nSchplace>>c>>p.nSchComm>>c>>p.nSchNeighborhood>>c>>
p.x_coor_school>>c>>p.y_coor_school;
	if(iss.eof()) break;
		pvec.push_back(p);
		
	}
	iss.close();
	//pvec.pop_back();
//convert from subzone to zone
	for(std::vector<Person>::iterator it = pvec.begin();it!=pvec.end();it++){
		Person &p=*it;
		for(int i=1;i<608;i=i+2){
			if(p.nHomeComm==zone_subzone[i]){
				p.zone_home=zone_subzone[i-1];
				break;
			}
		}
	}

	for(std::vector<Person>::iterator it = pvec.begin();it!=pvec.end();it++)
	{
		Person &p = *it;
		for(int i=1;i<608;i=i+2){
				if(p.subzone_work==zone_subzone[i]&&p.subzone_work!=-1){
					p.zone_work=zone_subzone[i-1];
					break;
			}
		}
	}
	for(std::vector<Person>::iterator it = pvec.begin();it!=pvec.end();it++){
		Person &p=*it;
		if(p.exactAge>=0 && p.exactAge<=4)
			p.age=0;
		else if(p.exactAge>=5 && p.exactAge<=18)
			p.age=1;
		else if(p.exactAge>=19 && p.exactAge<=29)
			p.age=2;
		else if(p.exactAge>=30 && p.exactAge<=64)
			p.age=3;
		else 
			p.age=4;

	}
	std::vector< int > zone_work_id;
	for(std::vector<Person>::iterator it = pvec.begin();it!=pvec.end();it++){
		Person &p = *it;
		zone_work_id.push_back(p.postcode_work);
		zone_work_id.push_back(p.nSchNeighborhood);
	}
	sort(zone_work_id.begin(),zone_work_id.end());
	zone_work_id.erase(unique(zone_work_id.begin(),zone_work_id.end()),zone_work_id.end());
	for(int i=0;i<zone_work_id.size();i++)		
		for(std::vector<Person>::iterator it = pvec.begin();it!=pvec.end();it++){
			Person &p=*it;
			if(p.postcode_work==zone_work_id[i] && p.nSchNeighborhood==-1){
				p.nWorkplace=i;
			}
			if(p.postcode_work==-1 && p.nSchNeighborhood==zone_work_id[i]){
				p.nWorkplace=i;
			}
			if(p.postcode_work==-1 && p.nSchNeighborhood==-1)
				p.nWorkplace=-1;
		}
	std::string Revise=name+"-Revised input.csv";
	std::ofstream file(Revise.c_str());
file<<"id"<<","<<"exactAge"<<","<<"age"<<","<<"gender"<<","<<"race"<<","<<"citizen"<<","<<"marital"<<","<<"number_children"<<","<<"nSchoolType"<<","<<"qualification"<<","<<"work"<<","<<"income"<<","<<"headhouse"<<","<<"occupation"<<","<<"religion"<<","<<"partnerid"<<","<<"family"<<","<<"dwelling"<<","<<"hasmaid"<<","<<"mobility"<<","<<"nhours"<<","<<"industry"<<","<<"NS"<<","<<"trans_mode"<<","<<"trans_time"<<","<<"oversea"<<","<<"nFamilySize"<<","<<"f_nuclei"<<","<<"generation"<<","<<"youngest"<<","<<"numwork"<<","<<"hincome"<<","<<"momid"<<","<<"dadid"<<","<<"unit_id"<<","<<"zone_home"<<","<<"nHomeComm"<<","<<"postcode_home"<<","<<"x_coor_home"<<","<<"y_coor_home"<<","<<"nWorkplace"<<","<<"zone_work"<<","<<"subzone_work"<<","<<"postcode_work"<<","<<"x_coor_work"<<","<<"y_coor_work"<<","<<"nSchplace"<<","<<"nSchComm"<<","<<"nSchNeighborhood"<<","<<"x_coor_school"<<","<<"y_coor_school"<<std::endl;
	for(std::vector<Person>::iterator it = pvec.begin();it!=pvec.end();it++){
		Person &p=*it;
file<<p.id<<","<<p.exactAge<<","<<(int)p.age<<","<<p.gender<<","<<p.race<<","<<p.citizen<<","<<p.marital<<","<<p.number_children<<","<<p.nSchoolType<<","<<p.qualification<<","<<p.work<<","<<p.
income<<","<<p.headhouse<<","<<p.occupation<<","<<p.religion<<","<<p.partnerid<<","<<p.family<<","<<p.dwelling<<","<<p.hasmaid<<","<<p.mobility<<","<<p.nhours<<","<<p.
industry<<","<<p.NS<<","<<p.trans_mode<<","<<p.trans_time<<","<<p.oversea<<","<<(int)p.nFamilySize<<","<<p.f_nuclei<<","<<p.generation<<","<<p.youngest<<","<<p.numwork<<","<<p.
hincome<<","<<p.momid<<","<<p.dadid<<","<<p.unit_id<<","<<p.zone_home<<","<<p.nHomeComm<<","<<p.postcode_home<<","<<p.x_coor_home<<","<<p.y_coor_home<<","<<p.nWorkplace<<","<<p.zone_work<<
","<<p.subzone_work<<","<<p.postcode_work<<","<<p.x_coor_work<<","<<p.y_coor_work<<","<<p.nSchplace<<","<<p.nSchComm<<","<<p.nSchNeighborhood<<","<<p.x_coor_school<<","<<p.y_coor_school<<std::endl;
	}
	file.close();
//for file
/*	for(std::vector<Person>::iterator it = pvec.begin();it!=pvec.end();it++)
	{
		Person &p = *it;
		zone_values.push_back(p.zone_home);
	}
	sort(zone_values.begin(),zone_values.end());
	zone_values.erase(unique(zone_values.begin(),zone_values.end()),zone_values.end());

	for(int i= 0;i<zone_values.size();i++)
	{
		Zone z;
		z.id=zone_values[i];
		z.pop=0;
		z.employable=0;
		z.employed=0;
		z.x=0;
		z.y=0;
		std::vector<int> temp_z;
		for(std::vector<Person>::iterator it=pvec.begin();it!=pvec.end();it++)
		{
			Person &p = *it;
			if(p.zone_home==z.id)
			{
				z.pop++;
				z.x+=p.x_coor_home;
				z.y+=p.y_coor_home;
				if(p.exactAge>=16 && p.exactAge<=65)
				{
					z.employable++;
					if(p.occupation!=0)
					{
						z.employed++;
						z.workzones.push_back(p.zone_work);
					}
				}
				
			}
		}
		temp_z=z.workzones;
		sort(z.workzones.begin(),z.workzones.end());
		z.workzones.erase(unique(z.workzones.begin(),z.workzones.end()),z.workzones.end());
		for(int i=0;i<z.workzones.size();i++)
		{
			z.workpop.push_back(std::count(temp_z.begin(),temp_z.end(),z.workzones[i]));
		}
		z.x=z.x/z.pop;
		z.y=z.y/z.pop;
		zvec.push_back(z);
	}

/*
//writing into input files
	std::string tract_file=name+"-tracts.dat";
	std::string employment_file=name+"-employment.dat";
	std::string wf_file=name+"-wf.dat";
	std::ofstream file(tract_file.c_str());
	for(std::vector<Zone>::iterator it=zvec.begin();it!=zvec.end();it++)
	{
		Zone &z=*it;
		file<<1<<","<<std::setfill('0')<<std::setw(3)<<1<<","<<std::setfill('0')<<std::setw(6)<<z.id<<","<<z.pop<<","<<z.x<<","<<z.y<<endl;
	}
	file.close();
	file.open(employment_file.c_str());
	for(std::vector<Zone>::iterator it=zvec.begin();it!=zvec.end();it++)
	{
		Zone &z=*it;
		file<<1<<" "<<std::setfill('0')<<std::setw(3)<<1<<" "<<std::setfill('0')<<std::setw(6)<<z.id<<" "<<z.employed<<" "<<z.employable<<endl;
	}
	file.close();
	file.open(wf_file.c_str());
	for(std::vector<Zone>::iterator it=zvec.begin();it!=zvec.end();it++)
	{
		Zone &z=*it;
		for(int i=0;i<z.workzones.size();i++)
			file<<1<<" "<<std::setfill('0')<<std::setw(3)<<1<<" "<<std::setfill('0')<<std::setw(6)<<z.id<<" "<<1<<" "<<std::setfill('0')<<std::setw(3)<<1<<" "<<std::setfill('0')<<std::setw(6)<<z.workzones[i]<<" "<<z.workpop[i]<<endl;
	}
	file.close();*/
}
