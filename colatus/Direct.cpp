//
//  Direct.cpp
//  NBodySim
//
//  Created by José Ricardo da Silva  Júnior on 06/02/12.
//  Modified by Diego H. Stalder Diaz   18/09/12
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

/*#define By_s 3.15569252e16 //Bys/s*m/Mpc
#define m_Mpc  (3.0856776e22) //Bys/s*m/Mpc
#define _10km_m 10000 //Bys/s*m/Mpc
#define msol (1.9885e30*1e10) //Bys/s*m/Mpc
#define Grav 6.67384e-11 //Bys/s*m/Mpc
*/
#include "Direct.h"
int fromChar(char* value){
	int i;
	std::istringstream is((char*)value);
	is >> i;
	return i;
}
double fromCharF(char* value){
	double i;
	std::istringstream is((char*)value);
	is >> i;

	return i;
}

std::string fromInt(int value){
	std::stringstream out;
	out << value;
	return out.str();
}
std::string fromDouble(double value){
	std::stringstream out;
	out << value;
	return out.str();
}
double lenght(double3 vec){
	return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

namespace NBody {
namespace Simulators {

Direct::Direct(char *filename) :Pos(NULL), Vel(NULL), acc(NULL) {
	// Carregar arquivo
	sim_cfgFilename=std::string((char*)filename);
	if (!Load_sim_parameters() == 0){
		fprintf(stderr, "Error opening parameters file");
	}

	if (!LoadParticleData() == 0){
		fprintf(stderr, "Error opening IC file");
	}
}

Direct::~Direct(){
	if (Pos != NULL) free(Pos);
	if (Vel != NULL) free(Vel);
	if (acc != NULL) free(acc);
}


void Direct::Go(){

if (use_gpu){

	printf("\n  GPU ");
	DoGPUSimulation();
}
else {

	printf("\n Only CPU ");
	DoCPUSimulation();
}
write_snapshots_info();
liberar_memoriaCPU();

}

void Direct::write_snapshots_info()
{
	//char buf[60];
	char buf1[60];

	sprintf(buf, "%s.xml", mFileOut.c_str());
	xmlTextWriterPtr writer = xmlNewTextWriterFilename(buf, 0);
	            xmlTextWriterStartDocument(writer, NULL, NULL, NULL);
	            {
	            	xmlTextWriterStartElement(writer, BAD_CAST "NBody_Sim_out");
	                {
					int i;
					printf("\n Saving xml...");
					sprintf(buf, "%d", mNumBodies);
					xmlTextWriterWriteAttribute(writer, BAD_CAST "NumBodies", xmlCharStrdup(buf));
					 sprintf(buf, "%f", L);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Box", xmlCharStrdup(buf));
					 sprintf(buf, "%ssnap",mFileOut.c_str());
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Output_base_filename", xmlCharStrdup(buf));
					 sprintf(buf, "%d", snap);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Snapshots", xmlCharStrdup(buf));
					if(mBinOutput)
						sprintf(buf, "Yes");
					else
						sprintf(buf, "No");
					xmlTextWriterWriteAttribute(writer, BAD_CAST "File_binary", xmlCharStrdup(buf));

					 sprintf(buf, "%e",tcom);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Com_time", xmlCharStrdup(buf));
					 sprintf(buf, "%e",tproc);
                     xmlTextWriterWriteAttribute(writer, BAD_CAST "Proc_time", xmlCharStrdup(buf));
                     sprintf(buf, "%d",it_old);
                     xmlTextWriterWriteAttribute(writer, BAD_CAST "Time_steps", xmlCharStrdup(buf));
                     for (i=0;i<snap;i++)
					     {
						 sprintf(buf1,"Snap_%d",i);
						 sprintf(buf, "%2.4f", z_s[i]);
						 xmlTextWriterWriteAttribute(writer,BAD_CAST  xmlCharStrdup(buf1), xmlCharStrdup(buf));
					     }
					 sprintf(buf, "%s.energy.dat",mFileOut.c_str());
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Out_file_energy", xmlCharStrdup(buf));
					 sprintf(buf, "%s.i.dat",mFileOut.c_str());
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Out_file_tracking", xmlCharStrdup(buf));
					 sprintf(buf, "%s.a.dat",mFileOut.c_str());
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Out_file_scaler_factor", xmlCharStrdup(buf));
					 sprintf(buf, "%s",mICFilename.c_str());
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "IC_info", xmlCharStrdup(buf));
					 sprintf(buf, "%s", sim_cfgFilename.c_str());
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Simulation_config", xmlCharStrdup(buf));
					 printf("...done\n ");

					 xmlTextWriterEndElement(writer);

	                }
				xmlTextWriterEndElement(writer);
	            }
			xmlTextWriterEndDocument(writer);
		xmlFreeTextWriter(writer);

}


int Direct::Load_sim_parameters(){
	// Ler p do xml
	int flag;
	xmlDocPtr doc;
	xmlNodePtr cur;
	//printf("\n %s",sim_cfgFilename.c_str());

	doc = xmlParseFile(sim_cfgFilename.c_str());
		if (doc == NULL){
		fprintf(stderr, "Error reading file!");
		return -1;
	}

	cur = xmlDocGetRootElement(doc);
	if (cur == NULL){
		fprintf(stderr, "File is empty!");
		xmlFreeDoc(doc);
		return -1;
	}

	// Read data
	if (!xmlStrcmp(cur->name, (const xmlChar*) "NBody_Sim_Config")){

		xmlChar* r_buf;		std::string parameter;

		r_buf=xmlGetProp(cur, (const xmlChar*) "Device_to_run");
		parameter=std::string(( char*)r_buf);
		use_gpu=((parameter.compare("GPU"))==0);
	    parameter.empty();
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Benchmark");
	    parameter=std::string((char*)r_buf);
	    Benchmark=((parameter.compare("Yes"))== 0);
	    parameter.clear();
	    xmlFree(r_buf);


        r_buf=xmlGetProp(cur, (const xmlChar*) "Compute_energy");
	    parameter=std::string((char*)r_buf);
		sav_ener=((parameter.compare("Yes"))== 0);
	    parameter.clear();
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Save_snapshots");
	    parameter=std::string((char*)r_buf);
	    sav_snap=((parameter.compare("Yes"))==0);
	    parameter.clear();
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Max_Snapshots");
	    max_snapshots= fromChar((char*) r_buf);
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Softering_epsilon");
	    //printf("\n epsilon %s",r_buf);
	    epsilon= fromCharF((char*) r_buf);
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Max_time_steps");
	    max_steps = fromChar((char*)r_buf);
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Time_Step_Factor_Eta");
	    eta= fromCharF((char*) r_buf);
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Time_Linearization_Alpha");
	    alpha= fromCharF((char*) r_buf);
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Redshift_final");
	    zf= fromCharF((char*) r_buf);
	    xmlFree(r_buf);


	    r_buf=xmlGetProp(cur, (const xmlChar*) "Tracking_Particules");
		parameter=std::string(( char*)r_buf);
		track_on=((parameter.compare("Yes"))==0);
		parameter.empty();

	    r_buf=xmlGetProp(cur, (const xmlChar*) "HowManyParticulesTrack");
	   	int ni_snap= fromChar((char*) r_buf);
	   	xmlFree(r_buf);


	    r_buf=xmlGetProp(cur, (const xmlChar*) "IC_file");
        mICFilename=std::string((char*)r_buf);
        xmlFree(r_buf);

        r_buf=xmlGetProp(cur, (const xmlChar*) "Binary_output");
	    parameter=std::string((char*)r_buf);
	    mBinOutput=((parameter.compare("Yes"))==0);
	    parameter.clear();
	    xmlFree(r_buf);

	    r_buf=xmlGetProp(cur, (const xmlChar*) "Out_file");
	    mFileOut=std::string((char*)r_buf);

		xmlFree(r_buf);

	}
	xmlFree(doc);
	return(0);
}
int Direct::LoadParticleData(){

	// Ler dados do xml
	int flag;
	xmlDocPtr doc;
	xmlNodePtr cur;
	printf("\n Lendo parametros..");
	doc = xmlParseFile(mICFilename.c_str());


	if (doc == NULL){
		fprintf(stderr, "Error reading file!");
		return -1;
	}

	cur = xmlDocGetRootElement(doc);
	if (cur == NULL){
		fprintf(stderr, "File is empty!");
		xmlFreeDoc(doc);
		return -1;
	}

	// Read data
	if (!xmlStrcmp(cur->name, (const xmlChar*) "N-GenIC")){
		xmlChar* r_buf;

	    r_buf= xmlGetProp(cur, (const xmlChar*) "NumBodies");
	    mNumBodies= fromChar((char*) r_buf);
		xmlFree(r_buf);

	    r_buf= xmlGetProp(cur, (const xmlChar*) "Redshift");
	   	zin= fromCharF((char*) r_buf);
		xmlFree(r_buf);

	   	r_buf= xmlGetProp(cur, (const xmlChar*) "Box");
	   	L= fromCharF((char*) r_buf);
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "OmegaMatter");
		omegaM= fromCharF((char*) r_buf);
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "OmegaLambda");
		omegaL= fromCharF((char*) r_buf);
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "OmegaBaryon");
		omegaBM= fromCharF((char*) r_buf);
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "HubbleParam_10km_s_Mpc_h");
		H0= fromCharF((char*) r_buf);//10km/Mpc_h/s
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "UnitLength_in_m_Mpc_h");
    	m_Mpc= fromCharF((char*) r_buf);
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "UnitMass_in_kg_1e10_Msol_h");
		msol= fromCharF((char*) r_buf);
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "UnitVelocity_in_m_s_10Km_s");
		_10km_m= fromCharF((char*) r_buf);
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "Gravitational_Constant_m3_kg_s2");
		G= fromCharF((char*) r_buf)*msol/(m_Mpc)/_10km_m/_10km_m;//10km*10km*(Mpc_h)/1e10Ms_h/s/s;
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "Binary_file");
		std::string	parameter=std::string((char*)r_buf);
		mBinInput=((parameter.compare("Yes"))==0);
		parameter.clear();
		xmlFree(r_buf);

		r_buf= xmlGetProp(cur, (const xmlChar*) "Data_filename");
		mICFilename.clear();
		mICFilename=std::string((char*)r_buf);

		xmlFree(r_buf);
	}

	xmlFree(doc);
	printf("..particulas");

   int n;
	barionic_matter=false;


    Pos = (double3*) malloc(sizeof(double3) * mNumBodies);
	Vel = (double3*) malloc(sizeof(double3) * mNumBodies);
	acc = (double3*) malloc(sizeof(double3) * mNumBodies);

	printf("\n Reading_IC %d %s...",  mNumBodies,mICFilename.c_str());



	FILE* _pdata = fopen(mICFilename.c_str(), "rb");
	if (_pdata == NULL) {
		fprintf(stderr, "Error opening binary data!");
		free(Pos);
		free(Vel);
		free(acc);
		Pos= NULL;
		Vel=NULL;
		acc=NULL;
		return -1;
	}
	rewind(_pdata);

// Ler dados
	double3 pm, pos,vel;
	double d;
	int k=0,j,i;


	i_snap = (int*) calloc(ni_snap,sizeof(int));
	pm.x=L/2;
	pm.y=L/2;
	pm.z=L/2;
	j=0;

   for (int ii = 0; ii < (mNumBodies); ii++) {
		i=ii;
		i=i-j;
		flag=fread(&m, sizeof(double), 1, _pdata);
		flag=fread(&pos, sizeof(double), 3, _pdata);
		flag=fread(&vel, sizeof(double), 3, _pdata);

			d=(pm.x-pos.x)*(pm.x-pos.x)+(pm.y-pos.y)*(pm.y-pos.y)+(pm.z-pos.z)*(pm.z-pos.z);
			if(d<9 and k < ni_snap)
			{
				i_snap[k]=i;
				k++;
			}
		Pos[i].x=pos.x; 		Pos[i].y=pos.y;		Pos[i].z=pos.z;
     	Vel[i].x=vel.x;  		Vel[i].y=vel.y;   	Vel[i].z=vel.z;
   	}

	ni_snap=k;

	printf(" %d...done ",ni_snap);

	fclose(_pdata);
return 0;
}

void Direct::setSaveSnapshotRedshift(){
int i;
	for(int i=0;i<(max_snapshots+1);i++)
	{
			z_s[i]=pow(p0+dp_min*max_steps/(max_snapshots-1)*i ,1/alpha);
			if(z_s[i]<(1/(10+1)*1.05) and z_s[i]>(1/(10+1)*0.95))
					z_s[i]=1/11;
				if(z_s[i]<(1/(5+1)*1.05) and z_s[i]>(1/(5+1)*0.95))
							z_s[i]=1/(6);
				if(z_s[i]<(1/(3+1)*1.05) and z_s[i]>(1/(3+1)*0.95))
						z_s[i]=1/(4);
				if(z_s[i]<(1/(1+1)*1.05) and z_s[i]>(1/(1+1)*0.95))
						z_s[i]=1/(2);
				if(z_s[i]<(1/(0.5+1)*1.05) and z_s[i]>(1/(0.5+1)*0.95))
						z_s[i]=1/(0.5+1);
				if(z_s[i]<(1/(0.3+1)*1.05) and z_s[i]>(1/(0.3+1)*0.95))
							z_s[i]=1/(0.3+1);
				if(z_s[i]<(1/(1e-1+1)*1.05) and z_s[i]>(1/(1e-1+1)*0.95))
							z_s[i]=1/(1e-1+1);
		}

}
void Direct::DoCPUSimulation(){
	allocar_memoriaCPU();
	struct timeval t0,t1,tb0,tb1,t00,t01;
	//variable para el calculo del tiempo
	snap = 0;	it = 0;
	eps_sq=epsilon * epsilon;
	a[0]=a0=1.0/(zin+1); af=1.0/(zf+1);
	//Intervalo de de tempo da simulacao  para p p=a^alpha(parametrizado);
	p0=pow(a0,alpha); p=p0;  dp_min=(pow(af,alpha)-p0)/max_steps;


	setSaveSnapshotRedshift();
	z=1/a0-1;ek=0;
	double dTau=Compute_DTau();
	dtau_min=dTau/max_steps;
	dp_min=dp_min*eta;dp=dp_min;
	dtau_min=dtau_min*eta;	dtau=dtau_min;

    inicialize_leap_frog_CPU();
	if (sav_ener && !Benchmark) compute_energyCPU();
	save_snapshot();
	gettimeofday(&t0,NULL); // inicio da medida
	// Efeutar a simulacao
	it_old=0;			double c1,c2,c3;
int s=1;
	for (int i = 0; i < max_steps; i++) {
		it=i-it_old;   z=1/a[it]-1;

		if(sav_ener&&!Benchmark)integrate_edelta();
		if(Benchmark&&(i==1||i==101||i==201))gettimeofday(&tb0,NULL);

		da_dtau=sqrt(1+omegaM*(1/a[it]-1)+omegaR*(1/a[it]/a[it]-1)+omegaL*(a[it]*a[it]-1));
		d2a_dtau2=-da_dtau*da_dtau/a[it]*(omegaM/2-omegaL-omegaR);
				dp_dt=(alpha*p*da_dtau)/a[it];
				A=(1+alpha+d2a_dtau2*a[it]/(da_dtau*da_dtau))/(2*alpha*p);
				B=1/(dp_dt*H0*a[it]*a[it]*a[it]);
				compute_acc_stepsize(B,-2*A,eps_sq/a[it]/a[it]);
						printf("\n acc_max^2 %f vel_max^2 %f",d[0],d[1]);
                    	d[0]=eta/sqrt(sqrt(d[0]));				d[1]=eta/sqrt(d[1])/dp_dt/H0;
						printf(" d1 %f  d2 %f ",d[0],d[1]);
						dp=d[0]<d[1]?d[0]:d[1];
						dp=dp<dp_min?dp_min:dp;
						printf("\n Interacao: %d , a=%f z=%f p=%f dp=%f  A %f B %f 1/dp_dt %f \n", i,a[it],z,p,dp,A,B,dp_dt);
						c1=1;//vel
						c2=	dp/(1+A*dp);//acc
						c3=dp/dp_dt/H0;// vel
						integrate(c1,c2,c3);
				a[it+1]=pow(p,1.0/alpha);
		if(Benchmark&&(i==101||i==201||i==301)){
			gettimeofday(&tb1,NULL);// fim da medida
			printf("\n Tempo de procesamento de 100 passos %f secs \n",calcular_tempo(tb1,tb0));
			if (i==301) break;
		}

		if (((a[it] <=a_up) and (a[it] >=a_dow)and sav_snap and !Benchmark)or (a[it+1]>=(af))or(i==(max_steps-1))){
			gettimeofday(&t00,NULL);//inicio da medida
			save_particle_info();
			if (sav_ener) compute_energyCPU();
			save_snapshot();
			gettimeofday(&t01,NULL);// fim da medida
			tcom=tcom+calcular_tempo(t01,t00);
		}



		// Verificar se deve ser gravado um snapshot
		if (a[it+1]>=(af))break;
		else
			if (a[it]>a_up)
				{
				a_up=1.02*z_s[snap+s];
				a_dow=a_up/1.02*.98;
				s++;
				}
	}
	gettimeofday(&t1,NULL);// fim da medidaHtau[i]
	tproc=calcular_tempo(t1,t0);
	printf("\n Tempo de procesamento %f secs \n",tproc-tcom);
	printf("\n Tempo de comunicacao e calculo da energia em %f secs \n",tcom);
	save_particle_info();
    if (sav_ener)save_energy();

}

void Direct::DoGPUSimulation(){

	struct timeval t0,t1,t00,t01;
    //variable para el calculo del tiempo
	flip = false;
	allocar_memoriaCPU();
	printf("\n Particulas de materia escura %d de barions ",mNumBodies);
	// Copiar dados para a memoria de GPU
	allocar_memoriaGPU();


	gettimeofday(&t0,NULL); // inicio da medida
	copy_cputogpu();
	gettimeofday(&t1,NULL);// fim da medida
	tcom+=calcular_tempo(t1,t0);
	printf("\n Finalizado em %f secs \n",tcom);
	//constantes da GPU
	setSofteningValue(epsilon); setL(L);  setm(m); setmG(m*G);
	setisnap(i_snapg,ni_snap);

	setNumBodies(mNumBodies);   seteta(eta);
	set_address(Pos_g_a,Vel_g_a,Pos_g_b,Vel_g_b,acc_g,norm_acc_g,norm_vel_g,pot_g,d_g);
	set_addressi(Pos_ig0,Vel_ig0,Pos_ig1,Vel_ig1,Pos_ig2,Vel_ig2,Pos_ig3,Vel_ig3,Pos_ig4,Vel_ig4,Pos_ig5,Vel_ig5,Pos_ig6,Vel_ig6,Pos_ig7,Vel_ig7);
	eps_sq=epsilon*epsilon;
	snap=0;it=0;
	a[0]=1.0/(zin+1);		 af=1.0/(zf+1);
	a0=a[0];
	z=1/a[0]-1;ek=0.0;
	//Intervalo de de tempo da simulacao  para p p=a^alpha(parametrizado);
	p0=pow(a0,alpha); 
	p=p0;
	dp_min=(pow(af,alpha)-p)/max_steps;

	setSaveSnapshotRedshift();



	double dTau=Compute_DTau();
	dtau_min=dTau/max_steps;
	dp=dp_min;
	dtau=dtau_min;
	dp_min=dp_min*eta;
	dtau_min=dtau_min*eta;
			inicialize_leap_frog_GPU();
	copy_gputocpuP();
	if (sav_ener)compute_energyGPU();
	save_snapshot();
it_old=0;			double c1,c2,c3;
int s=1;
	gettimeofday(&t0,NULL); // inicio da medida
		for (int i = 0; i < max_steps; i++) {
			it=i-it_old;
			if(sav_ener )integrate_edelta();
			da_dtau=sqrt(1+omegaM*(1/a[it]-1)+omegaR*(1/a[it]/a[it]-1)+omegaL*(a[it]*a[it]-1));
			d2a_dtau2=-da_dtau*da_dtau/a[it]*(omegaM/2-omegaL-omegaR);
			z=1/a[it]-1;
					Htau=da_dtau/a[it];
							c1=1/H0/a[it]/a[it];//sum
							c2=	-Htau;//vel
							printf(" c1 %f  c2 %f ",c1,c2);
							g_compute_acc(c1,c2,eps_sq/a[it]/a[it],mNumBodies);cudaThreadSynchronize();
							g_compute_stepsize(mNumBodies);
							cudaThreadSynchronize();
							cudaMemcpy(d,d_g,sizeof(double)*2, cudaMemcpyDeviceToHost);cudaThreadSynchronize();
							printf("\n acc_max^2 %f vel_max^2 %f",d[0],d[1]);
							d[0]=eta*sqrt(a[it]/sqrt(d[0]))*a[it];               d[1]=eta/sqrt(d[1]);
								d[0]=d[0]/sqrt(a[it]);
								d[1]=d[1]*sqrt(a[it]);
							printf(" d1 %f  d2 %f ",d[0],d[1]);
							dtau=d[0]<d[1]?d[0]:d[1];				dtau=dtau<dtau_min?dtau_min:dtau;
							c1=1.0;
							c2=dtau/(1+Htau/2*dtau);//acc
							c3=dtau/H0/a[it];// vel
						g_integrate(it,c1,c2,c3,mNumBodies);

						//	printf("\n c1=%f c2=%f c3=%f ", c1,c2,c3);
					a[it+1]=a[it]+da_dtau*dtau+d2a_dtau2*(dtau*dtau)/2;
					printf("\n Interacao: %d it= %d, a=%f a[it+1]=%f z=%f p=%f dtau_min %f dtau=%f \n", i,it,a[it],a[it+1],z,p,dtau_min,dtau);

			tau[it+1]=tau[it+1]+dtau;
			cudaThreadSynchronize();

			g_save_isnap(ni_snap,it);

			// Verificar se deve ser gravado um snapshot
			cudaThreadSynchronize();

			if(!flip)
				flip_pos(Pos_g_b,Vel_g_b,Pos_g_a,Vel_g_a);
			else
				flip_pos(Pos_g_a,Vel_g_a,Pos_g_b,Vel_g_b);
			flip = !flip;

			if ((a[it] <=a_up) and (a[it] >=a_dow)and(sav_snap)or a[it+1]>=(af)or(i==(max_steps-1))){
					gettimeofday(&t00,NULL); // inicio da medida
					copy_gputocpu_particules();
					save_particle_info();
					copy_gputocpuP();
			       	if (sav_ener)compute_energyGPU();
					save_snapshot();
					gettimeofday(&t01,NULL);// fim da medida
					tcom=tcom+calcular_tempo(t01,t00);
					}

			if (a[it+1]>=af) break;
			else
					if (a[it]>a_up)
						{
						a_up=1.02*z_s[snap+s];
						a_dow=a_up/1.02*.98;
						s++;
						}
			//Salir no redshift sufientemente pequenho
			// Inverter variaveis
		}

		gettimeofday(&t1,NULL);// fim da medida
		tproc=calcular_tempo(t1,t0)-tcom;
		printf("\n Finalizado Tempo da gpu %f secs \n",tproc);
		printf("\n Tempo de comunicacao+energia %f secs \n",tcom);
		if (sav_ener)save_energy();
		liberar_memoria_gpu();



}

void Direct::allocar_memoriaCPU()
{
	printf("\n Allocate memory RAM...");
	if(sav_ener)
		{
		pot = (double*) malloc(sizeof(double)*mNumBodies);
		ekin = (double*) malloc(sizeof(double)*(max_snapshots+4));
		epot = (double*) malloc(sizeof(double)*(max_snapshots+4));
		etotal = (double*) malloc(sizeof(double)*(max_snapshots+4));
		edelta = (double*) malloc(sizeof(double)*(max_snapshots+4));
		if (pot==NULL or ekin==NULL or epot==NULL or etotal==NULL or edelta==NULL ){
		printf("\n Not enough memory for save snapshots");
	    exit (1);}
		edelta[0]=0.0;
		}
	z_s = (double*) malloc(sizeof(double)*(max_snapshots+4));
    buf = (char*) malloc(sizeof(char)*(50));
	a = (double*) malloc(sizeof(double)*(int)(max_steps/5));
	tau = (double*) malloc(sizeof(double)*(int)(max_steps/5));


if (ni_snap>0){
	Pos_i= (double3**) malloc(sizeof(double3*)*ni_snap);
	Vel_i= (double3**) malloc(sizeof(double3*)*ni_snap);
  for (int i = 0; i < (ni_snap); i++) {
	  Pos_i[i]= (double3*) malloc(sizeof(double3)*(int)max_steps/5);
	  Vel_i[i]= (double3*) malloc(sizeof(double3)*(int)max_steps/5);
	}
}
	if ( z_s==NULL or a==NULL or buf==NULL )
	{
	printf("\n Not enough memory");
		exit (1);
	}
	printf("...done\n");
}
void Direct::liberar_memoriaCPU()
{
	printf("\n Free RAM...");
	if(sav_ener)
		{
		free(pot);
		free(ekin);
		free(epot);
		free(etotal);
		free(edelta);
		}

//free(z_s);
	free(a);
	free(tau);
	free(buf);
	if (ni_snap>0){
	for (int i = 0; i < (ni_snap); i++) {
	free(Pos_i[i]);
	free(Vel_i[i]);
	}
	free(Pos_i);
	free(Vel_i);
	}
	printf("...done\n");
}
void Direct::allocar_memoriaGPU()
{
	printf("\n Allocate GPU memory ");
	cudaMalloc(&d_g, sizeof(double)*2);
	cudaMalloc(&norm_acc_g, sizeof(double) * mNumBodies);
	cudaMalloc(&norm_vel_g, sizeof(double) * mNumBodies);
	cudaMalloc(&acc_g, sizeof(double3) * mNumBodies);
	cudaMalloc(&Pos_g_a, sizeof(double3) * mNumBodies);
	cudaMalloc(&Pos_g_b, sizeof(double3) * mNumBodies);
	cudaMalloc(&Vel_g_a, sizeof(double3) * mNumBodies);
	cudaMalloc(&Vel_g_b, sizeof(double3)* mNumBodies);


	int nn=(int)max_steps/5;
	cudaMalloc(&Pos_ig0, sizeof(double3*)*nn);
	cudaMalloc(&Vel_ig0, sizeof(double3*)*nn);
	cudaMalloc(&Pos_ig1, sizeof(double3*)*nn);
	cudaMalloc(&Vel_ig1, sizeof(double3*)*nn);
	cudaMalloc(&Pos_ig2, sizeof(double3*)*nn);
	cudaMalloc(&Vel_ig2, sizeof(double3*)*nn);
	cudaMalloc(&Pos_ig3, sizeof(double3*)*nn);
	cudaMalloc(&Vel_ig3, sizeof(double3*)*nn);
	cudaMalloc(&Pos_ig4, sizeof(double3*)*nn);
	cudaMalloc(&Vel_ig4, sizeof(double3*)*nn);
	cudaMalloc(&Pos_ig5, sizeof(double3*)*nn);
	cudaMalloc(&Vel_ig5, sizeof(double3*)*nn);
	cudaMalloc(&Pos_ig6, sizeof(double3*)*nn);
	cudaMalloc(&Vel_ig6, sizeof(double3*)*nn);
	cudaMalloc(&Pos_ig7, sizeof(double3*)*nn);
    cudaMalloc(&Vel_ig7, sizeof(double3*)*nn);
	cudaMalloc(&i_snapg, sizeof(int)* ni_snap);
	if(sav_ener)cudaMalloc(&pot_g, sizeof(double) * mNumBodies);
	printf("...done\n");
}
void Direct::Bodyinteration(int i,double eps){
	sum_mGrji_dij2=make_double3(0,0,0);
	#pragma omp parallel for
	for (int j = 0; j < mNumBodies; j++){
		if(i!=j){
				double3 r;
					r.x = fmod(Pos[j].x+L*1.5-Pos[i].x,L)-L*0.5;
					r.y = fmod(Pos[j].y+L*1.5-Pos[i].y,L)-L*0.5;
					r.z = fmod(Pos[j].z+L*1.5-Pos[i].z,L)-L*0.5;
					double d_rij = sqrt(r.x*r.x + r.y*r.y + r.z*r.z +eps);
					sum_mGrji_dij2.x += m*G*r.x /d_rij/d_rij/d_rij;
					sum_mGrji_dij2.y += m*G*r.y /d_rij/d_rij/d_rij;
					sum_mGrji_dij2.z += m*G*r.z /d_rij/d_rij/d_rij;
		}
		}
}
void Direct::save_i_snap() {
	#pragma omp parallel for
	for (int i = 0; i < (ni_snap); i++) {
    	Pos_i[i][it]=Pos[i_snap[i]];
	    Vel_i[i][it]=Vel[i_snap[i]];
	}
}
void Direct::compute_acc(double c1,double c2,double eps_sq_c) {
double acc_max=0.0,n_acc,vel_max=0.0,n_vel;
#pragma omp parallel for
for (int i = 0; i < mNumBodies; i++){
	Bodyinteration(i,eps_sq_c);
	acc[i].x =sum_mGrji_dij2.x*c1+Vel[i].x*c2;
	acc[i].y =sum_mGrji_dij2.y*c1+Vel[i].y*c2;
	acc[i].z =sum_mGrji_dij2.z*c1+Vel[i].z*c2;
}

}
void Direct::compute_acc_stepsize(double c1,double c2,double eps_sq_c) {
double acc_max=0.0,n_acc,vel_max=0.0,n_vel;
#pragma omp parallel for
for (int i = 0; i < mNumBodies; i++){
	Bodyinteration(i,eps_sq_c);
	acc[i].x =sum_mGrji_dij2.x*c1+Vel[i].x*c2;
	acc[i].y =sum_mGrji_dij2.y*c1+Vel[i].y*c2;
	acc[i].z =sum_mGrji_dij2.z*c1+Vel[i].z*c2;
	n_acc=(acc[i].x*acc[i].x+acc[i].y*acc[i].y+acc[i].z*acc[i].z);
	n_vel=(Vel[i].x*Vel[i].x+Vel[i].y*Vel[i].y+Vel[i].z*Vel[i].z);
	if (n_acc > acc_max )
		acc_max=n_acc;
	if (n_vel > vel_max )
		vel_max=n_vel;

}
d[0]=acc_max;
d[1]=vel_max;

}
void Direct::integrate_lf_1_2(double c1,double c2,double c3) {
	double3 vel;
	#pragma omp parallel for
	for (int i = 0; i < mNumBodies; i++) {
		vel=Vel[i];
		Vel[i].x=c1*Vel[i].x;
		Vel[i].y=c1*Vel[i].y;
		Vel[i].z=c1*Vel[i].z;
		vel.x = vel.x*c1+ acc[i].x* c2;
		vel.y = vel.y*c1+ acc[i].y* c2;
		vel.z = vel.z*c1+ acc[i].z* c2;
		Pos[i].x = fmod(Pos[i].x+vel.x* c3+L,L);
		Pos[i].y = fmod(Pos[i].y+vel.y* c3+L,L);
		Pos[i].z = fmod(Pos[i].z+vel.z* c3+L,L);

	}
	save_i_snap();

}
void Direct::integrate(double c1,double c2,double c3) {
	#pragma omp parallel for
	for (int i = 0; i < mNumBodies; i++) {
		    Vel[i].x = Vel[i].x*c1 +acc[i].x*c2;
		    Vel[i].y = Vel[i].y*c1 +acc[i].y*c2;
		    Vel[i].z = Vel[i].z*c1 +acc[i].z*c2;
	    	Pos[i].x = fmod (Pos[i].x+Vel[i].x*c3+L,L);
	    	Pos[i].y = fmod (Pos[i].y+Vel[i].y*c3+L,L);
	    	Pos[i].z = fmod (Pos[i].z+Vel[i].z*c3+L,L);
	}
	if(track_on)save_i_snap();
}
void Direct::inicialize_leap_frog_GPU() {
double c1,c2,c3;it=0;
	printf("\n Inicialize Leap_frog GPU...");
	da_dtau=sqrt(1+omegaM*(1/a[it]-1)+omegaR*(1/a[it]/a[it]-1)+omegaL*(a[it]*a[it]-1));
	d2a_dtau2=-da_dtau*da_dtau/a[it]*(omegaM/2-omegaL-omegaR);

	       Htau=da_dtau/a[it];
			c1=	1/H0/a[it]/a[it];
			c2=-Htau*a[it];
			printf("\n c1 %f  c2%f ",c1,c2);

			//g_compute_acc_andtimestep(c1,c2,eps_sq/a[it]/a[it],mNumBodies);
			g_compute_acc(c1,c2,eps_sq/a[it]/a[it],mNumBodies);
				cudaThreadSynchronize();

				g_compute_stepsize(mNumBodies);
				cudaThreadSynchronize();
				cudaMemcpy(d,d_g,sizeof(double)*2, cudaMemcpyDeviceToHost);
				cudaThreadSynchronize();
			    printf("\nv2 acc_max^2 %f vel_max^2 %f",d[0],d[1]*a[it]*a[it]);
				d[0]=eta*sqrt(a[it]/sqrt(d[0]))*a[it]/sqrt(a[it]);               d[1]=eta/sqrt(d[1])*sqrt(a[it]);
				printf(" d1 %f  d2 %f ",d[0],d[1]);
				dtau=d[0]<d[1]?d[0]:d[1];
				dtau=dtau<dtau_min?dtau_min:dtau;

			printf("\n Interacao: %d , a=%f z=%f p=%f dtau=%f \n", it,a[it],z,p,dtau);
			c1=a[it];
			c2=dtau/4;//vel
			c3=dtau/H0/2/a[it];//acc
			g_integrate_lf_1_2(it,c1,c2,c3,mNumBodies);

			//		  printf("\n c1 %f  c2%f ",c1,c2);
			a[it]=a[it]+da_dtau*dtau/2+d2a_dtau2*(dtau*dtau)/2/4;

	tau[it]=dtau/2;
	cudaThreadSynchronize();
	g_save_isnap(ni_snap,it);
	cudaThreadSynchronize();

	flip_pos(Pos_g_b,Vel_g_b,Pos_g_a,Vel_g_a);
	flip = !flip;

	printf("...done\n");
}
void Direct::inicialize_leap_frog_CPU() {
	double c1,c2,c3;
it=0;
	printf("\n Inicialize Leap_frog...");
	da_dtau=sqrt(1+omegaM*(1/a[it]-1)+omegaR*(1/a[it]/a[it]-1)+omegaL*(a[it]*a[it]-1));
	d2a_dtau2=-da_dtau*da_dtau/a[it]*(omegaM/2-omegaL-omegaR);
				Htau=da_dtau/a[it];

				c1=	1/H0/a[it]/a[it];
				c2=-Htau*a[it];
					printf("\n c1 %f  c2%f ",c1,c2);
				compute_acc_stepsize(c1,c2,eps_sq/a[it]/a[it]);
			printf("\n v acc_max^2 %f vel_max^2 %f",d[0],d[1]*a[it]*a[it]);
				d[0]=eta*sqrt(a[it]/sqrt(d[0]))*a[it]/sqrt(a[it]);               d[1]=eta/sqrt(d[1])*sqrt(a[it]);
				printf(" d1 %f  d2 %f ",d[0],d[1]);
				dtau=d[0]<d[1]?d[0]:d[1];
				dtau=dtau<dtau_min?dtau_min:dtau;
			printf("\n Interacao: %d , a=%f z=%f p=%f dtau=%f \n", it,a[it],z,p,dtau);
			c1=a[it];
			c2=dtau/4;//vel
			c3=dtau/H0/2/a[it];//acc
			integrate_lf_1_2(c1,c2,c3);
			printf("\n c1 %f  c2%f ",c1,c2);

			a[it]=a[it]+da_dtau*dtau/2+d2a_dtau2*(dtau*dtau)/2/4;
	tau[it]=dtau/2;

printf("...done\n");
}
double Direct::EKin(){
	double _ekin = 0.0;
	for (int i = 0; i < mNumBodies; i++)
		_ekin +=  (Vel[i].x*Vel[i].x+Vel[i].y*Vel[i].y+Vel[i].z*Vel[i].z);
	return 0.5*_ekin/mNumBodies;
}
double Direct::EPot(){
	double _epot = 0.0;
	for (int i = 0; i < mNumBodies; i++) {
		#pragma omp parallel for
		for (int j = 0; j < mNumBodies; j++){
			if (i != j){
		      double3 r;
		      double L_2=L/2;
				r.x = fmod(Pos[j].x+L_2-Pos[i].x,L)-L_2;
				r.y = fmod(Pos[j].y+L_2-Pos[i].y,L)-L_2;
				r.z = fmod(Pos[j].z+L_2-Pos[i].z,L)-L_2;
				double d = (r.x*r.x+r.y*r.y+r.z*r.z);
				_epot += -m*G / sqrt(d);
			}
		}

	}
	return 0.5*_epot/mNumBodies;
}
void Direct::integrate_edelta(){
	double ek_1=ek;
	printf("\n Integral a*aWda...");

	if(use_gpu)
	{
		printf(" Kinetic GPU");
		g_EKin(mNumBodies);
		cudaThreadSynchronize();
		cudaMemcpy(pot, pot_g,sizeof(double) * mNumBodies, cudaMemcpyDeviceToHost);
		cudaThreadSynchronize();
		ek=0.0;
		#pragma omp parallel for
		for (int j = 0; j < mNumBodies; j++)
			ek += pot[j];
    	ek=0.5*ek/mNumBodies;
    	printf("..it %d snap %d..%f",it,snap,ek);
	}
	else
	{
		printf(" Kinetic CPU");
		ek=EKin();
		printf("..it %d snap %d..%f",it,snap,ek);

	}

	if(it>=1)
	 edelta[snap]+=(a[it]*a[it]*ek+a[it-1]*a[it-1]*ek_1)*(a[it]-a[it-1])/2;

	printf("...done");

}
void Direct::compute_energyGPU(){
	double ek_1=ek,edelta_it=0.0;
	printf("\n Computing energy %d ...",snap);
	printf("Potential gpu ...");
	g_EPot(mNumBodies);
	cudaThreadSynchronize();
	cudaMemcpy(pot, pot_g,sizeof(double) * mNumBodies, cudaMemcpyDeviceToHost);
	cudaThreadSynchronize();
	epot[snap]=0.0;
	#pragma omp parallel for
	for (int j = 0; j < mNumBodies; j++)
		epot[snap] += pot[j];
	epot[snap] =-0.5*epot[snap]/mNumBodies ;
	printf("...ok");
	ek=EKin();

    printf("..kin %f  pot %f",ek,epot[snap]);
	ekin[snap]=ek;
	edelta[snap+1]=edelta[snap];
	if(snap==0)	a[it+1]=a[it];
	ekin[snap] = a[it+1]*a[it+1]*a[it+1]*ek;
	if(snap>0)edelta_it=(a[it+1]*a[it+1]*ek+a[it]*a[it]*ek_1)*(a[it+1]-a[it])/2;
	edelta[snap]+=edelta_it;	ek=ek_1;

	etotal[snap] = epot[snap]+ekin[snap]+edelta[snap];
	printf("\n snap %d W %f +T %f =  eTotal %f etotal+1 %f",snap,ekin[snap],epot[snap],etotal[snap],edelta[snap+1]);
	printf("...done \n");

	}
void Direct::compute_energyCPU(){
	printf("\n Compute energy %d ...",snap);
	double ek_1=ek,edelta_it=0.0;
		ek=EKin();
		ekin[snap]=ek;
	epot[snap] = EPot();
	edelta[snap+1]=edelta[snap];
	if(snap==0)	a[it+1]=a[it];
	printf("..kin %f  pot %f  a %f",ek,epot[snap],a[it+1]);

    ekin[snap] = a[it+1]*a[it+1]*a[it+1]*ek;
	if(snap>0)edelta_it=(a[it+1]*a[it+1]*ek+a[it]*a[it]*ek_1)*(a[it+1]-a[it])/2;
	edelta[snap]+=edelta_it; ek=ek_1;
	etotal[snap] = epot[snap]+ekin[snap]+edelta[snap];
	printf("\n snap %d W %f +T %f =  eTotal %f etotal+1 %f",snap,ekin[snap],epot[snap],etotal[snap],edelta[snap+1]);

	printf("...done\n");

	}


void Direct::save_snapshot(){
FILE *F;
	printf("\n Saving snapshot %d ....",snap);
	a_up=1.02*z_s[snap+1];
	a_dow=a_up/1.02*.98;
    z_s[snap]=z;
	if(snap<(max_snapshots))printf("\n  Next %f %f",1/a_up-1,1/a_dow-1);
				// Arquivo
	double3 vel;
	 sprintf(buf,"%ssnap%d.dat",mFileOut.c_str(),snap);
	 			if (mBinOutput) f = fopen( buf, "wb");
				else f = fopen(  buf, "w");
         		if (!mBinOutput)
					for (int i = 0; i < mNumBodies; i++){
						vel.x=10*Vel[i].x;	            vel.y=10*Vel[i].y;                 vel.z=10*Vel[i].z;
						fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n", m,Pos[i].x, Pos[i].y, Pos[i].z,vel.x,vel.y,vel.z);
					}
				else
				{
					printf("...in binary...");
					for (int i = 0; i < mNumBodies; i++) {
						vel.x=10*Vel[i].x;	 vel.y=10*Vel[i].y;        vel.z=10*Vel[i].z;
						fwrite((void*)&m, sizeof(double), 1, f);
						fwrite((void*)&Pos[i], sizeof(double3), 1, f);
						fwrite((void*)&vel, sizeof(double3), 1, f);
					}

					printf("...done\n");
				}
				fclose(f);

				if (((a[it+1]>=(af))or(it+it_old==(max_steps-1))) and mBinOutput)
				{
				printf("...last in ascci...");
				 sprintf(buf,"%ssnap%d.dat.ascii",mFileOut.c_str(),snap);
					F = fopen(  buf, "w");
					for (int i = 0; i < mNumBodies; i++){
						vel.x=10*Vel[i].x;	 vel.y=10*Vel[i].y;        vel.z=10*Vel[i].z;
						fprintf(F, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n", m,Pos[i].x, Pos[i].y, Pos[i].z,vel.x,vel.y, vel.z);
					}
					fclose(F);

				}
				printf("...done\n");
				snap++;
			}
void Direct::save_particle_info(){
	printf("\n Saving tracking info (%d particules.. max %d)...",it,(int)(max_steps/5));
	if((it)>(max_steps/5))	printf("\n Error: Allocate more memory for scale factor");

	double3 vel;int st=1;
	sprintf(buf, "%s.i.dat",mFileOut.c_str());

	if (snap==1){
		f = fopen(buf, "w");
		st=0;
	}
	else
		f = fopen(buf, "a");
	for (int i = st; i <= it; i++)
	{
		fprintf(f, "%d %f %f %f ",i,tau[i],a[i],1.0/a[i]-1.0);
		for (int j = 0; j < (ni_snap); j++) {
		fprintf(f, "%d ",j);

				vel.x=10*Vel_i[j][i].x;	 vel.y=10*Vel_i[j][i].y;        vel.z=10*Vel_i[j][i].z;
			fprintf(f, "%f %f %f %f %f %f ",Pos_i[j][i].x,Pos_i[j][i].y,Pos_i[j][i].z,vel.x,vel.y,vel.z);
		}
		fprintf(f, "\n");

	}
	fclose(f);
	printf("\n...scale factor...");


	sprintf(buf, "%s.a.dat",mFileOut.c_str());
	if (snap==1){
	  f = fopen(buf, "w");
	}
	else
	  f = fopen(buf, "a");

	for (int i = st; i <= it; i++)
	{
		fprintf(f, "%d  %f %f \n",i+it_old,tau[i],a[i],(1.0/a[i]-1.0));
	}

	fclose(f);
		printf("...done\n");
it_old=it_old+it;
a[0]=a[it];
a[1]=a[it+1];
it=0;
}
void Direct::save_energy(){

	if (sav_ener)
			{
				printf("\n Saving energy...");
				 sprintf(buf, "%s.energy.dat",mFileOut.c_str());
				f = fopen(buf, "w");
				for (int i = 0; i < snap; i++){
						fprintf(f, "%f %f %f %f %f \n",ekin[i], epot[i],edelta[i],etotal[i],z_s[snap]);
//						fprintf(f, "%f %f %f %f \n",ekin[i], epot[i],etotal[i],z_s[snap]);

				}

			fclose(f);
				printf("...done\n");
			}
}
void Direct::copy_gputocpuP(){
	printf("\n Copiando da GPU para CPU snapshot %d iteracao%d z= %f...",snap,it,z);
	if (flip){
		cudaMemcpy(Pos,Pos_g_b,sizeof(double3) * mNumBodies, cudaMemcpyDeviceToHost);
		cudaMemcpy(Vel,Vel_g_b,sizeof(double3) * mNumBodies, cudaMemcpyDeviceToHost);
	} else {
		cudaMemcpy(Pos,Pos_g_a,sizeof(double3) * mNumBodies, cudaMemcpyDeviceToHost);
		cudaMemcpy(Vel,Vel_g_a,sizeof(double3) * mNumBodies, cudaMemcpyDeviceToHost);
	}
	cudaThreadSynchronize();
	printf("...done\n");

}
void Direct::copy_cputogpu(){
	printf("\n Copy CPU to GPU ...");
		cudaMemcpy(Vel_g_a, Vel, sizeof(double3) * mNumBodies, cudaMemcpyHostToDevice);
		cudaMemcpy(Pos_g_a, Pos, sizeof(double3) * mNumBodies, cudaMemcpyHostToDevice);
		cudaMemcpy(Pos_g_a, Pos, sizeof(double3) * mNumBodies, cudaMemcpyHostToDevice);
		cudaMemcpy(i_snapg,i_snap, sizeof(int) * ni_snap, cudaMemcpyHostToDevice);


		printf("...done\n");
}

void Direct::copy_gputocpu_particules()
{
	printf("\n Copy GPU to CPU(particules) ...");
	cudaMemcpy(Pos_i[0],Pos_ig0,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaMemcpy(Vel_i[0],Vel_ig0,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);

	cudaMemcpy(Pos_i[1],Pos_ig1,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaMemcpy(Vel_i[1],Vel_ig1,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);

	cudaMemcpy(Pos_i[2],Pos_ig2,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaMemcpy(Vel_i[2],Vel_ig2,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);

	cudaMemcpy(Pos_i[3],Pos_ig3,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaMemcpy(Vel_i[3],Vel_ig3,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);

	cudaMemcpy(Pos_i[4],Pos_ig4,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaMemcpy(Vel_i[4],Vel_ig4,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);

	cudaMemcpy(Pos_i[5],Pos_ig5,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaMemcpy(Vel_i[5],Vel_ig5,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);

	cudaMemcpy(Pos_i[6],Pos_ig6,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaMemcpy(Vel_i[6],Vel_ig6,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);

	cudaMemcpy(Pos_i[7],Pos_ig7,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaMemcpy(Vel_i[7],Vel_ig7,sizeof(double3) * (it+1), cudaMemcpyDeviceToHost);
	cudaThreadSynchronize();
	printf("...done\n");

}
void Direct::liberar_memoria_gpu()
{
	printf("\n Free gpu memory ...");

	cudaFree(d_g);
	cudaFree(norm_acc_g);
	cudaFree(norm_vel_g);
	cudaFree(acc_g);
	cudaFree(Pos_g_a);
	cudaFree(Vel_g_a);
	cudaFree(Pos_g_b);
	cudaFree(Vel_g_b);


	cudaFree(Vel_ig0);
	cudaFree(Pos_ig0);

	cudaFree(Vel_ig1);
	cudaFree(Pos_ig1);

	cudaFree(Vel_ig2);
	cudaFree(Pos_ig2);

	cudaFree(Vel_ig3);
	cudaFree(Pos_ig3);

	cudaFree(Vel_ig4);
	cudaFree(Pos_ig4);

	cudaFree(Vel_ig5);
	cudaFree(Pos_ig5);

	cudaFree(Vel_ig6);
	cudaFree(Pos_ig6);

	cudaFree(Vel_ig7);
	cudaFree(Pos_ig7);


	if(sav_ener)cudaFree(pot_g);

	printf("...done\n");
}

double Direct::calcular_tempo(struct timeval t1, struct timeval t0){
	double t_sec,t_usec,sec;
	t_sec  = (double)  (t1.tv_sec - t0.tv_sec);
	t_usec = (double)  (t1.tv_usec - t0.tv_usec);
	sec=t_sec + t_usec/1.0e+6;
	printf("\n Finalizado em %f secs \n",sec);

	return sec;
}
double Direct::Compute_DTau()
{
	double result, abserr;
	para par;
	par.omegaM=omegaM;
	par.omegaL=omegaL;
	par.omegaR=omegaR;

	void *params_ptr=&par;
	gsl_integration_workspace *workspace;
	gsl_function F;

	workspace = gsl_integration_workspace_alloc(WORKSIZE);

	F.function = &fun;
	F.params = params_ptr;
	gsl_integration_qag(&F, 1/(zin+1), 1/(zf+1), 0, 1.0e-12, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

	gsl_integration_workspace_free(workspace);

	return (double)result;
}

extern "C"{
double fun(double a ,void* par)
		{
	    para p=*(para*) par;
		return 1/a/sqrt(p.omegaM / (a * a * a) + p.omegaR / (a * a * a * a) + p.omegaL);
		}
}
}


}


