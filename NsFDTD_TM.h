#ifndef _NS_FDTD_TM_
#define _NS_FDTD_TM_
#include "FDTD_TM.h"

class NsFDTD_TM: public FDTD_TM{
	typedef FDTD_TM super;
	double R_P, R_M;
public:
	NsFDTD_TM();
	~NsFDTD_TM();
	bool calc();
	void field();
	void PMLfield();
private:

	complex<double> Dx2_n(complex<double> *p, int i, int j, int t){		//差分演算子d'2の成分
		return (p[index(i,j+1, t)] + p[index(i, j-1, t)] - p[index(i-1,j+1, t)] - p[index(i-1, j-1, t)])/2.0;
	};

	complex<double> Dy2_n(complex<double> *p, int i, int j, int t){
		return (p[index(i+1,j, t)] + p[index(i-1, j, t)] - p[index(i+1,j-1, t)] - p[index(i-1, j-1, t)])/2.0;
	}

	complex<double> Dx0_n(complex<double> *p, int i, int j, int t) {	//差分演算子d'0の成分
		complex<double> dx = p[index(i+1, j, t)] - p[index(i, j, t)];
		return (dx + R_M/2 * DxDy2(p, i, j, 0));
	};

	complex<double> Dy0_n(complex<double> *p, int i, int j, int t) {
		complex<double> dy = p[index(i, j+1, t)] - p[index(i, j, t)];
		return (dy + R_M/2 * DxDy2(p, i, j, 0));
	};

	//todo 境界近傍にS-FDTDを使ってない
	void CalcE() {
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
		//電界の計算
		for (int i = 1; i < mField->getNpx() - 1; i++) {
			for (int j = 1; j < mField->getNpy() - 1; j++) {
				EZ(i, j, +1) = CEZ(i, j)*EZ(i, j, 0)
					+ CEZLH(i, j)*(HY(i, j, 0) - HY(i - 1, j, 0) - (HX(i, j, 0) - HX(i, j - 1, 0))	//dx(Hy) - dy(Hy)
						);
			}
		}
	}


	//吸収境界はEzにしか適用しないから,Hは領域の端も普通に計算する(できる分は)
	void CalcH(){	//todo 計算領域 i=0?, 1? j=0?, 1?
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
		{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=1; i<mField->getNpx()-1; i++){
			for(int j=0; j<mField->getNpy()-1; j++){
				HX(i, j, +1) = HX(i, j, 0) - CHXLYNS(i, j)*(R_P*(EZ(i, j + 1, +1) - EZ(i, j, +1)) + R_M*Dy2_n(Ez, i, j + 1, +1));	//Hxの計算 Hx(i, j+1/2) -> Hx[i,j]
			}
		}
				
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=0; i<mField->getNpx()-1; i++){
			for(int j=1; j<mField->getNpy()-1; j++){
				HY(i, j, +1) = HY(i, j, 0) + CHYLXNS(i, j)*(R_P*(EZ(i + 1, j, +1) - EZ(i, j, +1)) + R_M*Dx2_n(Ez, i + 1, j, +1));	//Hyの計算 Hy(i+1/2, j) -> Hy[i,j]
			}
		}
			
		}
	};

	void absorbing();	//吸収境界

	// Ns_PML
	void CalcE_PML(){		
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
		{
		#ifdef _OPENMP
		#pragma omp for
		#endif
			for (int i = 1; i < mField->getNpx() - 1; i++) {
				for (int j = 1; j < mField->getNpy() - 1; j++) {
					EZX(i, j, +1) = BEZXP(i, j)*BEZXM(i, j)*EZX(i, j, 0)
						+ BEZXP(i, j)*CEZLH(i, j)*(HY(i, j, 0) - HY(i - 1, j, 0));
				}
			}

		#ifdef _OPENMP
		#pragma omp for
		#endif
			for (int i = 1; i < mField->getNpx() - 1; i++) {
				for (int j = 1; j < mField->getNpy() - 1; j++) {
					EZY(i, j, +1) = BEZYP(i, j)*BEZYM(i, j)*EZY(i, j, 0)
						- BEZYP(i, j)*CEZLH(i, j)*(HX(i, j, 0) - HX(i, j - 1, 0));
				}
			}

		#ifdef _OPENMP
		#pragma omp for
		#endif
			for (int i = 1; i < mField->getNpx() - 1; i++)
				for (int j = 1; j < mField->getNpy() - 1; j++)
					EZ(i, j, +1) = EZX(i, j, +1) + EZY(i, j, +1);
		}
	}

	void CalcH_PML(){	//一部Ez=0を参照しているが、配列定義外参照はないからStは使わなくていい？
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
		{
		#ifdef _OPENMP
		#pragma omp for
		#endif
			for (int i = 1; i<mField->getNpx()-1; i++) {
				for (int j = 0; j<mField->getNpy() - 1; j++) {
//					if (i == 1 || i == mField->getNpx() - 2 || j == 0 || j == mField->getNpy() - 2)
//						HX(i, j, +1) = CHX(i, j)*HX(i, j, 0)
//						- CHXLY(i, j)*(EZX(i, j + 1, +1) - EZX(i, j, +1) + EZY(i, j + 1, +1) - EZY(i, j, +1));
//					else
						HX(i, j, +1) = BHYP(i, j) * BHYM(i, j) * HX(i, j, 0) 
							- BHYP(i, j)*CHXLYNS(i, j)*(R_P*(EZ(i, j+1, +1) - EZ(i, j, +1)) + R_M*Dy2_n(Ez, i, j+1, +1));	//Hxの計算 Hx(i, j+1/2) -> Hx[i,j]
													 //*Dy0_n(Ez, i, j, +1);
				}
			}

		#ifdef _OPENMP
		#pragma omp for
		#endif
			for (int i = 0; i<mField->getNpx() - 1; i++) {
				for (int j = 1; j<mField->getNpy()-1; j++) {
//					if (j == 1 || i == mField->getNpy() - 2 || i==0 || i==mField->getNpx()-2)
//						HY(i, j, +1) = CHY(i, j)*HY(i, j, 0)
//						+ CHYLX(i, j)*(EZX(i + 1, j, +1) - EZX(i, j, +1) + EZY(i + 1, j, +1) - EZY(i, j, +1));
//					else
						HY(i, j, +1) = BHXP(i, j) * BHXM(i, j) * HY(i, j, 0) 
							+ BHXP(i, j)*CHYLXNS(i, j)*(R_P*(EZ(i+1, j, +1) - EZ(i, j, +1)) + R_M*Dx2_n(Ez, i+1, j, +1));	//Hyの計算 Hy(i+1/2, j) -> Hy[i,j]
													 //*Dx0_n(Ez, i, j, +1);
				}
			}

		}
	};

	bool EndTask();		//1回のシミュレーションが終わったときの処理

	void ReStart(){
		super::Initialize();
		field();
	}
	void printDebugCalcE(int,int,int,int);
	void printDebugCalcHx(int,int,int,int);
	void printDebugCalcHy(int,int,int,int);
};

//プリントデバッグ
inline void NsFDTD_TM::printDebugCalcE(int i, int j, int _i, int _j){
	if(i != _i || j!=_j) return;
	cout << "===========pringDebug CalcE TM Mode===============" << endl; 
	cout << "(" + to_s(i) + "," + to_s(j) + ")" << endl;
	cout << "EZ(i,j,+1) = " << EZ(i,j,+1) << endl;
	cout << "EZ(i,j, 0) = " << EZ(i,j,0) << endl;
	cout << "EZ(i,j,-1) = " << EZ(i,j,-1) << endl;
	cout << "CE(i,j) = " << CEZ(i,j) << endl;
	cout << "CEZLH(i,j) = " <<  CEZLH(i,j) << endl;
	cout << "HY(i  ,j  , 0) = " << HY(i,j,0) << endl;
	cout << "HY(i-1,j  , 0) = " << HY(i-1,j,0) << endl;
	cout << "HY(i  ,j+1, 0) = " << HY(i,j+1,0) << endl;
	cout << "HY(i  ,j-1, 0) = " << HY(i,j-1,0) << endl;
	cout << "HY(i-1,j+1, 0) = " << HY(i-1,j+1,0) << endl;
	cout << "HY(i-1,j-1, 0) = " << HY(i-1,j-1,0) << endl;

	cout << "HX(i,j,0) = " << HX(i,j,0) << endl;

	cout << "HX(i  ,j-1,0) = " << HY(i,j-1,0) << endl;
	cout << "HX(i+1,j  ,0) = " << HX(i+1,j,0) << endl;
	cout << "HX(i-1,j  ,0) = " << HY(i-1,j,0) << endl;
	cout << "HX(i+1,j-1,0) = " << HX(i+1,j-1,0) << endl;
	cout << "HX(i-1,j-1,0) = " << HY(i-1,j-1,0) << endl;
	cout << "R_P = " << R_P << endl;
	cout << "R_M = " << R_M << endl;
	cout << "==================================================" << endl; 
};

//プリントデバッグ
inline void NsFDTD_TM::printDebugCalcHx(int i, int j, int _i, int _j){
	if(i != _i || j!=_j) return;
	cout << "===========pringDebug CalcHx TM Mode===============" << endl; 
	cout << "(" + to_s(i) + "," + to_s(j) + ")" << endl;
	cout << "HXi,j,+1) = "  << HX(i,j,+1) << endl;
	cout << "HX(i,j, 0) = " << HX(i,j,0) << endl;
	cout << "HX(i,j,-1) = " << HX(i,j,-1) << endl;

	cout << "CHXLY(i,j) = " << CHXLY(i,j) << endl;
	cout << "EZ(i  ,j  , 0) = " << EZ(i,j,0) << endl;
	cout << "EZ(i  ,j+1, 0) = " << EZ(i,j+1,0) << endl;
	cout << "==================================================" << endl; 
}

//プリントデバッグ
inline void NsFDTD_TM::printDebugCalcHy(int i, int j, int _i, int _j){
	if(i != _i || j!=_j) return;
	cout << "===========pringDebug CalcHy TM Mode===============" << endl; 
	cout << "(" + to_s(i) + "," + to_s(j) + ")" << endl;
	cout << "HYi,j,+1) = "  << HY(i,j,+1) << endl;
	cout << "HY(i,j, 0) = " << HY(i,j,0) << endl;
	cout << "HY(i,j,-1) = " << HY(i,j,-1) << endl;

	cout << "CHYLX(i,j) = " << CHYLX(i,j) << endl;
	cout << "EZ(i  ,j  , 0) = " << EZ(i,j,0) << endl;
	cout << "EZ(i+1,j  , 0) = " << EZ(i+1,j,0) << endl;
	cout << "==================================================" << endl; 
}


#endif

/*
				EZ(i,j,+1) = CEZ(i,j)*EZ(i,j, 0)
					+ CEZLH(i,j)*(
					R_P*( HY(i+1,j, 0) - HY(i,j, 0) - ( HX(i,j+1, 0) - HX(i,j, 0) )  )	//dx(Hy) - dy(Hy)
					+R_M*(  Dx2_n(Hy,i,j, 0) - Dy2_n(Hx,i,j, 0) )
					);
					*/
//HX(i,j, +1) = HX(i,j, 0) - CHXLY(i,j)*( EZ(i,j, +1) - EZ(i,j-1, +1) );	//Hxの計算 Hx(i, j+1/2) -> Hx[i,j]
			
//HY(i,j, +1) = HY(i,j, 0) + CHYLX(i,j)*( EZ(i,j, +1) - EZ(i-1,j, +1) );	//Hyの計算 Hy(i+1/2, j) -> Hy[i,j]