#ifndef _MODEL_H
#define _MODEL_H
#include"PhysicalConstant.h"
#include<string>
#include<direct.h>
#include"Object.h"
#include"Field.h"
using namespace std;
enum INTEG{
	D_XY = 0,
	D_X,
	D_Y
};

class FazzyModel{
protected:
	Field *mField;
public:
	FazzyModel(Field *field){
		mField = field;
	};
	virtual string mkdir(string root) = 0;
	virtual double calcEPS(const double&, const double&, enum INTEG = D_XY) = 0;
	virtual double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY) = 0;
	virtual bool update(int) = 0;
	virtual void Initialize()=0;
};

class FazzySlabModel :public FazzyModel{
	const double ep1, ep2;
	const int width1, width2;
public:
	FazzySlabModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	bool update(int){
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyMieModel :public FazzyModel{
	double ep;
	double r;
public:
	FazzyMieModel(Field* f, double _r);
	string mkdir(string root);	//ディレクトリ作成

	//誘電率計算
	double calcEPS(const double& x, const double& y, enum INTEG f);
	bool update(int a){
		  return true;
	}
	void Initialize()
	{
	}
};

class FazzyHair_incidenceModel :public FazzyModel {
	const double ep1, ep2;
	const int alpha, cwidth, r;
	double alphaR, length, ln, lx, ly, rn;
public:
	FazzyHair_incidenceModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_normalModel :public FazzyModel {
	const double ep1, ep2, e;
	const int r;
	double rn, ax, by;
public:
	FazzyHair_normalModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG f);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_NONcuticleModel :public FazzyModel {
	const double ep1, ep2;
	const int r;
	double rn;
public:
	FazzyHair_NONcuticleModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

enum STRUCTURE{
	SHELF = true,
	NONSHELF = false
};

class FazzyMorphoModel:public FazzyModel{
protected:
	int height[2];
	int min,max;
	int width;
	double e;
	int num;
	double ep[2];
	enum STRUCTURE shelf;	//互い違いかどうか
public:
	FazzyMorphoModel(Field*, double _h0, double _h1, enum STRUCTURE kind);

	string mkdir(string root);

	double calcEPS(const double &x, const double &y, enum INTEG f);

	bool update(int dh);

	void Initialize(){
		height[1] = min;
		height[0] = min;
	}
};

class FazzyNoModel:public FazzyModel{
public:
	FazzyNoModel(Field* f)
		:FazzyModel(f)
	{
	}

	  string mkdir(string root);

	   double calcEPS(const double &x, const double &y, enum INTEG f){
		   return EPSILON_0_S;
	   }

	   bool update(int a){
		   return true;
	   }

	   void Initialize()
	   {
		  
	   }

};

#endif //_MODEL_H