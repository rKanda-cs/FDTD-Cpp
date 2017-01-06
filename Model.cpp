#include"Model.h"
#include"Field.h"

/*---------------------------------------------*/
/*--------------�~��Mie�U��--------------------*/
/*---------------------------------------------*/
FazzyMieModel::FazzyMieModel(Field *f, double _r):
FazzyModel(f),r(_r)
{
	ep = 1.6*1.6*EPSILON_0_S;			//�U�d�� = (���ܗ�)^2
		cout << "r=" + to_s((int)mField->cellToNano(r)) << endl;
}

string FazzyMieModel::mkdir(string root){
	_mkdir((root + "Mie").c_str());

	string name = "Mie/" + to_s((int)(mField->cellToNano(r))) +"nm,"+ mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

double FazzyMieModel::calcEPS(const double& x, const double& y, enum INTEG f){

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)�����_�ɃV�t�g
	double _y = my - 0.5*mField->getNy();

	//���S�Ƃ̋��������a+��2/2�Z���ȏ�Ȃ�, ���S�ɔ}���̊O(�O�̂���, ���a+1 �ȏォ���ׂĂ���)
	if(_x*_x + _y*_y >= pow(r+1, 2.0))

		return EPSILON_0_S;

	//���S�Ƃ̋�����, ���a-��2/2�Z���ȓ��Ȃ�, ���S�ɔ}���̊O
	if(_x*_x + _y*_y <= pow(r-1, 2.0))
		return ep;

	double s=0;

	double a = 1.0,b=1.0;
	if(f == D_X) b = 0;
	if(f == D_Y) a = 0;
	for(double i=-16+0.5; i<16; i+=1)
		for(double j=-16+0.5; j<16; j+=1)
			if(pow(_x+a*i/32.0, 2.0) + pow(_y+b*j/32.0, 2.0) <= r*r)
				s+=1; 
	s /= 32.0*32.0;
	return ep*s + EPSILON_0_S*(1-s);
}

/*---------------------------------------------*/
/*--------------���w��-------------------------*/
/*---------------------------------------------*/
FazzySlabModel::FazzySlabModel(Field* f):
FazzyModel(f), ep1(2.0*2.0*EPSILON_0_S), ep2(EPSILON_0_S), width1(250), width2(50)
{
}

double FazzySlabModel::calcEPS(const double& x, const double& y, enum INTEG f){
//��100nm����,250nm�Ԋu��50nm�̃X���u�����Ă���  **��250nm����(L70.71)10nm�X���u�ɕύX(L73)
//���w��
	
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();

	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	int k    = (int)(mField->cellToNano(mx) - 250)%250;
	double l =      (mField->cellToNano(mx) - 250)/250;

	if( k > 0 && k <=10 && l < 5)
		return ep1;
	else
		return ep2;

}

string FazzySlabModel::mkdir(string root){
	_mkdir((root + "SlabModel").c_str());

	string name = "SlabModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

/*---------------------------------------------*/
/*---------------�є�--------------------------*/
/*---------------------------------------------*/
/* �c�f�� */
FazzyHair_incidenceModel::FazzyHair_incidenceModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), alpha(3), cwidth(1), r(32+(8-cwidth))
	//alpha:�L���[�e�B�N���̊p�x(deg)  length:�L���[�e�B�N���̕�(��m)  r:�т̔��a(��m)(���a+�L���[�e�B�N�����d�Ȃ�̈�)
{
	alphaR = alpha * PI / 180;
	length = cwidth / sin(alphaR);
	cout << "�L���[�e�B�N���̊p�x : " + to_s(alpha) + "deg" << endl;
	cout << "�L���[�e�B�N���� : " + to_s(cwidth) + "micro" << endl;
	cout << "�L���[�e�B�N��1���̘I�o�� : " + to_s(length) + "micro" << endl;
}

double FazzyHair_incidenceModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	rn = mField->nanoToCell(r * 1000);
	
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;
	
	double h = mField->nanoToCell(0*1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�
/*
	if (alpha == 0.0) {		//�L���[�e�B�N���Ȃ��̏ꍇ
		if (my <= rn + cy)	return ep1;
		else  return ep2;
	}
*/
	int c = mField->getNx() / lx + 1;		//�v�Z�͈͓��̃L���[�e�B�N���̐�
	for (int i = 0; i < c; i++) {
		if (mx > i * lx + h && mx < (i + 1) * lx + h && mx < mField->getNx() - h) {
			//			if (my > tan(alphaR) * (mx - lx*i) + cy + rn)	return ep2;
			//			else return ep1;		//Fuzzy�Ȃ�(Staircase���f��)

			double dy1 = my - (tan(alphaR) * (mx - lx*i - h) + cy + rn);
			double dy2 = my - (tan(alphaR) * ((mx - lx*i - h) + 1) + cy + rn);
			double s;
			if (dy1 > 0 && dy2 > 0) return ep2;		//�L���[�e�B�N�������̊O�� (1)
			if (fabs(dy1) > 1 && fabs(dy2) > 1) return ep1;		//�L���[�e�B�N�������̓��� (2)

			if (dy1 <= 0 && dy2 <= 0) {
				if (fabs(dy1) <= 1 && fabs(dy2) <= 1) {
					s = (fabs(dy1) + fabs(dy2)) * 1.0 / 2.0;
					return ep1 * s + ep2 * (1 - s);		// (3)
				}
				if (fabs(dy1) < 1 && fabs(dy2) > 1) {
					s = (1 - fabs(dy1)) * ((my - cy - rn) / tan(alphaR) - (mx - lx*i - h)) / 2;
					return ep2 * s + ep1 * (1 - s);		// (4)
				}
			}
			if (dy1 > 0 && dy2 < 0) {
				s = fabs(dy2) * (((mx - lx*i - h) + 1) - (my - cy - rn) / tan(alphaR)) / 2;
				return ep1 * s + ep2 * (1 - s);		// (5)
			}
		}
		else
			continue;
			//break;
	}
	return ep2;
}

double FazzyHair_incidenceModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//���F�̐F�f
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_incidenceModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
/*
	_mkdir((root + "HairModel/incidenceplane").c_str());				//�z���W���Ȃ��̏ꍇ
	string name = "HairModel/incidenceplane/" + mField->getStringCellInfo();
*/
	_mkdir((root + "HairModel/incidenceplane_withSig").c_str());		//�z���W������̏ꍇ
	string name = "HairModel/incidenceplane_withSig/" + mField->getStringCellInfo();
	
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	
	return name + "/";
}


/* �c�f��(�L���[�e�B�N���Ȃ�) */
FazzyHair_NONcuticleModel::FazzyHair_NONcuticleModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), r(32)
	//r:�т̔��a(��m)
{
	cout << "�L���[�e�B�N�� : �Ȃ�"<< endl;
}

double FazzyHair_NONcuticleModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�

	if (my <= rn + cy)	return ep1;
	else  return ep2;
}

double FazzyHair_NONcuticleModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//���F�̐F�f
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_NONcuticleModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
/*
	_mkdir((root + "HairModel/NONcuticle").c_str());				//�z���W���Ȃ��̏ꍇ
	string name = "HairModel/NONcuticle/" + mField->getStringCellInfo();
*/
	_mkdir((root + "HairModel/NONcuticle_withSig").c_str());		//�z���W������̏ꍇ
	string name = "HairModel/NONcuticle_withSig/" + mField->getStringCellInfo();

	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}


/* ���f�� */
FazzyHair_normalModel::FazzyHair_normalModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), e(0.6), r(32)
	//a:���S��  r:�т̔��a(��m)
{
	cout << "�ȉ~�̗��S�� = " + to_s((double)e) << endl;
}

double FazzyHair_normalModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);
	ax = rn;
	by = ax * sqrt(1 - e*e);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML�w

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)�����_�ɃV�t�g
	double _y = my - 0.5*mField->getNy();

	//���S�Ƃ̋��������a+��2/2�Z���ȏ�Ȃ�, ���S�ɔ}���̊O(�O�̂���, ���a+1 �ȏォ���ׂĂ���)
	double _ax = ax+1, _by = by+1;
	if ((_x*_x) / (_ax*_ax) + (_y*_y) / (_by*_by) >= 1)
		return ep2;

	_ax = ax - 1;
	_by = by - 1;
	//���S�Ƃ̋�����, ���a-��2/2�Z���ȓ��Ȃ�, ���S�ɔ}���̊O
	if ((_x*_x) / (_ax*_ax) + (_y*_y) / (_by*_by) <= 1)
		return ep1;

	double s = 0;

	double a = 1.0, b = 1.0;
	if (f == D_X) b = 0;
	if (f == D_Y) a = 0;
	for (double i = -16 + 0.5; i<16; i += 1)
		for (double j = -16 + 0.5; j<16; j += 1)
			if (pow(_x + a*i / 32.0, 2.0) / (ax*ax) + pow(_y + b*j / 32.0, 2.0) / (by*by) <= 1)
				s += 1;
	s /= 32.0*32.0;
	return ep1*s + ep2*(1 - s);
}

string FazzyHair_normalModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	_mkdir((root + "HairModel/normalplane").c_str());
	
	string name = "HairModel/normalplane/e=" + to_s((double)e);
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	
	name = "HairModel/normalplane/e=" + to_s((double)e) + "/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

/*---------------------------------------------*/
/*--------------�����t�H��--------------------*/
/*---------------------------------------------*/
FazzyMorphoModel::FazzyMorphoModel(Field* f, double _h0, double _h1, enum STRUCTURE kind):
FazzyModel(f), shelf(kind)
{
	num = 8;				//�ςݏd�˂鐔
	// ep[1] = 3.5*3.5*EPSILON_0_S;	//�U�d��
	// ep[0] = 1.45*1.45*EPSILON_0_S;
	ep[1] = 1.56*1.56*EPSILON_0_S;	//�U�d��
	ep[0] = 1.0*1.0*EPSILON_0_S;
	width = mField->nanoToCell(300);	//������150�ŌŒ�  **�_���ł�300

	min = mField->nanoToCell(120);
	max = mField->nanoToCell(120);
	height[1] = min;
	height[0] = min;
	cout << min << endl;
	cout << max << endl;
}

double FazzyMorphoModel::calcEPS(const double &x, const double &y, enum INTEG f){
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	//N_X/2�𒆐S����,�����`�����Ɍ݂��Ⴂ�ɔz�u
	int dis = height[0] + height[1];
	int oy = (mField->getNy() - num*dis)/2.0;
	int ox = mField->getNx()/2.0;
	double _x = mx-ox;	//ox,oy�����W�̌��_��
	double _y = my-oy;
	
	//���f���̍��E���㉺�ɗ���Ă���Δ}���̊O
	if(abs(_x) > width+1 || abs(_y-1.0*num*dis/2.0) > num*dis/2.0+1)
		return EPSILON_0_S;

	double s[2]={0,0};
	double a = 1.0,b=1.0;
	if(f == D_X) b = 0;
	if(f == D_Y) a = 0;
	for(double i=-16+0.5; i<16; i+=1){
		for(double j=-16+0.5; j<16; j+=1){
			double sx = _x + a*i/32.0;
			double sy = _y + b*j/32.0;
			if( abs(sx) > width  || abs(sy-1.0*num*dis/2.0) > 1.0*num*dis/2.0) continue;

			bool k = ((int)sy%dis+0.5) > height[0];	//���E��Ŕ�ׂȂ��悤��0.5�����Ă���(0���傫��1�����Ȃ牽�ł�����)
			//bool k =  (floor(sy/ dis)*dis < sy) && ( sy < floor(sy/ dis)*dis + height[0]);

			if (sx < 0 && shelf)
				k = !k;		//���E�Ŕ��], �݂��Ⴂ�łȂ������甽�]���Ȃ�

			s[k] +=1;
		}
	}
	s[0] /= 32.0*32.0;
	s[1] /= 32.0*32.0;
	return EPSILON_0_S*(1-s[0]-s[1]) + ep[0]*s[0] + ep[1]*s[1];
}

string FazzyMorphoModel::mkdir(string root){
	string label = "Morpho(" + to_s(sqrt(ep[0])) + "," + to_s(sqrt(ep[1])) + ")M=" + to_s(num);
	//  string label = "Morpho";
	_mkdir((root + label).c_str());
	string name;

	if(shelf)
		name = label + "/" + to_s((int)(mField->cellToNano(height[0]))) + "nm" + mField->getStringCellInfo() ;
	else
		name = label + "/" + to_s((int)(mField->cellToNano(height[0]))) + "nm(nonShelf)" + mField->getStringCellInfo();

	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}


/*--------------------------------*/
/*-----------���f���Ȃ�-----------*/
/*--------------------------------*/
bool FazzyMorphoModel::update(int dh){
	height[0] += (int) mField->nanoToCell(dh);
	height[1] += (int) mField->nanoToCell(dh);

	if(height[0] > max)
		return false;

	return true;
}

string FazzyNoModel::mkdir(string root){
	_mkdir((root + "NoModel").c_str());

	string name = "NoModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}