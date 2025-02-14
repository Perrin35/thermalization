OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.74855411) q[0];
sx q[0];
rz(-1.3286123) q[0];
sx q[0];
rz(-2.7207029) q[0];
rz(2.4731877) q[1];
sx q[1];
rz(-1.547812) q[1];
sx q[1];
rz(-1.9903543) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52510696) q[0];
sx q[0];
rz(-0.49607222) q[0];
sx q[0];
rz(0.43323364) q[0];
x q[1];
rz(1.0082863) q[2];
sx q[2];
rz(-1.9676625) q[2];
sx q[2];
rz(-2.4728554) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9734263) q[1];
sx q[1];
rz(-0.81422537) q[1];
sx q[1];
rz(2.3553215) q[1];
x q[2];
rz(0.26845308) q[3];
sx q[3];
rz(-0.9359979) q[3];
sx q[3];
rz(0.47273794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0737754) q[2];
sx q[2];
rz(-2.1442118) q[2];
sx q[2];
rz(-1.0214405) q[2];
rz(1.3218309) q[3];
sx q[3];
rz(-0.50522155) q[3];
sx q[3];
rz(-2.5652313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3385056) q[0];
sx q[0];
rz(-1.0271238) q[0];
sx q[0];
rz(-0.3013674) q[0];
rz(1.1400247) q[1];
sx q[1];
rz(-0.81914425) q[1];
sx q[1];
rz(1.0637306) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6607644) q[0];
sx q[0];
rz(-1.4505634) q[0];
sx q[0];
rz(1.3025492) q[0];
x q[1];
rz(2.930141) q[2];
sx q[2];
rz(-0.7728979) q[2];
sx q[2];
rz(-2.7711902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96460198) q[1];
sx q[1];
rz(-0.97787217) q[1];
sx q[1];
rz(-1.9511392) q[1];
rz(1.0349991) q[3];
sx q[3];
rz(-1.4685653) q[3];
sx q[3];
rz(-0.4680948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86338824) q[2];
sx q[2];
rz(-2.3953343) q[2];
sx q[2];
rz(-2.9928652) q[2];
rz(-1.1778098) q[3];
sx q[3];
rz(-1.1921459) q[3];
sx q[3];
rz(-1.8671794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3026368) q[0];
sx q[0];
rz(-1.9479072) q[0];
sx q[0];
rz(1.0021915) q[0];
rz(-1.3305371) q[1];
sx q[1];
rz(-1.1794773) q[1];
sx q[1];
rz(2.5340714) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1992545) q[0];
sx q[0];
rz(-2.1954241) q[0];
sx q[0];
rz(2.6074431) q[0];
rz(0.98812466) q[2];
sx q[2];
rz(-2.6651504) q[2];
sx q[2];
rz(0.63240766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.62251112) q[1];
sx q[1];
rz(-1.4591328) q[1];
sx q[1];
rz(-0.51736781) q[1];
x q[2];
rz(-2.420531) q[3];
sx q[3];
rz(-2.0264056) q[3];
sx q[3];
rz(-0.30649116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9017631) q[2];
sx q[2];
rz(-2.2679057) q[2];
sx q[2];
rz(2.152781) q[2];
rz(1.4766988) q[3];
sx q[3];
rz(-1.3380545) q[3];
sx q[3];
rz(0.57507676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3636417) q[0];
sx q[0];
rz(-2.2626484) q[0];
sx q[0];
rz(-1.9212035) q[0];
rz(1.8266504) q[1];
sx q[1];
rz(-1.5053791) q[1];
sx q[1];
rz(0.13994089) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5410239) q[0];
sx q[0];
rz(-1.9697088) q[0];
sx q[0];
rz(-2.5068796) q[0];
rz(-pi) q[1];
x q[1];
rz(2.432517) q[2];
sx q[2];
rz(-1.1458414) q[2];
sx q[2];
rz(0.72689012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32711682) q[1];
sx q[1];
rz(-0.45683858) q[1];
sx q[1];
rz(-1.4903699) q[1];
rz(-pi) q[2];
rz(1.1871376) q[3];
sx q[3];
rz(-1.5042437) q[3];
sx q[3];
rz(-2.6714747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6949029) q[2];
sx q[2];
rz(-1.0580772) q[2];
sx q[2];
rz(-3.0042082) q[2];
rz(1.8047699) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(-3.1409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3985652) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(0.16125691) q[0];
rz(-2.297961) q[1];
sx q[1];
rz(-0.7470986) q[1];
sx q[1];
rz(-1.8341433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6661412) q[0];
sx q[0];
rz(-1.8384116) q[0];
sx q[0];
rz(2.3580736) q[0];
rz(-1.0867811) q[2];
sx q[2];
rz(-2.7325411) q[2];
sx q[2];
rz(-1.3889927) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54164179) q[1];
sx q[1];
rz(-2.5516918) q[1];
sx q[1];
rz(-1.7322503) q[1];
rz(-2.3639115) q[3];
sx q[3];
rz(-0.84319226) q[3];
sx q[3];
rz(0.2524337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1696986) q[2];
sx q[2];
rz(-2.7104388) q[2];
sx q[2];
rz(-0.92740721) q[2];
rz(1.6124604) q[3];
sx q[3];
rz(-1.1011139) q[3];
sx q[3];
rz(0.35759887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2138432) q[0];
sx q[0];
rz(-2.1723211) q[0];
sx q[0];
rz(-1.7359605) q[0];
rz(-1.7270145) q[1];
sx q[1];
rz(-2.3443293) q[1];
sx q[1];
rz(2.3419211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9093577) q[0];
sx q[0];
rz(-0.063061558) q[0];
sx q[0];
rz(-1.7574278) q[0];
x q[1];
rz(0.39787103) q[2];
sx q[2];
rz(-2.8074773) q[2];
sx q[2];
rz(-2.1159005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0015687167) q[1];
sx q[1];
rz(-2.713361) q[1];
sx q[1];
rz(-2.2111441) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42678435) q[3];
sx q[3];
rz(-2.7023661) q[3];
sx q[3];
rz(-2.3595485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8863525) q[2];
sx q[2];
rz(-2.4031874) q[2];
sx q[2];
rz(-0.8026455) q[2];
rz(2.8211527) q[3];
sx q[3];
rz(-1.3405864) q[3];
sx q[3];
rz(-1.505544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7771512) q[0];
sx q[0];
rz(-0.027712263) q[0];
sx q[0];
rz(0.030315422) q[0];
rz(1.3069356) q[1];
sx q[1];
rz(-1.0825284) q[1];
sx q[1];
rz(2.511715) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8960285) q[0];
sx q[0];
rz(-1.3828074) q[0];
sx q[0];
rz(2.7134368) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0391049) q[2];
sx q[2];
rz(-1.9699735) q[2];
sx q[2];
rz(0.98083996) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67757183) q[1];
sx q[1];
rz(-1.7742762) q[1];
sx q[1];
rz(-2.8427441) q[1];
x q[2];
rz(-1.685018) q[3];
sx q[3];
rz(-2.0124314) q[3];
sx q[3];
rz(0.87678443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.031875413) q[2];
sx q[2];
rz(-0.66706359) q[2];
sx q[2];
rz(-1.1472222) q[2];
rz(2.4540497) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(-1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3654093) q[0];
sx q[0];
rz(-2.9054346) q[0];
sx q[0];
rz(0.48015204) q[0];
rz(-2.7393553) q[1];
sx q[1];
rz(-1.6419342) q[1];
sx q[1];
rz(0.38280907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3881809) q[0];
sx q[0];
rz(-2.7482902) q[0];
sx q[0];
rz(-2.1476782) q[0];
x q[1];
rz(0.11522861) q[2];
sx q[2];
rz(-1.8602716) q[2];
sx q[2];
rz(0.25184271) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3814712) q[1];
sx q[1];
rz(-0.49182004) q[1];
sx q[1];
rz(2.2342938) q[1];
rz(0.14842378) q[3];
sx q[3];
rz(-0.19834861) q[3];
sx q[3];
rz(-2.655011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.39168921) q[2];
sx q[2];
rz(-0.73035556) q[2];
sx q[2];
rz(2.7833617) q[2];
rz(0.41424888) q[3];
sx q[3];
rz(-2.1130989) q[3];
sx q[3];
rz(2.7481368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8527894) q[0];
sx q[0];
rz(-1.8599334) q[0];
sx q[0];
rz(-0.417867) q[0];
rz(-1.4990384) q[1];
sx q[1];
rz(-2.7211029) q[1];
sx q[1];
rz(-0.61161673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8369524) q[0];
sx q[0];
rz(-1.3654183) q[0];
sx q[0];
rz(2.9083328) q[0];
x q[1];
rz(2.8496004) q[2];
sx q[2];
rz(-0.21459178) q[2];
sx q[2];
rz(-2.8923182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7193209) q[1];
sx q[1];
rz(-0.55972717) q[1];
sx q[1];
rz(-2.1156963) q[1];
x q[2];
rz(0.34682746) q[3];
sx q[3];
rz(-1.4937823) q[3];
sx q[3];
rz(1.6517757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.186782) q[2];
sx q[2];
rz(-1.6839226) q[2];
sx q[2];
rz(-0.35624722) q[2];
rz(0.48480836) q[3];
sx q[3];
rz(-1.7907413) q[3];
sx q[3];
rz(0.029732186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.818882) q[0];
sx q[0];
rz(-2.0368545) q[0];
sx q[0];
rz(-1.6792962) q[0];
rz(-2.7541584) q[1];
sx q[1];
rz(-2.2144364) q[1];
sx q[1];
rz(-0.80518728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062267508) q[0];
sx q[0];
rz(-2.1342417) q[0];
sx q[0];
rz(-0.84765537) q[0];
rz(-pi) q[1];
rz(2.7323805) q[2];
sx q[2];
rz(-1.1537038) q[2];
sx q[2];
rz(-0.52601782) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8987052) q[1];
sx q[1];
rz(-1.4243898) q[1];
sx q[1];
rz(-2.8553748) q[1];
rz(-1.8191387) q[3];
sx q[3];
rz(-1.2183918) q[3];
sx q[3];
rz(3.0737511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8055973) q[2];
sx q[2];
rz(-0.62312859) q[2];
sx q[2];
rz(-1.222329) q[2];
rz(2.3313816) q[3];
sx q[3];
rz(-0.92456341) q[3];
sx q[3];
rz(1.5626102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41351086) q[0];
sx q[0];
rz(-1.5039197) q[0];
sx q[0];
rz(1.529827) q[0];
rz(-1.275508) q[1];
sx q[1];
rz(-0.38560148) q[1];
sx q[1];
rz(1.7332981) q[1];
rz(-1.1861943) q[2];
sx q[2];
rz(-1.6976962) q[2];
sx q[2];
rz(-3.1023956) q[2];
rz(-0.44966251) q[3];
sx q[3];
rz(-0.95316124) q[3];
sx q[3];
rz(1.1464473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
