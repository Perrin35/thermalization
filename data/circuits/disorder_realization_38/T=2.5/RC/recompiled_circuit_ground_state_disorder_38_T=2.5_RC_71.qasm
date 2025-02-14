OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.030169686) q[0];
sx q[0];
rz(4.5017894) q[0];
sx q[0];
rz(5.9202249) q[0];
rz(1.9650004) q[1];
sx q[1];
rz(-1.3624374) q[1];
sx q[1];
rz(-1.6972313) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1189277) q[0];
sx q[0];
rz(-1.3937409) q[0];
sx q[0];
rz(-1.4757968) q[0];
rz(1.6955201) q[2];
sx q[2];
rz(-1.2850003) q[2];
sx q[2];
rz(0.054626183) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.23420424) q[1];
sx q[1];
rz(-1.475913) q[1];
sx q[1];
rz(2.0070875) q[1];
x q[2];
rz(-1.5115755) q[3];
sx q[3];
rz(-2.539336) q[3];
sx q[3];
rz(1.8730375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.968367) q[2];
sx q[2];
rz(-2.3543365) q[2];
sx q[2];
rz(-3.0464029) q[2];
rz(1.3683052) q[3];
sx q[3];
rz(-0.94822001) q[3];
sx q[3];
rz(0.29909721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66265166) q[0];
sx q[0];
rz(-2.3421685) q[0];
sx q[0];
rz(-0.69315243) q[0];
rz(0.33723801) q[1];
sx q[1];
rz(-1.8353029) q[1];
sx q[1];
rz(-2.3487924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0017298) q[0];
sx q[0];
rz(-0.81714367) q[0];
sx q[0];
rz(-0.344622) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9117965) q[2];
sx q[2];
rz(-1.0655128) q[2];
sx q[2];
rz(1.2838703) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9321654) q[1];
sx q[1];
rz(-1.8575728) q[1];
sx q[1];
rz(-2.3580103) q[1];
x q[2];
rz(-0.81962193) q[3];
sx q[3];
rz(-1.2030798) q[3];
sx q[3];
rz(-1.468827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1482131) q[2];
sx q[2];
rz(-2.9714163) q[2];
sx q[2];
rz(1.4519838) q[2];
rz(-2.4827237) q[3];
sx q[3];
rz(-1.4631319) q[3];
sx q[3];
rz(-2.7448867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0186998) q[0];
sx q[0];
rz(-2.1227699) q[0];
sx q[0];
rz(-3.0320211) q[0];
rz(-0.84596363) q[1];
sx q[1];
rz(-2.1956317) q[1];
sx q[1];
rz(-0.74839655) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7331241) q[0];
sx q[0];
rz(-1.6604024) q[0];
sx q[0];
rz(-1.0138847) q[0];
x q[1];
rz(1.0399299) q[2];
sx q[2];
rz(-2.0874839) q[2];
sx q[2];
rz(0.45848192) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33967623) q[1];
sx q[1];
rz(-1.5775562) q[1];
sx q[1];
rz(1.5077782) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.844039) q[3];
sx q[3];
rz(-0.58876172) q[3];
sx q[3];
rz(1.8901643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1154068) q[2];
sx q[2];
rz(-1.4226961) q[2];
sx q[2];
rz(2.4288948) q[2];
rz(-1.6648434) q[3];
sx q[3];
rz(-1.289117) q[3];
sx q[3];
rz(3.0570928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0067714) q[0];
sx q[0];
rz(-1.1399784) q[0];
sx q[0];
rz(-0.086932927) q[0];
rz(-0.25761071) q[1];
sx q[1];
rz(-2.0517709) q[1];
sx q[1];
rz(1.2923406) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0213762) q[0];
sx q[0];
rz(-2.3142186) q[0];
sx q[0];
rz(2.8143456) q[0];
x q[1];
rz(0.33165641) q[2];
sx q[2];
rz(-1.7229967) q[2];
sx q[2];
rz(-0.70144049) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7304157) q[1];
sx q[1];
rz(-2.431793) q[1];
sx q[1];
rz(-1.0214772) q[1];
rz(0.35859006) q[3];
sx q[3];
rz(-1.7106283) q[3];
sx q[3];
rz(2.7519325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8452235) q[2];
sx q[2];
rz(-1.951428) q[2];
sx q[2];
rz(0.84687084) q[2];
rz(1.752468) q[3];
sx q[3];
rz(-1.4877157) q[3];
sx q[3];
rz(-0.59051591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95319372) q[0];
sx q[0];
rz(-1.1253091) q[0];
sx q[0];
rz(0.97287792) q[0];
rz(-2.1994622) q[1];
sx q[1];
rz(-0.82754389) q[1];
sx q[1];
rz(0.00076248893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036569422) q[0];
sx q[0];
rz(-0.75398406) q[0];
sx q[0];
rz(-2.3466909) q[0];
rz(-pi) q[1];
rz(-0.63456391) q[2];
sx q[2];
rz(-1.3009614) q[2];
sx q[2];
rz(0.25411221) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3734596) q[1];
sx q[1];
rz(-2.4985857) q[1];
sx q[1];
rz(-2.0483584) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0263264) q[3];
sx q[3];
rz(-0.59820931) q[3];
sx q[3];
rz(2.0245352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3804649) q[2];
sx q[2];
rz(-1.7385812) q[2];
sx q[2];
rz(-0.38499704) q[2];
rz(0.72832406) q[3];
sx q[3];
rz(-2.5326122) q[3];
sx q[3];
rz(-2.716632) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7321135) q[0];
sx q[0];
rz(-0.80687579) q[0];
sx q[0];
rz(-0.4581067) q[0];
rz(1.7181646) q[1];
sx q[1];
rz(-2.3511395) q[1];
sx q[1];
rz(3.0480393) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3781965) q[0];
sx q[0];
rz(-2.3468821) q[0];
sx q[0];
rz(2.2346157) q[0];
rz(0.92624591) q[2];
sx q[2];
rz(-2.8684596) q[2];
sx q[2];
rz(-0.63970837) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1653891) q[1];
sx q[1];
rz(-1.7862583) q[1];
sx q[1];
rz(-0.54241753) q[1];
x q[2];
rz(-0.5597975) q[3];
sx q[3];
rz(-1.5889394) q[3];
sx q[3];
rz(-1.960235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88390049) q[2];
sx q[2];
rz(-2.7187686) q[2];
sx q[2];
rz(-0.60919961) q[2];
rz(-0.60112634) q[3];
sx q[3];
rz(-1.6060035) q[3];
sx q[3];
rz(-2.6570184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5179829) q[0];
sx q[0];
rz(-1.4731151) q[0];
sx q[0];
rz(1.3129039) q[0];
rz(0.081722109) q[1];
sx q[1];
rz(-1.3879958) q[1];
sx q[1];
rz(2.0770226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3400187) q[0];
sx q[0];
rz(-0.94807762) q[0];
sx q[0];
rz(-2.4452371) q[0];
rz(-2.1871952) q[2];
sx q[2];
rz(-1.5845044) q[2];
sx q[2];
rz(0.15132667) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.21614046) q[1];
sx q[1];
rz(-1.6419159) q[1];
sx q[1];
rz(2.1468162) q[1];
x q[2];
rz(-0.96780583) q[3];
sx q[3];
rz(-1.9120312) q[3];
sx q[3];
rz(1.6399217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.99819034) q[2];
sx q[2];
rz(-0.4021796) q[2];
sx q[2];
rz(-1.1581988) q[2];
rz(0.68425933) q[3];
sx q[3];
rz(-1.3254157) q[3];
sx q[3];
rz(-1.5023331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96712464) q[0];
sx q[0];
rz(-3.0376349) q[0];
sx q[0];
rz(-0.52484584) q[0];
rz(1.9606494) q[1];
sx q[1];
rz(-0.8371822) q[1];
sx q[1];
rz(2.4139074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0200266) q[0];
sx q[0];
rz(-1.354573) q[0];
sx q[0];
rz(-0.89631594) q[0];
x q[1];
rz(-2.0725756) q[2];
sx q[2];
rz(-2.2282218) q[2];
sx q[2];
rz(-1.9819966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39036738) q[1];
sx q[1];
rz(-1.8341244) q[1];
sx q[1];
rz(0.45231426) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7102881) q[3];
sx q[3];
rz(-0.65912853) q[3];
sx q[3];
rz(-2.0122676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9119447) q[2];
sx q[2];
rz(-1.0934528) q[2];
sx q[2];
rz(1.7952807) q[2];
rz(2.4800269) q[3];
sx q[3];
rz(-1.1450359) q[3];
sx q[3];
rz(0.0865817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10053703) q[0];
sx q[0];
rz(-2.4424398) q[0];
sx q[0];
rz(0.93851411) q[0];
rz(1.8427294) q[1];
sx q[1];
rz(-2.2524736) q[1];
sx q[1];
rz(-0.13359698) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.035048) q[0];
sx q[0];
rz(-1.4770916) q[0];
sx q[0];
rz(-2.8117715) q[0];
rz(-pi) q[1];
rz(0.57564013) q[2];
sx q[2];
rz(-1.8274022) q[2];
sx q[2];
rz(-0.20028608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5028959) q[1];
sx q[1];
rz(-1.2112482) q[1];
sx q[1];
rz(-1.4877968) q[1];
rz(-0.64188285) q[3];
sx q[3];
rz(-1.8699904) q[3];
sx q[3];
rz(-3.0282216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9962697) q[2];
sx q[2];
rz(-1.0066373) q[2];
sx q[2];
rz(2.2553196) q[2];
rz(0.26850548) q[3];
sx q[3];
rz(-1.0849181) q[3];
sx q[3];
rz(1.6296384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4292384) q[0];
sx q[0];
rz(-2.7463284) q[0];
sx q[0];
rz(-2.9048753) q[0];
rz(-0.16353823) q[1];
sx q[1];
rz(-1.5312559) q[1];
sx q[1];
rz(-0.18855655) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8343413) q[0];
sx q[0];
rz(-2.3545958) q[0];
sx q[0];
rz(-2.045162) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8256906) q[2];
sx q[2];
rz(-1.1045278) q[2];
sx q[2];
rz(2.3026274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1592387) q[1];
sx q[1];
rz(-0.51430632) q[1];
sx q[1];
rz(-1.4532754) q[1];
x q[2];
rz(-1.7018407) q[3];
sx q[3];
rz(-1.1342955) q[3];
sx q[3];
rz(2.1646162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4208372) q[2];
sx q[2];
rz(-1.6381702) q[2];
sx q[2];
rz(-2.3756964) q[2];
rz(0.01865538) q[3];
sx q[3];
rz(-1.3466287) q[3];
sx q[3];
rz(0.53001058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6388549) q[0];
sx q[0];
rz(-1.1549594) q[0];
sx q[0];
rz(1.4422944) q[0];
rz(0.28932183) q[1];
sx q[1];
rz(-2.0850291) q[1];
sx q[1];
rz(-2.6883968) q[1];
rz(-1.7062023) q[2];
sx q[2];
rz(-2.1056021) q[2];
sx q[2];
rz(3.0816002) q[2];
rz(2.0404242) q[3];
sx q[3];
rz(-1.508699) q[3];
sx q[3];
rz(-1.1796622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
