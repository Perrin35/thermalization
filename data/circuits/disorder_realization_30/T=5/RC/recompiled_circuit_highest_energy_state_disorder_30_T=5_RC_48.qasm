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
rz(-2.8415866) q[0];
sx q[0];
rz(-2.0553148) q[0];
sx q[0];
rz(-2.3334184) q[0];
rz(-1.1440682) q[1];
sx q[1];
rz(-0.77163982) q[1];
sx q[1];
rz(1.5626524) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90955665) q[0];
sx q[0];
rz(-2.252451) q[0];
sx q[0];
rz(-2.2907592) q[0];
rz(-pi) q[1];
rz(-1.3035266) q[2];
sx q[2];
rz(-1.7198945) q[2];
sx q[2];
rz(1.9097415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11400819) q[1];
sx q[1];
rz(-1.6817087) q[1];
sx q[1];
rz(-1.0618889) q[1];
rz(-pi) q[2];
rz(3.0959356) q[3];
sx q[3];
rz(-1.6049269) q[3];
sx q[3];
rz(-1.0505763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.943104) q[2];
sx q[2];
rz(-2.9297332) q[2];
sx q[2];
rz(-1.0342106) q[2];
rz(-0.69283038) q[3];
sx q[3];
rz(-1.0537909) q[3];
sx q[3];
rz(-1.0606162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8387872) q[0];
sx q[0];
rz(-0.028258709) q[0];
sx q[0];
rz(0.57648188) q[0];
rz(-0.020596404) q[1];
sx q[1];
rz(-0.3978022) q[1];
sx q[1];
rz(2.1131262) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3555455) q[0];
sx q[0];
rz(-1.3554327) q[0];
sx q[0];
rz(0.26305612) q[0];
rz(-pi) q[1];
rz(2.2510301) q[2];
sx q[2];
rz(-0.46698585) q[2];
sx q[2];
rz(-0.61636965) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0513269) q[1];
sx q[1];
rz(-1.169186) q[1];
sx q[1];
rz(0.36220596) q[1];
rz(-0.64739689) q[3];
sx q[3];
rz(-1.2522373) q[3];
sx q[3];
rz(-0.35435521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2241406) q[2];
sx q[2];
rz(-2.1805306) q[2];
sx q[2];
rz(1.2646328) q[2];
rz(0.97359109) q[3];
sx q[3];
rz(-1.4507989) q[3];
sx q[3];
rz(0.26315954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6384386) q[0];
sx q[0];
rz(-1.3042903) q[0];
sx q[0];
rz(-2.2880182) q[0];
rz(1.0517993) q[1];
sx q[1];
rz(-2.167326) q[1];
sx q[1];
rz(-3.1265756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1117843) q[0];
sx q[0];
rz(-1.5620323) q[0];
sx q[0];
rz(0.86017227) q[0];
rz(-pi) q[1];
rz(-0.096108561) q[2];
sx q[2];
rz(-1.5920361) q[2];
sx q[2];
rz(0.48890314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89755631) q[1];
sx q[1];
rz(-1.4970395) q[1];
sx q[1];
rz(2.4270127) q[1];
rz(-1.5068568) q[3];
sx q[3];
rz(-1.2815164) q[3];
sx q[3];
rz(2.6562986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0217648) q[2];
sx q[2];
rz(-1.652635) q[2];
sx q[2];
rz(-0.28142288) q[2];
rz(1.6875387) q[3];
sx q[3];
rz(-1.931793) q[3];
sx q[3];
rz(2.37229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5467095) q[0];
sx q[0];
rz(-1.2458206) q[0];
sx q[0];
rz(2.967714) q[0];
rz(-1.8100544) q[1];
sx q[1];
rz(-1.7200108) q[1];
sx q[1];
rz(-1.4283659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3417018) q[0];
sx q[0];
rz(-0.19289324) q[0];
sx q[0];
rz(1.4961722) q[0];
rz(-2.8160353) q[2];
sx q[2];
rz(-2.7635305) q[2];
sx q[2];
rz(-0.69272536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.169379) q[1];
sx q[1];
rz(-1.6092818) q[1];
sx q[1];
rz(0.65757081) q[1];
x q[2];
rz(-2.3263116) q[3];
sx q[3];
rz(-2.6297556) q[3];
sx q[3];
rz(-2.3533604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7017158) q[2];
sx q[2];
rz(-0.82650799) q[2];
sx q[2];
rz(0.41713777) q[2];
rz(-0.16962984) q[3];
sx q[3];
rz(-0.093234213) q[3];
sx q[3];
rz(0.72181845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193264) q[0];
sx q[0];
rz(-2.7665783) q[0];
sx q[0];
rz(-0.30935031) q[0];
rz(-2.3720692) q[1];
sx q[1];
rz(-2.0517495) q[1];
sx q[1];
rz(1.5187029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5734539) q[0];
sx q[0];
rz(-0.72069695) q[0];
sx q[0];
rz(1.9573523) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3974472) q[2];
sx q[2];
rz(-2.2501906) q[2];
sx q[2];
rz(-1.3199922) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57598439) q[1];
sx q[1];
rz(-0.19153015) q[1];
sx q[1];
rz(-1.7267141) q[1];
rz(-pi) q[2];
rz(-1.7047911) q[3];
sx q[3];
rz(-1.4729285) q[3];
sx q[3];
rz(2.0201473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7421444) q[2];
sx q[2];
rz(-0.98230201) q[2];
sx q[2];
rz(2.1461416) q[2];
rz(-0.71850145) q[3];
sx q[3];
rz(-1.813443) q[3];
sx q[3];
rz(0.79513454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28498483) q[0];
sx q[0];
rz(-1.4549078) q[0];
sx q[0];
rz(1.6108151) q[0];
rz(1.6948304) q[1];
sx q[1];
rz(-1.6981354) q[1];
sx q[1];
rz(-1.4379427) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73285028) q[0];
sx q[0];
rz(-2.7332531) q[0];
sx q[0];
rz(-1.6224003) q[0];
rz(-pi) q[1];
rz(-1.1968568) q[2];
sx q[2];
rz(-2.1293921) q[2];
sx q[2];
rz(0.93511673) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29026689) q[1];
sx q[1];
rz(-1.0630634) q[1];
sx q[1];
rz(1.6666056) q[1];
rz(3.1385057) q[3];
sx q[3];
rz(-0.75957889) q[3];
sx q[3];
rz(0.39955968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97856727) q[2];
sx q[2];
rz(-1.7820396) q[2];
sx q[2];
rz(-1.9095518) q[2];
rz(-1.1466522) q[3];
sx q[3];
rz(-1.8012643) q[3];
sx q[3];
rz(0.31057096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50800407) q[0];
sx q[0];
rz(-0.81729832) q[0];
sx q[0];
rz(2.001413) q[0];
rz(2.2445402) q[1];
sx q[1];
rz(-1.3373809) q[1];
sx q[1];
rz(0.030489771) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59647467) q[0];
sx q[0];
rz(-1.6870572) q[0];
sx q[0];
rz(1.3482987) q[0];
rz(-pi) q[1];
rz(1.8980128) q[2];
sx q[2];
rz(-0.95409617) q[2];
sx q[2];
rz(-1.1730574) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5329689) q[1];
sx q[1];
rz(-2.2734959) q[1];
sx q[1];
rz(-2.7547902) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6875911) q[3];
sx q[3];
rz(-2.2488065) q[3];
sx q[3];
rz(-2.0583722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9474779) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(-1.4878368) q[3];
sx q[3];
rz(-1.0349118) q[3];
sx q[3];
rz(-1.3207159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79539913) q[0];
sx q[0];
rz(-0.29618707) q[0];
sx q[0];
rz(-1.7561703) q[0];
rz(-1.9253383) q[1];
sx q[1];
rz(-1.6938035) q[1];
sx q[1];
rz(2.6522955) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74936825) q[0];
sx q[0];
rz(-2.4697692) q[0];
sx q[0];
rz(-3.1098614) q[0];
rz(1.9070315) q[2];
sx q[2];
rz(-1.6052941) q[2];
sx q[2];
rz(-2.2740342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9605105) q[1];
sx q[1];
rz(-2.0415293) q[1];
sx q[1];
rz(3.0357643) q[1];
x q[2];
rz(-1.9527444) q[3];
sx q[3];
rz(-2.8738465) q[3];
sx q[3];
rz(1.4827888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3936002) q[2];
sx q[2];
rz(-1.1425428) q[2];
sx q[2];
rz(-2.6753814) q[2];
rz(1.6097869) q[3];
sx q[3];
rz(-1.6377056) q[3];
sx q[3];
rz(-2.8209749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521249) q[0];
sx q[0];
rz(-0.89972275) q[0];
sx q[0];
rz(3.0925282) q[0];
rz(-2.9009254) q[1];
sx q[1];
rz(-1.0180232) q[1];
sx q[1];
rz(0.022445591) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81240679) q[0];
sx q[0];
rz(-2.5617449) q[0];
sx q[0];
rz(1.5270698) q[0];
rz(-pi) q[1];
rz(0.38752611) q[2];
sx q[2];
rz(-1.8948613) q[2];
sx q[2];
rz(1.9112196) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7514403) q[1];
sx q[1];
rz(-0.94822394) q[1];
sx q[1];
rz(2.1639813) q[1];
rz(-pi) q[2];
rz(2.4565036) q[3];
sx q[3];
rz(-1.9307319) q[3];
sx q[3];
rz(1.0459096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13751328) q[2];
sx q[2];
rz(-1.235032) q[2];
sx q[2];
rz(2.1732886) q[2];
rz(2.3051895) q[3];
sx q[3];
rz(-0.29763779) q[3];
sx q[3];
rz(-1.77232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4150998) q[0];
sx q[0];
rz(-1.8510171) q[0];
sx q[0];
rz(0.37877628) q[0];
rz(-0.0060161034) q[1];
sx q[1];
rz(-2.6117987) q[1];
sx q[1];
rz(2.5078497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0519052) q[0];
sx q[0];
rz(-1.2335334) q[0];
sx q[0];
rz(0.71292065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1442398) q[2];
sx q[2];
rz(-1.4359983) q[2];
sx q[2];
rz(0.9465132) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4162916) q[1];
sx q[1];
rz(-1.0868727) q[1];
sx q[1];
rz(-2.121622) q[1];
rz(2.2489088) q[3];
sx q[3];
rz(-1.0859981) q[3];
sx q[3];
rz(-0.78692737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4358431) q[2];
sx q[2];
rz(-1.6195932) q[2];
sx q[2];
rz(-0.61326927) q[2];
rz(0.67522007) q[3];
sx q[3];
rz(-2.7628511) q[3];
sx q[3];
rz(-0.81775445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9153862) q[0];
sx q[0];
rz(-1.7788667) q[0];
sx q[0];
rz(-1.9052196) q[0];
rz(-0.18648237) q[1];
sx q[1];
rz(-1.7541371) q[1];
sx q[1];
rz(-1.6288155) q[1];
rz(-1.5908505) q[2];
sx q[2];
rz(-2.0303844) q[2];
sx q[2];
rz(-2.5302237) q[2];
rz(0.091681176) q[3];
sx q[3];
rz(-2.1909919) q[3];
sx q[3];
rz(1.7237678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
