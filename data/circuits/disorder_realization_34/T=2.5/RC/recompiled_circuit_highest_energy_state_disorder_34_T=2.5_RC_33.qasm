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
rz(-2.3930385) q[0];
sx q[0];
rz(-1.8129803) q[0];
sx q[0];
rz(-0.42088977) q[0];
rz(2.4731877) q[1];
sx q[1];
rz(-1.547812) q[1];
sx q[1];
rz(-1.9903543) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4822455) q[0];
sx q[0];
rz(-1.3696241) q[0];
sx q[0];
rz(2.6849999) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45999476) q[2];
sx q[2];
rz(-2.0850217) q[2];
sx q[2];
rz(-2.0005039) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1382153) q[1];
sx q[1];
rz(-2.1103519) q[1];
sx q[1];
rz(2.2141377) q[1];
x q[2];
rz(2.8731396) q[3];
sx q[3];
rz(-2.2055948) q[3];
sx q[3];
rz(-2.6688547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0737754) q[2];
sx q[2];
rz(-2.1442118) q[2];
sx q[2];
rz(1.0214405) q[2];
rz(-1.8197618) q[3];
sx q[3];
rz(-0.50522155) q[3];
sx q[3];
rz(-2.5652313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8030871) q[0];
sx q[0];
rz(-2.1144688) q[0];
sx q[0];
rz(0.3013674) q[0];
rz(-2.001568) q[1];
sx q[1];
rz(-0.81914425) q[1];
sx q[1];
rz(1.0637306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50142215) q[0];
sx q[0];
rz(-0.29336624) q[0];
sx q[0];
rz(1.9984931) q[0];
x q[1];
rz(1.3688886) q[2];
sx q[2];
rz(-2.3221952) q[2];
sx q[2];
rz(-2.47987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82599431) q[1];
sx q[1];
rz(-1.2578405) q[1];
sx q[1];
rz(0.62773825) q[1];
rz(2.1065936) q[3];
sx q[3];
rz(-1.4685653) q[3];
sx q[3];
rz(-2.6734979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2782044) q[2];
sx q[2];
rz(-2.3953343) q[2];
sx q[2];
rz(-0.14872742) q[2];
rz(-1.9637828) q[3];
sx q[3];
rz(-1.1921459) q[3];
sx q[3];
rz(-1.2744133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3026368) q[0];
sx q[0];
rz(-1.9479072) q[0];
sx q[0];
rz(-1.0021915) q[0];
rz(-1.3305371) q[1];
sx q[1];
rz(-1.1794773) q[1];
sx q[1];
rz(-0.60752121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1992545) q[0];
sx q[0];
rz(-2.1954241) q[0];
sx q[0];
rz(0.53414957) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9776865) q[2];
sx q[2];
rz(-1.8259138) q[2];
sx q[2];
rz(1.4681548) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5190815) q[1];
sx q[1];
rz(-1.4591328) q[1];
sx q[1];
rz(-0.51736781) q[1];
rz(-pi) q[2];
rz(2.420531) q[3];
sx q[3];
rz(-1.115187) q[3];
sx q[3];
rz(-0.30649116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2398296) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3636417) q[0];
sx q[0];
rz(-2.2626484) q[0];
sx q[0];
rz(-1.2203891) q[0];
rz(1.8266504) q[1];
sx q[1];
rz(-1.5053791) q[1];
sx q[1];
rz(-3.0016518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6005687) q[0];
sx q[0];
rz(-1.9697088) q[0];
sx q[0];
rz(2.5068796) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5342388) q[2];
sx q[2];
rz(-0.80728856) q[2];
sx q[2];
rz(-1.2918778) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23754696) q[1];
sx q[1];
rz(-1.1155459) q[1];
sx q[1];
rz(-0.039467875) q[1];
rz(-0.071752944) q[3];
sx q[3];
rz(-1.9535616) q[3];
sx q[3];
rz(2.0677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6949029) q[2];
sx q[2];
rz(-1.0580772) q[2];
sx q[2];
rz(-3.0042082) q[2];
rz(1.8047699) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(0.00060753879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7430275) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(2.9803357) q[0];
rz(-0.84363168) q[1];
sx q[1];
rz(-2.3944941) q[1];
sx q[1];
rz(-1.8341433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8358993) q[0];
sx q[0];
rz(-0.81863716) q[0];
sx q[0];
rz(-0.37055117) q[0];
x q[1];
rz(-1.0867811) q[2];
sx q[2];
rz(-0.40905158) q[2];
sx q[2];
rz(1.3889927) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5999509) q[1];
sx q[1];
rz(-2.5516918) q[1];
sx q[1];
rz(1.4093424) q[1];
rz(0.67477711) q[3];
sx q[3];
rz(-2.1222847) q[3];
sx q[3];
rz(-1.2434208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1696986) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(0.92740721) q[2];
rz(1.5291322) q[3];
sx q[3];
rz(-1.1011139) q[3];
sx q[3];
rz(2.7839938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9277495) q[0];
sx q[0];
rz(-2.1723211) q[0];
sx q[0];
rz(-1.7359605) q[0];
rz(-1.7270145) q[1];
sx q[1];
rz(-2.3443293) q[1];
sx q[1];
rz(-0.79967156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1522923) q[0];
sx q[0];
rz(-1.5591027) q[0];
sx q[0];
rz(1.508827) q[0];
rz(-pi) q[1];
rz(-2.7437216) q[2];
sx q[2];
rz(-0.33411538) q[2];
sx q[2];
rz(2.1159005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.68467632) q[1];
sx q[1];
rz(-1.231319) q[1];
sx q[1];
rz(2.875332) q[1];
x q[2];
rz(-1.762885) q[3];
sx q[3];
rz(-1.9682765) q[3];
sx q[3];
rz(-0.31645889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.25524011) q[2];
sx q[2];
rz(-2.4031874) q[2];
sx q[2];
rz(2.3389471) q[2];
rz(-2.8211527) q[3];
sx q[3];
rz(-1.8010062) q[3];
sx q[3];
rz(-1.505544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36444148) q[0];
sx q[0];
rz(-3.1138804) q[0];
sx q[0];
rz(-0.030315422) q[0];
rz(-1.3069356) q[1];
sx q[1];
rz(-2.0590643) q[1];
sx q[1];
rz(2.511715) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8960285) q[0];
sx q[0];
rz(-1.7587852) q[0];
sx q[0];
rz(-0.4281559) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81923072) q[2];
sx q[2];
rz(-2.5360245) q[2];
sx q[2];
rz(-1.8964963) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83104937) q[1];
sx q[1];
rz(-1.2782974) q[1];
sx q[1];
rz(-1.7834375) q[1];
x q[2];
rz(-2.9050352) q[3];
sx q[3];
rz(-0.45522296) q[3];
sx q[3];
rz(2.5270324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1097172) q[2];
sx q[2];
rz(-0.66706359) q[2];
sx q[2];
rz(1.9943705) q[2];
rz(0.68754292) q[3];
sx q[3];
rz(-0.59194982) q[3];
sx q[3];
rz(-1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3654093) q[0];
sx q[0];
rz(-2.9054346) q[0];
sx q[0];
rz(-2.6614406) q[0];
rz(2.7393553) q[1];
sx q[1];
rz(-1.6419342) q[1];
sx q[1];
rz(-0.38280907) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3881809) q[0];
sx q[0];
rz(-2.7482902) q[0];
sx q[0];
rz(2.1476782) q[0];
rz(-pi) q[1];
rz(-1.9391962) q[2];
sx q[2];
rz(-2.830626) q[2];
sx q[2];
rz(-2.5045365) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3814712) q[1];
sx q[1];
rz(-2.6497726) q[1];
sx q[1];
rz(0.90729883) q[1];
x q[2];
rz(1.6005101) q[3];
sx q[3];
rz(-1.374657) q[3];
sx q[3];
rz(-0.33523424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39168921) q[2];
sx q[2];
rz(-2.4112371) q[2];
sx q[2];
rz(-2.7833617) q[2];
rz(-0.41424888) q[3];
sx q[3];
rz(-1.0284938) q[3];
sx q[3];
rz(2.7481368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8527894) q[0];
sx q[0];
rz(-1.8599334) q[0];
sx q[0];
rz(2.7237256) q[0];
rz(-1.4990384) q[1];
sx q[1];
rz(-2.7211029) q[1];
sx q[1];
rz(2.5299759) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1660991) q[0];
sx q[0];
rz(-2.8320441) q[0];
sx q[0];
rz(0.73348372) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29199227) q[2];
sx q[2];
rz(-0.21459178) q[2];
sx q[2];
rz(0.24927441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4222718) q[1];
sx q[1];
rz(-0.55972717) q[1];
sx q[1];
rz(1.0258963) q[1];
rz(-pi) q[2];
rz(2.9183594) q[3];
sx q[3];
rz(-2.7866552) q[3];
sx q[3];
rz(-2.8508773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.186782) q[2];
sx q[2];
rz(-1.6839226) q[2];
sx q[2];
rz(0.35624722) q[2];
rz(-2.6567843) q[3];
sx q[3];
rz(-1.3508513) q[3];
sx q[3];
rz(-0.029732186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271069) q[0];
sx q[0];
rz(-2.0368545) q[0];
sx q[0];
rz(1.4622965) q[0];
rz(-0.38743427) q[1];
sx q[1];
rz(-0.92715627) q[1];
sx q[1];
rz(-0.80518728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1925404) q[0];
sx q[0];
rz(-2.1645191) q[0];
sx q[0];
rz(-0.70021061) q[0];
x q[1];
rz(2.3024955) q[2];
sx q[2];
rz(-0.57575127) q[2];
sx q[2];
rz(-2.8481399) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.78811121) q[1];
sx q[1];
rz(-2.821021) q[1];
sx q[1];
rz(-2.660257) q[1];
rz(-1.8191387) q[3];
sx q[3];
rz(-1.2183918) q[3];
sx q[3];
rz(-0.067841522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8055973) q[2];
sx q[2];
rz(-2.5184641) q[2];
sx q[2];
rz(1.222329) q[2];
rz(-0.81021106) q[3];
sx q[3];
rz(-2.2170292) q[3];
sx q[3];
rz(-1.5626102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41351086) q[0];
sx q[0];
rz(-1.6376729) q[0];
sx q[0];
rz(-1.6117657) q[0];
rz(1.275508) q[1];
sx q[1];
rz(-2.7559912) q[1];
sx q[1];
rz(-1.4082946) q[1];
rz(-1.8985844) q[2];
sx q[2];
rz(-0.40400436) q[2];
sx q[2];
rz(-1.8346471) q[2];
rz(2.2386407) q[3];
sx q[3];
rz(-1.2085689) q[3];
sx q[3];
rz(2.989789) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
