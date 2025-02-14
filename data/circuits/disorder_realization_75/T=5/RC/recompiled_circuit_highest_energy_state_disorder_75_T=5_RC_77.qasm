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
rz(-1.8259814) q[0];
sx q[0];
rz(-2.313518) q[0];
sx q[0];
rz(-2.2928884) q[0];
rz(-1.6500213) q[1];
sx q[1];
rz(-2.4529011) q[1];
sx q[1];
rz(-0.42660776) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4147581) q[0];
sx q[0];
rz(-1.4953185) q[0];
sx q[0];
rz(-0.10087905) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6831065) q[2];
sx q[2];
rz(-1.4698311) q[2];
sx q[2];
rz(-1.8835406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1060488) q[1];
sx q[1];
rz(-1.1263761) q[1];
sx q[1];
rz(1.4900833) q[1];
x q[2];
rz(2.3412104) q[3];
sx q[3];
rz(-2.0315779) q[3];
sx q[3];
rz(-2.9886101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9244869) q[2];
sx q[2];
rz(-1.3873528) q[2];
sx q[2];
rz(0.039487751) q[2];
rz(-2.110179) q[3];
sx q[3];
rz(-1.02905) q[3];
sx q[3];
rz(1.2379117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51204387) q[0];
sx q[0];
rz(-1.5271674) q[0];
sx q[0];
rz(1.7426096) q[0];
rz(1.6276739) q[1];
sx q[1];
rz(-0.84295034) q[1];
sx q[1];
rz(-1.4167851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2606974) q[0];
sx q[0];
rz(-1.7923549) q[0];
sx q[0];
rz(0.76789121) q[0];
x q[1];
rz(-2.3624647) q[2];
sx q[2];
rz(-1.7567524) q[2];
sx q[2];
rz(0.27006876) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9263404) q[1];
sx q[1];
rz(-1.7071586) q[1];
sx q[1];
rz(0.3991613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7339277) q[3];
sx q[3];
rz(-1.614794) q[3];
sx q[3];
rz(3.0394579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17025718) q[2];
sx q[2];
rz(-2.038326) q[2];
sx q[2];
rz(0.7106759) q[2];
rz(-3.1273048) q[3];
sx q[3];
rz(-0.24737869) q[3];
sx q[3];
rz(1.8054731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22064848) q[0];
sx q[0];
rz(-0.9592239) q[0];
sx q[0];
rz(2.7253819) q[0];
rz(1.8572218) q[1];
sx q[1];
rz(-1.4764079) q[1];
sx q[1];
rz(-2.823337) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98781768) q[0];
sx q[0];
rz(-1.1889699) q[0];
sx q[0];
rz(2.8433958) q[0];
rz(-0.84789611) q[2];
sx q[2];
rz(-1.0506786) q[2];
sx q[2];
rz(1.073357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-15/(11*pi)) q[1];
sx q[1];
rz(-2.4726082) q[1];
sx q[1];
rz(-2.9218052) q[1];
rz(-pi) q[2];
rz(1.8793568) q[3];
sx q[3];
rz(-1.3342382) q[3];
sx q[3];
rz(-2.6738338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2661813) q[2];
sx q[2];
rz(-2.3616932) q[2];
sx q[2];
rz(1.2904588) q[2];
rz(0.053704638) q[3];
sx q[3];
rz(-1.3196245) q[3];
sx q[3];
rz(-1.7558338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5586435) q[0];
sx q[0];
rz(-0.68840331) q[0];
sx q[0];
rz(-2.131856) q[0];
rz(0.91521493) q[1];
sx q[1];
rz(-1.5292294) q[1];
sx q[1];
rz(-0.42542747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64176004) q[0];
sx q[0];
rz(-1.7374037) q[0];
sx q[0];
rz(-2.2918022) q[0];
x q[1];
rz(3.0680823) q[2];
sx q[2];
rz(-0.82907721) q[2];
sx q[2];
rz(-2.4501767) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3784619) q[1];
sx q[1];
rz(-1.1647071) q[1];
sx q[1];
rz(2.2927845) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0056267319) q[3];
sx q[3];
rz(-1.4079908) q[3];
sx q[3];
rz(2.4037647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90011826) q[2];
sx q[2];
rz(-1.590531) q[2];
sx q[2];
rz(-3.0305064) q[2];
rz(0.19242081) q[3];
sx q[3];
rz(-2.7399053) q[3];
sx q[3];
rz(0.69453159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38839328) q[0];
sx q[0];
rz(-2.41112) q[0];
sx q[0];
rz(-0.86517349) q[0];
rz(2.6490037) q[1];
sx q[1];
rz(-1.0961696) q[1];
sx q[1];
rz(-1.385484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2619721) q[0];
sx q[0];
rz(-0.79085717) q[0];
sx q[0];
rz(-1.3778995) q[0];
rz(-pi) q[1];
rz(-0.4200468) q[2];
sx q[2];
rz(-0.85424747) q[2];
sx q[2];
rz(1.4619712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4069479) q[1];
sx q[1];
rz(-2.4666767) q[1];
sx q[1];
rz(-0.65924834) q[1];
rz(-2.288743) q[3];
sx q[3];
rz(-0.88969031) q[3];
sx q[3];
rz(1.8913325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2546786) q[2];
sx q[2];
rz(-0.46923894) q[2];
sx q[2];
rz(2.4349507) q[2];
rz(-2.8142269) q[3];
sx q[3];
rz(-1.8492161) q[3];
sx q[3];
rz(0.65739003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3051598) q[0];
sx q[0];
rz(-1.7921472) q[0];
sx q[0];
rz(2.7850372) q[0];
rz(-2.5794079) q[1];
sx q[1];
rz(-1.1327344) q[1];
sx q[1];
rz(-2.1551989) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0092457744) q[0];
sx q[0];
rz(-2.2725687) q[0];
sx q[0];
rz(2.470507) q[0];
x q[1];
rz(-2.7964727) q[2];
sx q[2];
rz(-1.2212996) q[2];
sx q[2];
rz(2.0728759) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74979177) q[1];
sx q[1];
rz(-1.5967249) q[1];
sx q[1];
rz(2.3830929) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3273193) q[3];
sx q[3];
rz(-1.3946574) q[3];
sx q[3];
rz(-1.0276664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76546136) q[2];
sx q[2];
rz(-2.4262846) q[2];
sx q[2];
rz(2.4410655) q[2];
rz(-2.1721407) q[3];
sx q[3];
rz(-1.0337831) q[3];
sx q[3];
rz(-0.9303003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9483865) q[0];
sx q[0];
rz(-2.8128615) q[0];
sx q[0];
rz(-2.9902003) q[0];
rz(-1.5104177) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(1.4600533) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23510012) q[0];
sx q[0];
rz(-2.0564046) q[0];
sx q[0];
rz(2.4686345) q[0];
rz(-pi) q[1];
rz(0.362927) q[2];
sx q[2];
rz(-1.0725601) q[2];
sx q[2];
rz(1.0403596) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9402079) q[1];
sx q[1];
rz(-1.0321069) q[1];
sx q[1];
rz(-2.687665) q[1];
x q[2];
rz(1.7655108) q[3];
sx q[3];
rz(-1.7140183) q[3];
sx q[3];
rz(-2.9778874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2363362) q[2];
sx q[2];
rz(-2.4589804) q[2];
sx q[2];
rz(-2.7899817) q[2];
rz(-0.44175276) q[3];
sx q[3];
rz(-0.92919246) q[3];
sx q[3];
rz(2.9898804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0996967) q[0];
sx q[0];
rz(-1.8459039) q[0];
sx q[0];
rz(1.1610485) q[0];
rz(-2.4940122) q[1];
sx q[1];
rz(-1.3787965) q[1];
sx q[1];
rz(-2.3191998) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1294382) q[0];
sx q[0];
rz(-2.4732051) q[0];
sx q[0];
rz(0.8968312) q[0];
x q[1];
rz(-0.47745173) q[2];
sx q[2];
rz(-2.3465354) q[2];
sx q[2];
rz(0.24518046) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97473652) q[1];
sx q[1];
rz(-0.8203764) q[1];
sx q[1];
rz(-1.1715164) q[1];
x q[2];
rz(-1.137072) q[3];
sx q[3];
rz(-1.9334643) q[3];
sx q[3];
rz(0.39946242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9245727) q[2];
sx q[2];
rz(-1.9968888) q[2];
sx q[2];
rz(1.8355969) q[2];
rz(-0.71803391) q[3];
sx q[3];
rz(-2.3170203) q[3];
sx q[3];
rz(0.048862783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43593916) q[0];
sx q[0];
rz(-2.0581364) q[0];
sx q[0];
rz(-0.017024592) q[0];
rz(1.1964993) q[1];
sx q[1];
rz(-2.7462609) q[1];
sx q[1];
rz(0.33504018) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8340613) q[0];
sx q[0];
rz(-2.1387707) q[0];
sx q[0];
rz(3.1295883) q[0];
rz(-pi) q[1];
rz(-0.52729221) q[2];
sx q[2];
rz(-2.4278473) q[2];
sx q[2];
rz(-1.1251768) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0373203) q[1];
sx q[1];
rz(-1.3824711) q[1];
sx q[1];
rz(2.2467747) q[1];
rz(-pi) q[2];
rz(-0.7850168) q[3];
sx q[3];
rz(-2.5015273) q[3];
sx q[3];
rz(-0.16798108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42387858) q[2];
sx q[2];
rz(-0.7462036) q[2];
sx q[2];
rz(-2.8625873) q[2];
rz(1.772607) q[3];
sx q[3];
rz(-1.5149346) q[3];
sx q[3];
rz(0.92931187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5170711) q[0];
sx q[0];
rz(-0.14685024) q[0];
sx q[0];
rz(2.3745234) q[0];
rz(-0.37297878) q[1];
sx q[1];
rz(-1.6423128) q[1];
sx q[1];
rz(-2.9708718) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0559346) q[0];
sx q[0];
rz(-1.5996337) q[0];
sx q[0];
rz(2.4066854) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4095979) q[2];
sx q[2];
rz(-1.9809763) q[2];
sx q[2];
rz(2.8267677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4258044) q[1];
sx q[1];
rz(-1.3223151) q[1];
sx q[1];
rz(1.2475509) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7307698) q[3];
sx q[3];
rz(-2.0334646) q[3];
sx q[3];
rz(1.8335952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3501222) q[2];
sx q[2];
rz(-0.56362027) q[2];
sx q[2];
rz(-0.0027837022) q[2];
rz(2.2400098) q[3];
sx q[3];
rz(-2.1193347) q[3];
sx q[3];
rz(0.80488747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053454178) q[0];
sx q[0];
rz(-2.8202941) q[0];
sx q[0];
rz(1.1711076) q[0];
rz(2.3114655) q[1];
sx q[1];
rz(-1.3273888) q[1];
sx q[1];
rz(-1.8524016) q[1];
rz(-0.89494643) q[2];
sx q[2];
rz(-2.4104626) q[2];
sx q[2];
rz(1.8117803) q[2];
rz(2.2318673) q[3];
sx q[3];
rz(-2.2561938) q[3];
sx q[3];
rz(1.6267226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
