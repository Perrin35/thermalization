OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6025699) q[0];
sx q[0];
rz(-0.56350001) q[0];
sx q[0];
rz(-2.6846057) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(1.2759804) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39014434) q[0];
sx q[0];
rz(-1.7062618) q[0];
sx q[0];
rz(-0.22060237) q[0];
x q[1];
rz(1.4957501) q[2];
sx q[2];
rz(-1.6610166) q[2];
sx q[2];
rz(0.98352945) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1162122) q[1];
sx q[1];
rz(-2.930021) q[1];
sx q[1];
rz(-1.8616574) q[1];
x q[2];
rz(2.2184578) q[3];
sx q[3];
rz(-1.2520408) q[3];
sx q[3];
rz(1.8766599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(-1.6312381) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(-0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(-0.045036137) q[0];
rz(-0.33915195) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(-2.3388458) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0720285) q[0];
sx q[0];
rz(-1.9575798) q[0];
sx q[0];
rz(-0.2441498) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0375508) q[2];
sx q[2];
rz(-2.0718956) q[2];
sx q[2];
rz(2.3061371) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74233494) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(1.4820815) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4071776) q[3];
sx q[3];
rz(-0.67577261) q[3];
sx q[3];
rz(2.8676652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37614432) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(-3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95214343) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.3396938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08683603) q[0];
sx q[0];
rz(-1.9924379) q[0];
sx q[0];
rz(2.3479168) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2604733) q[2];
sx q[2];
rz(-1.3291877) q[2];
sx q[2];
rz(-1.2348262) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9078319) q[1];
sx q[1];
rz(-1.7572174) q[1];
sx q[1];
rz(-0.16252653) q[1];
rz(-pi) q[2];
rz(0.6891059) q[3];
sx q[3];
rz(-1.2561241) q[3];
sx q[3];
rz(-2.2113707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(-1.408067) q[2];
rz(0.18313289) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6614439) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(-0.21425042) q[0];
rz(-2.7125773) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(-0.66657153) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65581215) q[0];
sx q[0];
rz(-2.862145) q[0];
sx q[0];
rz(0.73894545) q[0];
rz(-pi) q[1];
rz(2.4273657) q[2];
sx q[2];
rz(-1.7923317) q[2];
sx q[2];
rz(0.54745882) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2506634) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(-2.3198747) q[1];
rz(-pi) q[2];
rz(2.9922819) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(-0.028545054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8551222) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(0.66037035) q[2];
rz(1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-2.3045325) q[0];
rz(2.191026) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(-1.9821092) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75393049) q[0];
sx q[0];
rz(-1.2571063) q[0];
sx q[0];
rz(0.73648209) q[0];
rz(-1.3443089) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(-0.047705334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78427181) q[1];
sx q[1];
rz(-1.9958152) q[1];
sx q[1];
rz(2.1678863) q[1];
x q[2];
rz(0.67622185) q[3];
sx q[3];
rz(-0.75934221) q[3];
sx q[3];
rz(-0.29952213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(-1.8699899) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(-0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5746675) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(-1.1151474) q[0];
rz(0.043958157) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(-3.040722) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6504165) q[0];
sx q[0];
rz(-2.4917779) q[0];
sx q[0];
rz(-2.8894436) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8632554) q[2];
sx q[2];
rz(-1.1012226) q[2];
sx q[2];
rz(1.3172319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0724153) q[1];
sx q[1];
rz(-2.8848856) q[1];
sx q[1];
rz(-2.2290345) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40525873) q[3];
sx q[3];
rz(-1.7856981) q[3];
sx q[3];
rz(0.98291558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(2.960079) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.293175) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(2.4538453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6531658) q[0];
sx q[0];
rz(-1.7894206) q[0];
sx q[0];
rz(1.2290918) q[0];
x q[1];
rz(-2.9155832) q[2];
sx q[2];
rz(-1.8812211) q[2];
sx q[2];
rz(1.4357391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2076599) q[1];
sx q[1];
rz(-1.7001517) q[1];
sx q[1];
rz(-1.5473066) q[1];
rz(-pi) q[2];
rz(-1.5958105) q[3];
sx q[3];
rz(-2.5731312) q[3];
sx q[3];
rz(2.3908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8035651) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(0.98085105) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5279609) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(-3.0623073) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.1987196) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0543025) q[0];
sx q[0];
rz(-1.1371005) q[0];
sx q[0];
rz(-1.646423) q[0];
rz(-1.2077246) q[2];
sx q[2];
rz(-1.4903526) q[2];
sx q[2];
rz(2.9094839) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5115576) q[1];
sx q[1];
rz(-0.50543284) q[1];
sx q[1];
rz(0.10314421) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49097455) q[3];
sx q[3];
rz(-2.6262865) q[3];
sx q[3];
rz(-1.6083503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2856059) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(-0.78906995) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.566074) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(-2.5532706) q[0];
rz(-0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8342455) q[0];
sx q[0];
rz(-1.630162) q[0];
sx q[0];
rz(-0.10660118) q[0];
rz(-1.9503715) q[2];
sx q[2];
rz(-1.8827202) q[2];
sx q[2];
rz(-1.003007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0837299) q[1];
sx q[1];
rz(-1.607556) q[1];
sx q[1];
rz(2.4576393) q[1];
rz(-1.9282354) q[3];
sx q[3];
rz(-0.94944984) q[3];
sx q[3];
rz(-0.77222094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.062722) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(0.36994568) q[2];
rz(-0.89921078) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-2.488234) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(-2.9577589) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606124) q[0];
sx q[0];
rz(-1.4502589) q[0];
sx q[0];
rz(-1.6880077) q[0];
x q[1];
rz(0.75283639) q[2];
sx q[2];
rz(-0.92217731) q[2];
sx q[2];
rz(-0.93968771) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76845968) q[1];
sx q[1];
rz(-0.8994973) q[1];
sx q[1];
rz(-0.21233227) q[1];
x q[2];
rz(-3.1179908) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(0.49707801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.9188149) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(-2.9392021) q[2];
rz(-2.0683794) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595173) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(-0.36322414) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(-2.4167378) q[2];
sx q[2];
rz(-0.58314322) q[2];
sx q[2];
rz(-2.4287139) q[2];
rz(-0.98060676) q[3];
sx q[3];
rz(-1.4242314) q[3];
sx q[3];
rz(-0.84606708) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
