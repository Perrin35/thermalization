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
rz(-2.4515406) q[0];
sx q[0];
rz(-2.2909988) q[0];
sx q[0];
rz(2.9414862) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(5.7601647) q[1];
sx q[1];
rz(10.756607) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9560741) q[0];
sx q[0];
rz(-1.3393434) q[0];
sx q[0];
rz(-1.1197907) q[0];
x q[1];
rz(0.66352377) q[2];
sx q[2];
rz(-0.71394701) q[2];
sx q[2];
rz(-0.51807846) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9871949) q[1];
sx q[1];
rz(-1.618865) q[1];
sx q[1];
rz(-2.4632719) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0947919) q[3];
sx q[3];
rz(-1.4972357) q[3];
sx q[3];
rz(-1.5293763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0693822) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(1.3145831) q[2];
rz(-0.9032816) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8697934) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(-0.85579175) q[0];
rz(0.77404147) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(-2.3537297) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6928394) q[0];
sx q[0];
rz(-1.3200545) q[0];
sx q[0];
rz(3.1186799) q[0];
x q[1];
rz(-3.0237979) q[2];
sx q[2];
rz(-2.795145) q[2];
sx q[2];
rz(0.36855498) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98709244) q[1];
sx q[1];
rz(-2.5804511) q[1];
sx q[1];
rz(2.2461056) q[1];
rz(-pi) q[2];
rz(1.3608906) q[3];
sx q[3];
rz(-2.2271101) q[3];
sx q[3];
rz(-0.24649749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0756125) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(0.39609972) q[2];
rz(1.5710477) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1876672) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(3.0809825) q[0];
rz(-2.4568779) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(2.8025467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0051088) q[0];
sx q[0];
rz(-1.5500796) q[0];
sx q[0];
rz(2.5584975) q[0];
rz(0.18908638) q[2];
sx q[2];
rz(-1.4497309) q[2];
sx q[2];
rz(2.6775286) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78279692) q[1];
sx q[1];
rz(-0.41886815) q[1];
sx q[1];
rz(-2.2936506) q[1];
rz(-pi) q[2];
rz(-2.364758) q[3];
sx q[3];
rz(-1.5791487) q[3];
sx q[3];
rz(1.9296822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6591349) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(-2.6711312) q[2];
rz(2.2792234) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7793133) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(2.733574) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(-1.8203576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7119638) q[0];
sx q[0];
rz(-2.5416059) q[0];
sx q[0];
rz(-0.55498755) q[0];
x q[1];
rz(-0.086060103) q[2];
sx q[2];
rz(-2.2241908) q[2];
sx q[2];
rz(-1.6376405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3258354) q[1];
sx q[1];
rz(-1.6920631) q[1];
sx q[1];
rz(-1.497529) q[1];
rz(-1.9923237) q[3];
sx q[3];
rz(-1.9644004) q[3];
sx q[3];
rz(-0.43962653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9185751) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(0.72285405) q[2];
rz(0.84960788) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-2.6246922) q[0];
sx q[0];
rz(-2.5816259) q[0];
rz(0.2050744) q[1];
sx q[1];
rz(-1.8849025) q[1];
sx q[1];
rz(-2.8505039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2084853) q[0];
sx q[0];
rz(-1.4015084) q[0];
sx q[0];
rz(-0.27219682) q[0];
rz(-pi) q[1];
x q[1];
rz(0.061441378) q[2];
sx q[2];
rz(-0.80020088) q[2];
sx q[2];
rz(0.35129181) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4704125) q[1];
sx q[1];
rz(-2.414783) q[1];
sx q[1];
rz(1.2695168) q[1];
rz(-pi) q[2];
rz(1.4104615) q[3];
sx q[3];
rz(-2.3162127) q[3];
sx q[3];
rz(1.3813409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9186972) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(0.29829868) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(-1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30763141) q[0];
sx q[0];
rz(-1.6809604) q[0];
sx q[0];
rz(-2.5391915) q[0];
rz(-0.3715474) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(2.1098302) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1790947) q[0];
sx q[0];
rz(-2.966605) q[0];
sx q[0];
rz(1.6618769) q[0];
rz(-pi) q[1];
rz(0.40470064) q[2];
sx q[2];
rz(-1.3794624) q[2];
sx q[2];
rz(-1.9218685) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18373016) q[1];
sx q[1];
rz(-1.0938083) q[1];
sx q[1];
rz(-1.7836003) q[1];
x q[2];
rz(0.86182819) q[3];
sx q[3];
rz(-0.78360451) q[3];
sx q[3];
rz(2.6605822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3370257) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(1.5186914) q[2];
rz(0.8484146) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(-0.039610473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1161716) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(-0.99789944) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.7707228) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3034015) q[0];
sx q[0];
rz(-0.63119054) q[0];
sx q[0];
rz(3.1188008) q[0];
rz(-pi) q[1];
rz(-1.9657126) q[2];
sx q[2];
rz(-0.38287258) q[2];
sx q[2];
rz(-1.9667786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8805367) q[1];
sx q[1];
rz(-2.4789101) q[1];
sx q[1];
rz(0.39775325) q[1];
rz(2.3418535) q[3];
sx q[3];
rz(-2.1465786) q[3];
sx q[3];
rz(2.6588001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.85713282) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(0.029646309) q[2];
rz(0.61521411) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(-2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089559473) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(2.3585368) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(1.7436183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8496765) q[0];
sx q[0];
rz(-1.9789985) q[0];
sx q[0];
rz(1.085039) q[0];
rz(-pi) q[1];
rz(-0.88770788) q[2];
sx q[2];
rz(-2.7811433) q[2];
sx q[2];
rz(0.53149022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3191171) q[1];
sx q[1];
rz(-0.90293316) q[1];
sx q[1];
rz(-3.130359) q[1];
x q[2];
rz(1.3574202) q[3];
sx q[3];
rz(-0.7979352) q[3];
sx q[3];
rz(-2.7415581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29584259) q[2];
sx q[2];
rz(-2.5551899) q[2];
sx q[2];
rz(1.3207377) q[2];
rz(-2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(-2.1871908) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(1.0844213) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0823776) q[0];
sx q[0];
rz(-2.431916) q[0];
sx q[0];
rz(-1.6355455) q[0];
rz(-pi) q[1];
rz(2.886359) q[2];
sx q[2];
rz(-1.298578) q[2];
sx q[2];
rz(-1.24025) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3855648) q[1];
sx q[1];
rz(-1.4984908) q[1];
sx q[1];
rz(1.0977618) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74862759) q[3];
sx q[3];
rz(-1.4350865) q[3];
sx q[3];
rz(2.3627797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.083100975) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(-2.7044738) q[2];
rz(-1.3433749) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042353543) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(1.4105256) q[0];
rz(2.9208185) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-0.9332307) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8932874) q[0];
sx q[0];
rz(-2.3961005) q[0];
sx q[0];
rz(3.1214691) q[0];
x q[1];
rz(1.2104697) q[2];
sx q[2];
rz(-0.98305741) q[2];
sx q[2];
rz(0.34229842) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5367914) q[1];
sx q[1];
rz(-1.7030506) q[1];
sx q[1];
rz(-2.3459646) q[1];
rz(2.203412) q[3];
sx q[3];
rz(-1.4804192) q[3];
sx q[3];
rz(2.602102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.431388) q[2];
sx q[2];
rz(-2.7993671) q[2];
sx q[2];
rz(0.63146511) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5030293) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(-1.6497816) q[2];
sx q[2];
rz(-2.7597703) q[2];
sx q[2];
rz(-2.4673354) q[2];
rz(-2.9819103) q[3];
sx q[3];
rz(-0.57542141) q[3];
sx q[3];
rz(0.66948359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
