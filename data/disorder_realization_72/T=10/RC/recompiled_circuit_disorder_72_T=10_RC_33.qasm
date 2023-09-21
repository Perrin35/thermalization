OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(1.7680661) q[0];
sx q[0];
rz(11.058523) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(3.9197796) q[1];
sx q[1];
rz(9.9240886) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0190174) q[0];
sx q[0];
rz(-0.28261533) q[0];
sx q[0];
rz(-0.22960381) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7668031) q[2];
sx q[2];
rz(-0.73050806) q[2];
sx q[2];
rz(-1.0175878) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3589188) q[1];
sx q[1];
rz(-0.75190699) q[1];
sx q[1];
rz(2.9208675) q[1];
rz(-pi) q[2];
rz(-1.4673759) q[3];
sx q[3];
rz(-1.5442344) q[3];
sx q[3];
rz(2.6288222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(-2.1271558) q[2];
rz(-0.23400083) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(2.8570989) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(-1.6218119) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(1.0027592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53045814) q[0];
sx q[0];
rz(-0.76871745) q[0];
sx q[0];
rz(2.4918633) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1728671) q[2];
sx q[2];
rz(-1.7555408) q[2];
sx q[2];
rz(0.52344054) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1042852) q[1];
sx q[1];
rz(-0.98854317) q[1];
sx q[1];
rz(-2.9557455) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5428883) q[3];
sx q[3];
rz(-1.8667392) q[3];
sx q[3];
rz(-2.242089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(-2.9906452) q[2];
rz(-0.41444591) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(0.088236563) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8931483) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(0.88358203) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4602063) q[0];
sx q[0];
rz(-1.1707414) q[0];
sx q[0];
rz(-1.7494739) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29404624) q[2];
sx q[2];
rz(-1.9532734) q[2];
sx q[2];
rz(1.2750212) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5501432) q[1];
sx q[1];
rz(-1.5482229) q[1];
sx q[1];
rz(1.8557465) q[1];
rz(-pi) q[2];
rz(2.9900842) q[3];
sx q[3];
rz(-0.66473367) q[3];
sx q[3];
rz(0.96707771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(2.1515576) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(0.86301962) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54116762) q[0];
sx q[0];
rz(-0.74099243) q[0];
sx q[0];
rz(0.11359544) q[0];
rz(-pi) q[1];
rz(2.8407211) q[2];
sx q[2];
rz(-0.30756018) q[2];
sx q[2];
rz(3.0645264) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6039227) q[1];
sx q[1];
rz(-0.99891716) q[1];
sx q[1];
rz(-0.13725431) q[1];
rz(0.22927852) q[3];
sx q[3];
rz(-1.5203272) q[3];
sx q[3];
rz(0.63427395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.52175534) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(3.0299305) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6977285) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(2.2648947) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278496) q[0];
sx q[0];
rz(-0.62999524) q[0];
sx q[0];
rz(1.4907452) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5259597) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(0.63477883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1716869) q[1];
sx q[1];
rz(-1.415442) q[1];
sx q[1];
rz(-3.1093662) q[1];
rz(2.3239273) q[3];
sx q[3];
rz(-0.92901232) q[3];
sx q[3];
rz(1.6587917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-2.6780224) q[2];
rz(-0.56435895) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(0.92818964) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1148465) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-0.39658305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81178938) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(-2.0407709) q[0];
x q[1];
rz(-0.47607143) q[2];
sx q[2];
rz(-1.2928315) q[2];
sx q[2];
rz(-1.6518041) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0551758) q[1];
sx q[1];
rz(-2.2284818) q[1];
sx q[1];
rz(0.1187101) q[1];
rz(-pi) q[2];
rz(1.3271689) q[3];
sx q[3];
rz(-2.1414087) q[3];
sx q[3];
rz(2.7057196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(-1.5552103) q[2];
rz(1.4554626) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7431188) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(-0.78480762) q[0];
rz(-1.8709042) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(2.0152337) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9046017) q[0];
sx q[0];
rz(-0.31123529) q[0];
sx q[0];
rz(2.8471332) q[0];
x q[1];
rz(-1.8225841) q[2];
sx q[2];
rz(-2.6195824) q[2];
sx q[2];
rz(0.28282794) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4451051) q[1];
sx q[1];
rz(-1.3828053) q[1];
sx q[1];
rz(-2.6871215) q[1];
rz(-pi) q[2];
rz(1.2231636) q[3];
sx q[3];
rz(-1.4813444) q[3];
sx q[3];
rz(-1.4060494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.209098) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37041935) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(0.11288189) q[0];
rz(2.1408634) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(3.1088366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927133) q[0];
sx q[0];
rz(-1.8939051) q[0];
sx q[0];
rz(-2.8723641) q[0];
rz(-pi) q[1];
rz(0.92808) q[2];
sx q[2];
rz(-0.66814458) q[2];
sx q[2];
rz(1.1754787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2331088) q[1];
sx q[1];
rz(-1.9875803) q[1];
sx q[1];
rz(-3.0770739) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.950374) q[3];
sx q[3];
rz(-1.4879585) q[3];
sx q[3];
rz(-0.4656725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1774566) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-2.5194871) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416097) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(-0.73293066) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(-2.656235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1982614) q[0];
sx q[0];
rz(-1.7051464) q[0];
sx q[0];
rz(1.5166548) q[0];
rz(1.4114686) q[2];
sx q[2];
rz(-2.1630641) q[2];
sx q[2];
rz(1.9221905) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24385142) q[1];
sx q[1];
rz(-0.85695367) q[1];
sx q[1];
rz(2.7276911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58166196) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(-1.7285085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0316524) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(-1.4032646) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(1.2592978) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3132968) q[0];
sx q[0];
rz(-1.3643571) q[0];
sx q[0];
rz(0.22160991) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7190476) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(1.6809747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1640537) q[1];
sx q[1];
rz(-2.7643449) q[1];
sx q[1];
rz(1.8305199) q[1];
x q[2];
rz(1.7353021) q[3];
sx q[3];
rz(-1.4146155) q[3];
sx q[3];
rz(-0.7848878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(2.4895978) q[2];
rz(-2.5478798) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020637) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(-1.2364173) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(2.9700206) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(-2.0486352) q[3];
sx q[3];
rz(-2.408705) q[3];
sx q[3];
rz(-0.68791289) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
