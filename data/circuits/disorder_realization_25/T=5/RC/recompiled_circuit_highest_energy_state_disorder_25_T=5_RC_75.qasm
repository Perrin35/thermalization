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
rz(2.202232) q[0];
sx q[0];
rz(-1.1190177) q[0];
sx q[0];
rz(3.1079259) q[0];
rz(1.6839924) q[1];
sx q[1];
rz(-1.3352609) q[1];
sx q[1];
rz(-1.508498) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38070883) q[0];
sx q[0];
rz(-0.11822001) q[0];
sx q[0];
rz(1.8272754) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.066285) q[2];
sx q[2];
rz(-1.2774373) q[2];
sx q[2];
rz(-1.2874029) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7693387) q[1];
sx q[1];
rz(-0.70972356) q[1];
sx q[1];
rz(-2.918887) q[1];
x q[2];
rz(-2.331459) q[3];
sx q[3];
rz(-2.0647233) q[3];
sx q[3];
rz(2.877748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17026751) q[2];
sx q[2];
rz(-0.89189947) q[2];
sx q[2];
rz(0.83667052) q[2];
rz(-0.80638805) q[3];
sx q[3];
rz(-0.92525768) q[3];
sx q[3];
rz(-0.18536082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065501325) q[0];
sx q[0];
rz(-1.0668904) q[0];
sx q[0];
rz(1.2170894) q[0];
rz(-0.20054664) q[1];
sx q[1];
rz(-2.3417818) q[1];
sx q[1];
rz(-0.54380551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1233305) q[0];
sx q[0];
rz(-1.2549434) q[0];
sx q[0];
rz(-1.8118152) q[0];
x q[1];
rz(3.1214525) q[2];
sx q[2];
rz(-0.59505645) q[2];
sx q[2];
rz(-0.12660566) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75419352) q[1];
sx q[1];
rz(-0.48739823) q[1];
sx q[1];
rz(2.2226366) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1065349) q[3];
sx q[3];
rz(-2.2457099) q[3];
sx q[3];
rz(-1.1763151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66889888) q[2];
sx q[2];
rz(-0.74287477) q[2];
sx q[2];
rz(-1.6061858) q[2];
rz(-1.499048) q[3];
sx q[3];
rz(-0.58321548) q[3];
sx q[3];
rz(-1.8892052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1010308) q[0];
sx q[0];
rz(-0.56483785) q[0];
sx q[0];
rz(0.0095796883) q[0];
rz(1.0293695) q[1];
sx q[1];
rz(-1.5734438) q[1];
sx q[1];
rz(1.8052489) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9294465) q[0];
sx q[0];
rz(-0.98185724) q[0];
sx q[0];
rz(0.94561823) q[0];
rz(-2.1161386) q[2];
sx q[2];
rz(-0.57078123) q[2];
sx q[2];
rz(-0.88583835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.0071466999) q[1];
sx q[1];
rz(-1.916051) q[1];
sx q[1];
rz(2.9551639) q[1];
rz(-pi) q[2];
rz(-1.1824781) q[3];
sx q[3];
rz(-0.42436436) q[3];
sx q[3];
rz(-2.1774815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0336527) q[2];
sx q[2];
rz(-2.630271) q[2];
sx q[2];
rz(-1.8281724) q[2];
rz(-3.1016453) q[3];
sx q[3];
rz(-0.91244709) q[3];
sx q[3];
rz(-1.1279172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0082598) q[0];
sx q[0];
rz(-1.8513716) q[0];
sx q[0];
rz(2.7449352) q[0];
rz(2.2465514) q[1];
sx q[1];
rz(-1.7170693) q[1];
sx q[1];
rz(-0.76362124) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39786441) q[0];
sx q[0];
rz(-1.3258805) q[0];
sx q[0];
rz(-3.0586277) q[0];
x q[1];
rz(-1.5537797) q[2];
sx q[2];
rz(-1.6913171) q[2];
sx q[2];
rz(1.8486384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8473007) q[1];
sx q[1];
rz(-2.9025893) q[1];
sx q[1];
rz(1.4579178) q[1];
x q[2];
rz(1.607863) q[3];
sx q[3];
rz(-1.2801814) q[3];
sx q[3];
rz(0.99729482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.414173) q[2];
sx q[2];
rz(-1.5199993) q[2];
sx q[2];
rz(-2.011389) q[2];
rz(2.0027509) q[3];
sx q[3];
rz(-1.9500705) q[3];
sx q[3];
rz(1.1677008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30911699) q[0];
sx q[0];
rz(-2.583355) q[0];
sx q[0];
rz(0.48156893) q[0];
rz(2.8140977) q[1];
sx q[1];
rz(-1.7227252) q[1];
sx q[1];
rz(2.17735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010134546) q[0];
sx q[0];
rz(-1.9572963) q[0];
sx q[0];
rz(0.23614998) q[0];
rz(-pi) q[1];
rz(0.98924182) q[2];
sx q[2];
rz(-1.7316526) q[2];
sx q[2];
rz(-0.040078321) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89470071) q[1];
sx q[1];
rz(-0.27447005) q[1];
sx q[1];
rz(-1.1054705) q[1];
rz(-1.7751415) q[3];
sx q[3];
rz(-1.9540403) q[3];
sx q[3];
rz(2.1658773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1714736) q[2];
sx q[2];
rz(-2.0233266) q[2];
sx q[2];
rz(1.4946651) q[2];
rz(2.1286879) q[3];
sx q[3];
rz(-2.2154112) q[3];
sx q[3];
rz(-2.6265898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2038912) q[0];
sx q[0];
rz(-1.0045811) q[0];
sx q[0];
rz(1.0585693) q[0];
rz(0.92264289) q[1];
sx q[1];
rz(-2.1244996) q[1];
sx q[1];
rz(0.10291084) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71672463) q[0];
sx q[0];
rz(-2.1916336) q[0];
sx q[0];
rz(2.5607361) q[0];
rz(-1.2936959) q[2];
sx q[2];
rz(-0.99969802) q[2];
sx q[2];
rz(-1.4032316) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2051112) q[1];
sx q[1];
rz(-1.0432341) q[1];
sx q[1];
rz(2.9283306) q[1];
x q[2];
rz(-1.5779805) q[3];
sx q[3];
rz(-0.71643351) q[3];
sx q[3];
rz(-0.35182692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70204488) q[2];
sx q[2];
rz(-2.6921258) q[2];
sx q[2];
rz(-0.31583819) q[2];
rz(-0.99509197) q[3];
sx q[3];
rz(-1.8683878) q[3];
sx q[3];
rz(0.72527138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.134326) q[0];
sx q[0];
rz(-2.8747929) q[0];
sx q[0];
rz(0.63359487) q[0];
rz(-2.7569356) q[1];
sx q[1];
rz(-1.1793914) q[1];
sx q[1];
rz(0.42919174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4576526) q[0];
sx q[0];
rz(-2.5999617) q[0];
sx q[0];
rz(2.0787879) q[0];
rz(-pi) q[1];
x q[1];
rz(3.09893) q[2];
sx q[2];
rz(-0.68603078) q[2];
sx q[2];
rz(-1.3801284) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4218015) q[1];
sx q[1];
rz(-0.69703494) q[1];
sx q[1];
rz(-0.56682539) q[1];
x q[2];
rz(-2.2463838) q[3];
sx q[3];
rz(-0.67537245) q[3];
sx q[3];
rz(0.403756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0523494) q[2];
sx q[2];
rz(-2.3640335) q[2];
sx q[2];
rz(-2.8746129) q[2];
rz(-0.071768196) q[3];
sx q[3];
rz(-1.0177344) q[3];
sx q[3];
rz(1.4192386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778022) q[0];
sx q[0];
rz(-3.0402122) q[0];
sx q[0];
rz(-1.4005533) q[0];
rz(-1.9969223) q[1];
sx q[1];
rz(-2.0796227) q[1];
sx q[1];
rz(0.58715075) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82281821) q[0];
sx q[0];
rz(-0.94928926) q[0];
sx q[0];
rz(1.0448635) q[0];
rz(-1.6875504) q[2];
sx q[2];
rz(-1.1551305) q[2];
sx q[2];
rz(-2.1374201) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0441832) q[1];
sx q[1];
rz(-0.72260586) q[1];
sx q[1];
rz(0.34383066) q[1];
x q[2];
rz(2.1903992) q[3];
sx q[3];
rz(-1.8319329) q[3];
sx q[3];
rz(-1.2984683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2862386) q[2];
sx q[2];
rz(-1.9619532) q[2];
sx q[2];
rz(-0.41909763) q[2];
rz(1.6176443) q[3];
sx q[3];
rz(-2.2231299) q[3];
sx q[3];
rz(1.5036478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54290259) q[0];
sx q[0];
rz(-1.9200696) q[0];
sx q[0];
rz(-0.30938095) q[0];
rz(1.7912553) q[1];
sx q[1];
rz(-0.58809909) q[1];
sx q[1];
rz(-2.5877171) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2387008) q[0];
sx q[0];
rz(-2.6052778) q[0];
sx q[0];
rz(0.66030057) q[0];
rz(1.9543086) q[2];
sx q[2];
rz(-1.4330136) q[2];
sx q[2];
rz(-0.99064186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9123037) q[1];
sx q[1];
rz(-2.4472651) q[1];
sx q[1];
rz(-2.0635384) q[1];
rz(-2.7782544) q[3];
sx q[3];
rz(-0.73550341) q[3];
sx q[3];
rz(-0.21061646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89173335) q[2];
sx q[2];
rz(-1.6102108) q[2];
sx q[2];
rz(1.3472793) q[2];
rz(-0.76752082) q[3];
sx q[3];
rz(-1.9556421) q[3];
sx q[3];
rz(-2.2931113) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42221853) q[0];
sx q[0];
rz(-0.30000559) q[0];
sx q[0];
rz(-1.9268217) q[0];
rz(-1.0019373) q[1];
sx q[1];
rz(-1.5129713) q[1];
sx q[1];
rz(-0.92438662) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0262605) q[0];
sx q[0];
rz(-2.1450274) q[0];
sx q[0];
rz(2.6747245) q[0];
rz(-1.5255648) q[2];
sx q[2];
rz(-2.5965105) q[2];
sx q[2];
rz(-1.4807702) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6363057) q[1];
sx q[1];
rz(-1.5475379) q[1];
sx q[1];
rz(-0.22332286) q[1];
rz(-0.87189039) q[3];
sx q[3];
rz(-1.0950441) q[3];
sx q[3];
rz(-2.1995939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.34710106) q[2];
sx q[2];
rz(-0.93239409) q[2];
sx q[2];
rz(0.36716983) q[2];
rz(1.9764887) q[3];
sx q[3];
rz(-2.1737183) q[3];
sx q[3];
rz(2.2754106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871056) q[0];
sx q[0];
rz(-0.73910537) q[0];
sx q[0];
rz(0.11866971) q[0];
rz(0.26609303) q[1];
sx q[1];
rz(-0.91231822) q[1];
sx q[1];
rz(-1.4516861) q[1];
rz(-2.8546147) q[2];
sx q[2];
rz(-1.8331883) q[2];
sx q[2];
rz(2.2856648) q[2];
rz(-0.20826505) q[3];
sx q[3];
rz(-2.1622445) q[3];
sx q[3];
rz(1.6643659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
