OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(-2.196329) q[0];
sx q[0];
rz(0.52559108) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5358937) q[0];
sx q[0];
rz(-0.68471691) q[0];
sx q[0];
rz(-0.84795714) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9154539) q[2];
sx q[2];
rz(-1.7625426) q[2];
sx q[2];
rz(2.6297671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.875784) q[1];
sx q[1];
rz(-1.335653) q[1];
sx q[1];
rz(-1.9782515) q[1];
x q[2];
rz(-2.9526887) q[3];
sx q[3];
rz(-1.8763262) q[3];
sx q[3];
rz(0.73959914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(3.0453483) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(-0.76517117) q[0];
rz(1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5040341) q[0];
sx q[0];
rz(-0.090888977) q[0];
sx q[0];
rz(-2.3239517) q[0];
rz(-pi) q[1];
rz(-2.0775954) q[2];
sx q[2];
rz(-0.71841824) q[2];
sx q[2];
rz(1.6990627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6389097) q[1];
sx q[1];
rz(-2.0109634) q[1];
sx q[1];
rz(1.6619976) q[1];
rz(2.6547673) q[3];
sx q[3];
rz(-1.4273705) q[3];
sx q[3];
rz(-1.6008582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(0.32901397) q[2];
rz(0.66550231) q[3];
sx q[3];
rz(-2.9232959) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(-1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(0.51868784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1919353) q[0];
sx q[0];
rz(-1.9266832) q[0];
sx q[0];
rz(-2.3200672) q[0];
rz(-2.8027595) q[2];
sx q[2];
rz(-0.28140861) q[2];
sx q[2];
rz(-1.1178521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7171214) q[1];
sx q[1];
rz(-2.3470504) q[1];
sx q[1];
rz(2.8984927) q[1];
x q[2];
rz(1.5716001) q[3];
sx q[3];
rz(-1.5161848) q[3];
sx q[3];
rz(-2.2664203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(0.60244256) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(-0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(0.17424507) q[0];
rz(-0.53025591) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(0.51309103) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.355277) q[0];
sx q[0];
rz(-2.516541) q[0];
sx q[0];
rz(-0.13537188) q[0];
rz(-pi) q[1];
rz(-2.1999173) q[2];
sx q[2];
rz(-2.2721014) q[2];
sx q[2];
rz(0.47052449) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7375813) q[1];
sx q[1];
rz(-2.0239132) q[1];
sx q[1];
rz(-1.416942) q[1];
rz(-pi) q[2];
rz(0.030839132) q[3];
sx q[3];
rz(-1.3645002) q[3];
sx q[3];
rz(-1.7145969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47485581) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(-2.722548) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(0.45927799) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(2.4647734) q[0];
rz(-0.49304402) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(-0.61606032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13000935) q[0];
sx q[0];
rz(-1.5349712) q[0];
sx q[0];
rz(-1.623276) q[0];
rz(-pi) q[1];
rz(1.6904171) q[2];
sx q[2];
rz(-1.9336485) q[2];
sx q[2];
rz(-0.56064831) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38356009) q[1];
sx q[1];
rz(-2.5639503) q[1];
sx q[1];
rz(3.1314965) q[1];
x q[2];
rz(-2.2474399) q[3];
sx q[3];
rz(-1.3291306) q[3];
sx q[3];
rz(1.579293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(2.8862254) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(2.6155112) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(0.88476673) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677744) q[0];
sx q[0];
rz(-0.38823715) q[0];
sx q[0];
rz(-1.3740747) q[0];
rz(-pi) q[1];
rz(3.0753291) q[2];
sx q[2];
rz(-2.931086) q[2];
sx q[2];
rz(1.2103684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.51984519) q[1];
sx q[1];
rz(-1.8583082) q[1];
sx q[1];
rz(-1.9985755) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1690508) q[3];
sx q[3];
rz(-2.5707977) q[3];
sx q[3];
rz(-1.6789544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3970967) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-2.6122724) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41806528) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-3.080522) q[0];
rz(-0.04018499) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(-2.4087002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7948579) q[0];
sx q[0];
rz(-1.1746527) q[0];
sx q[0];
rz(-2.5651155) q[0];
x q[1];
rz(-1.0211208) q[2];
sx q[2];
rz(-0.78133821) q[2];
sx q[2];
rz(-2.6598425) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9961173) q[1];
sx q[1];
rz(-2.5158094) q[1];
sx q[1];
rz(-1.5596418) q[1];
rz(-pi) q[2];
rz(1.93768) q[3];
sx q[3];
rz(-0.11493472) q[3];
sx q[3];
rz(-0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0733033) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(-2.627009) q[2];
rz(1.9040646) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056203689) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-0.7094267) q[0];
rz(1.6363232) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(-0.27871305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3942791) q[0];
sx q[0];
rz(-0.25082591) q[0];
sx q[0];
rz(0.5745116) q[0];
rz(-pi) q[1];
rz(0.3785554) q[2];
sx q[2];
rz(-0.59213973) q[2];
sx q[2];
rz(2.2373667) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.41234327) q[1];
sx q[1];
rz(-1.3530429) q[1];
sx q[1];
rz(-1.2530243) q[1];
x q[2];
rz(3.0635733) q[3];
sx q[3];
rz(-1.3481513) q[3];
sx q[3];
rz(-1.5732461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(-3.0855132) q[2];
rz(2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(0.40518951) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99496019) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(0.96889281) q[0];
rz(2.6682207) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.1425346) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5794967) q[0];
sx q[0];
rz(-1.2764494) q[0];
sx q[0];
rz(-3.0466945) q[0];
rz(-3.1387024) q[2];
sx q[2];
rz(-0.70334607) q[2];
sx q[2];
rz(0.10290111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6856319) q[1];
sx q[1];
rz(-0.5545485) q[1];
sx q[1];
rz(-2.901652) q[1];
x q[2];
rz(0.90060602) q[3];
sx q[3];
rz(-1.4108676) q[3];
sx q[3];
rz(-0.98774324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.896495) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(2.1208105) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616515) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(-0.7014057) q[0];
rz(-2.2258863) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(1.2385626) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1322051) q[0];
sx q[0];
rz(-1.3552249) q[0];
sx q[0];
rz(-2.0995887) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.086834) q[2];
sx q[2];
rz(-2.3878532) q[2];
sx q[2];
rz(-2.4278305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.065579942) q[1];
sx q[1];
rz(-2.2231243) q[1];
sx q[1];
rz(0.14821649) q[1];
x q[2];
rz(0.11390399) q[3];
sx q[3];
rz(-0.98131991) q[3];
sx q[3];
rz(1.9437499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4548268) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(-2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713521) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(-3.1162221) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(-2.8364137) q[2];
sx q[2];
rz(-1.9532433) q[2];
sx q[2];
rz(0.80079186) q[2];
rz(-1.3958037) q[3];
sx q[3];
rz(-0.77944118) q[3];
sx q[3];
rz(-0.13164095) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
