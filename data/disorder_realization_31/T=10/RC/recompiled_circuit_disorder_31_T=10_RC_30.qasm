OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52656093) q[0];
sx q[0];
rz(-2.5685413) q[0];
sx q[0];
rz(-0.84258643) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(1.458118) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6095088) q[0];
sx q[0];
rz(-0.77021399) q[0];
sx q[0];
rz(2.8340333) q[0];
x q[1];
rz(-0.61383944) q[2];
sx q[2];
rz(-1.5547353) q[2];
sx q[2];
rz(1.2889372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1229413) q[1];
sx q[1];
rz(-0.39995799) q[1];
sx q[1];
rz(-2.8040228) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4487212) q[3];
sx q[3];
rz(-2.1616518) q[3];
sx q[3];
rz(-1.5199682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52790102) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(2.9620985) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(-0.32546145) q[0];
rz(-1.356396) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(1.1546086) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75630674) q[0];
sx q[0];
rz(-0.025408832) q[0];
sx q[0];
rz(2.4029762) q[0];
x q[1];
rz(-2.7484659) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(1.7413505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3719912) q[1];
sx q[1];
rz(-2.3725315) q[1];
sx q[1];
rz(-0.12188697) q[1];
rz(-pi) q[2];
rz(-1.3080018) q[3];
sx q[3];
rz(-1.7495219) q[3];
sx q[3];
rz(2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4521728) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(0.88341218) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.6261684) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-2.0498958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900345) q[0];
sx q[0];
rz(-1.0150195) q[0];
sx q[0];
rz(2.4843198) q[0];
rz(-pi) q[1];
rz(-1.0513564) q[2];
sx q[2];
rz(-2.4161077) q[2];
sx q[2];
rz(1.6348334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.162902) q[1];
sx q[1];
rz(-2.2802417) q[1];
sx q[1];
rz(2.6497926) q[1];
rz(-pi) q[2];
rz(-3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8213356) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(-0.88095218) q[2];
rz(-1.7679924) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110733) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(-2.7048892) q[0];
rz(-2.9084335) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(-2.8312347) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.045517) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(-1.4287352) q[0];
rz(-2.9878346) q[2];
sx q[2];
rz(-0.69090828) q[2];
sx q[2];
rz(-1.549987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0879285) q[1];
sx q[1];
rz(-1.9287319) q[1];
sx q[1];
rz(0.019836516) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6944828) q[3];
sx q[3];
rz(-1.5354772) q[3];
sx q[3];
rz(-0.24231054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(2.0641573) q[2];
rz(-0.056190101) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.5475387) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524277) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(1.5769618) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(-2.2713984) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10660431) q[0];
sx q[0];
rz(-2.7634794) q[0];
sx q[0];
rz(-0.58017054) q[0];
x q[1];
rz(2.444961) q[2];
sx q[2];
rz(-1.7156892) q[2];
sx q[2];
rz(-0.87755132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89228499) q[1];
sx q[1];
rz(-1.2586437) q[1];
sx q[1];
rz(0.73242964) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1327444) q[3];
sx q[3];
rz(-1.4268488) q[3];
sx q[3];
rz(2.2465003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(0.67374054) q[2];
rz(-0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(1.3195066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7917787) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(2.8836024) q[0];
rz(-0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(-1.4917096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8269577) q[0];
sx q[0];
rz(-1.1420113) q[0];
sx q[0];
rz(1.6786806) q[0];
x q[1];
rz(1.9906524) q[2];
sx q[2];
rz(-0.93517762) q[2];
sx q[2];
rz(-2.2395526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.7548435) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(0.90653231) q[1];
rz(-1.2061938) q[3];
sx q[3];
rz(-1.6657262) q[3];
sx q[3];
rz(1.3473542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(3.0498665) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82350746) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-0.033989865) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8695553) q[0];
sx q[0];
rz(-2.3453658) q[0];
sx q[0];
rz(1.3433334) q[0];
x q[1];
rz(-2.2074239) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(0.63461441) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.781144) q[1];
sx q[1];
rz(-2.2478659) q[1];
sx q[1];
rz(-3.0767246) q[1];
rz(1.4230698) q[3];
sx q[3];
rz(-0.96102321) q[3];
sx q[3];
rz(2.8019398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(-2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(2.7476655) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.6961018) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1720393) q[0];
sx q[0];
rz(-1.5656099) q[0];
sx q[0];
rz(1.1374377) q[0];
rz(2.1922853) q[2];
sx q[2];
rz(-1.0770123) q[2];
sx q[2];
rz(-2.2311503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.48965028) q[1];
sx q[1];
rz(-1.4133269) q[1];
sx q[1];
rz(0.55649346) q[1];
rz(-pi) q[2];
rz(-0.74650713) q[3];
sx q[3];
rz(-0.81614796) q[3];
sx q[3];
rz(-1.5796803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2945071) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(0.40714804) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(2.8318185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4315902) q[0];
sx q[0];
rz(-2.6012523) q[0];
sx q[0];
rz(-2.8938328) q[0];
x q[1];
rz(-1.3938815) q[2];
sx q[2];
rz(-1.5614911) q[2];
sx q[2];
rz(-2.6975346) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8621091) q[1];
sx q[1];
rz(-2.9276507) q[1];
sx q[1];
rz(-1.8474384) q[1];
rz(-0.59659776) q[3];
sx q[3];
rz(-1.784424) q[3];
sx q[3];
rz(-1.9006157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4391675) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.2058831) q[0];
rz(0.58569113) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(1.4996128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8382032) q[0];
sx q[0];
rz(-1.8614385) q[0];
sx q[0];
rz(-0.10586664) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6386547) q[2];
sx q[2];
rz(-2.3336764) q[2];
sx q[2];
rz(2.9615336) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2227576) q[1];
sx q[1];
rz(-1.3593874) q[1];
sx q[1];
rz(-1.097015) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5851192) q[3];
sx q[3];
rz(-2.9687772) q[3];
sx q[3];
rz(-0.9848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(2.1949027) q[2];
rz(-2.7729014) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
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
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(0.070925698) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(-2.6772066) q[2];
sx q[2];
rz(-1.1190363) q[2];
sx q[2];
rz(-2.6418532) q[2];
rz(2.0675038) q[3];
sx q[3];
rz(-0.91377331) q[3];
sx q[3];
rz(2.5627315) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
