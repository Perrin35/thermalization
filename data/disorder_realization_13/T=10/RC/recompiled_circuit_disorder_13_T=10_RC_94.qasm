OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(4.6306643) q[0];
sx q[0];
rz(10.319933) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(-1.4562343) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.668894) q[0];
sx q[0];
rz(-1.5510674) q[0];
sx q[0];
rz(1.701645) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37402447) q[2];
sx q[2];
rz(-1.0569388) q[2];
sx q[2];
rz(-0.49486578) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40741062) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(-0.68024866) q[1];
rz(2.8217444) q[3];
sx q[3];
rz(-1.5187129) q[3];
sx q[3];
rz(0.44714123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3866117) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(-0.1581986) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(-2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2125856) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(1.1520977) q[0];
rz(-1.2377897) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(2.6706085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5574871) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(1.6886061) q[0];
rz(-pi) q[1];
rz(-1.2330301) q[2];
sx q[2];
rz(-0.62440364) q[2];
sx q[2];
rz(-2.117346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4091332) q[1];
sx q[1];
rz(-1.0832548) q[1];
sx q[1];
rz(-3.0797466) q[1];
rz(-pi) q[2];
rz(3.0456411) q[3];
sx q[3];
rz(-2.4437685) q[3];
sx q[3];
rz(-2.5747091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.058078893) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(2.8857968) q[2];
rz(1.4852218) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(-2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8783022) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(0.27134744) q[0];
rz(-2.4052606) q[1];
sx q[1];
rz(-1.6059748) q[1];
sx q[1];
rz(-0.43930611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3210179) q[0];
sx q[0];
rz(-1.4903755) q[0];
sx q[0];
rz(0.029650173) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0053824) q[2];
sx q[2];
rz(-2.7779707) q[2];
sx q[2];
rz(-1.1656392) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9648793) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(1.8646851) q[1];
rz(-pi) q[2];
rz(0.4337173) q[3];
sx q[3];
rz(-2.0842413) q[3];
sx q[3];
rz(1.0182667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077483594) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0229605) q[0];
sx q[0];
rz(-1.4151787) q[0];
sx q[0];
rz(-0.96373425) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7165519) q[2];
sx q[2];
rz(-1.1894023) q[2];
sx q[2];
rz(2.8085453) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.306327) q[1];
sx q[1];
rz(-1.3439461) q[1];
sx q[1];
rz(-0.38362417) q[1];
x q[2];
rz(3.0365385) q[3];
sx q[3];
rz(-1.6451562) q[3];
sx q[3];
rz(0.33945938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0041634) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(-1.4271663) q[2];
rz(3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.78113294) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(0.57058913) q[0];
rz(-0.55496201) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(2.3983009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63440454) q[0];
sx q[0];
rz(-2.2162063) q[0];
sx q[0];
rz(-2.8651644) q[0];
rz(-pi) q[1];
rz(2.6402316) q[2];
sx q[2];
rz(-1.6622346) q[2];
sx q[2];
rz(2.3102592) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9651282) q[1];
sx q[1];
rz(-1.3511409) q[1];
sx q[1];
rz(-1.7996644) q[1];
rz(-pi) q[2];
rz(-2.3371313) q[3];
sx q[3];
rz(-1.3692229) q[3];
sx q[3];
rz(-1.3988914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-2.8732079) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(-0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.627581) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(1.0553029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5366093) q[0];
sx q[0];
rz(-0.6491001) q[0];
sx q[0];
rz(1.8761294) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91674532) q[2];
sx q[2];
rz(-1.142821) q[2];
sx q[2];
rz(-2.5452754) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9195337) q[1];
sx q[1];
rz(-2.2026081) q[1];
sx q[1];
rz(0.49948378) q[1];
x q[2];
rz(-2.8860693) q[3];
sx q[3];
rz(-0.75948411) q[3];
sx q[3];
rz(1.4042735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48173299) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(-0.64583889) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5053453) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(-2.3690467) q[0];
rz(1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(0.57377446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6360146) q[0];
sx q[0];
rz(-1.4609219) q[0];
sx q[0];
rz(1.7626324) q[0];
x q[1];
rz(-1.4622719) q[2];
sx q[2];
rz(-2.9402241) q[2];
sx q[2];
rz(0.078171922) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37807235) q[1];
sx q[1];
rz(-2.1280648) q[1];
sx q[1];
rz(-2.4850363) q[1];
rz(-2.9912234) q[3];
sx q[3];
rz(-1.339185) q[3];
sx q[3];
rz(-1.7600972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.27281) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(-1.1707206) q[0];
rz(2.6314578) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(-1.2957113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2670691) q[0];
sx q[0];
rz(-1.1133615) q[0];
sx q[0];
rz(1.9382856) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2175351) q[2];
sx q[2];
rz(-1.3152939) q[2];
sx q[2];
rz(3.0614292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.46591972) q[1];
sx q[1];
rz(-0.95341668) q[1];
sx q[1];
rz(0.27523756) q[1];
rz(2.0115764) q[3];
sx q[3];
rz(-1.545558) q[3];
sx q[3];
rz(-2.4079635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43508139) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(-2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4762964) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(0.67767674) q[0];
rz(0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(-0.83818865) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121376) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(1.8118993) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3269862) q[2];
sx q[2];
rz(-1.6199154) q[2];
sx q[2];
rz(2.5533822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14324489) q[1];
sx q[1];
rz(-1.3721826) q[1];
sx q[1];
rz(2.3193588) q[1];
rz(-pi) q[2];
rz(1.6543051) q[3];
sx q[3];
rz(-1.3132846) q[3];
sx q[3];
rz(-2.6058692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(-0.094749711) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(-0.23751968) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(0.231803) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7983539) q[0];
sx q[0];
rz(-0.51734561) q[0];
sx q[0];
rz(-1.377064) q[0];
rz(2.3583057) q[2];
sx q[2];
rz(-0.33595339) q[2];
sx q[2];
rz(0.17620262) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.17813645) q[1];
sx q[1];
rz(-1.3061211) q[1];
sx q[1];
rz(1.7174277) q[1];
rz(-pi) q[2];
rz(1.3123355) q[3];
sx q[3];
rz(-0.70445326) q[3];
sx q[3];
rz(2.7750912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0314177) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-2.0533662) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(1.5768455) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(1.8112524) q[2];
sx q[2];
rz(-1.8869055) q[2];
sx q[2];
rz(-1.4958924) q[2];
rz(-3.0392421) q[3];
sx q[3];
rz(-2.9526738) q[3];
sx q[3];
rz(-1.916193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
