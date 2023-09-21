OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(-2.8556813) q[0];
sx q[0];
rz(-0.51529348) q[0];
rz(1.3442858) q[1];
sx q[1];
rz(-2.9872515) q[1];
sx q[1];
rz(0.57758346) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56395036) q[0];
sx q[0];
rz(-0.9194153) q[0];
sx q[0];
rz(1.7786068) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0628113) q[2];
sx q[2];
rz(-1.0928109) q[2];
sx q[2];
rz(2.9302772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56596245) q[1];
sx q[1];
rz(-1.9735676) q[1];
sx q[1];
rz(-0.29213841) q[1];
rz(-2.3732784) q[3];
sx q[3];
rz(-0.81211219) q[3];
sx q[3];
rz(1.0031568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(2.4543767) q[2];
rz(2.1263188) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.8815536) q[0];
rz(1.0062224) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(-2.2959183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9189008) q[0];
sx q[0];
rz(-2.0881623) q[0];
sx q[0];
rz(-2.0738686) q[0];
rz(-pi) q[1];
rz(0.13375608) q[2];
sx q[2];
rz(-1.6998561) q[2];
sx q[2];
rz(-1.918902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7686292) q[1];
sx q[1];
rz(-1.0781204) q[1];
sx q[1];
rz(1.3586587) q[1];
x q[2];
rz(-0.91652292) q[3];
sx q[3];
rz(-2.407981) q[3];
sx q[3];
rz(-0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(-0.4804002) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.9539179) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(-2.3685266) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(-1.0864331) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7032996) q[0];
sx q[0];
rz(-2.6801531) q[0];
sx q[0];
rz(-1.5745387) q[0];
rz(-pi) q[1];
rz(1.7227206) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(2.8002847) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.8311027) q[1];
sx q[1];
rz(-1.5655787) q[1];
sx q[1];
rz(-0.28880854) q[1];
x q[2];
rz(-1.9589892) q[3];
sx q[3];
rz(-2.7079294) q[3];
sx q[3];
rz(2.7565623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(1.3712937) q[2];
rz(-2.7584372) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(-3.0932328) q[0];
rz(0.16391779) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(1.4455459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500344) q[0];
sx q[0];
rz(-1.7507179) q[0];
sx q[0];
rz(-0.66288373) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21638685) q[2];
sx q[2];
rz(-1.4828223) q[2];
sx q[2];
rz(0.83425922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7062263) q[1];
sx q[1];
rz(-1.5063783) q[1];
sx q[1];
rz(-2.9219342) q[1];
x q[2];
rz(-1.0501782) q[3];
sx q[3];
rz(-1.8225267) q[3];
sx q[3];
rz(-2.7057735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(1.6131489) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(-0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(-0.31907407) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(-2.2873926) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4695278) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(0.010849997) q[0];
rz(-2.1335667) q[2];
sx q[2];
rz(-2.3587583) q[2];
sx q[2];
rz(2.7340739) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90480587) q[1];
sx q[1];
rz(-0.85610897) q[1];
sx q[1];
rz(-0.76101117) q[1];
rz(-pi) q[2];
rz(2.3417926) q[3];
sx q[3];
rz(-0.721867) q[3];
sx q[3];
rz(-2.0197899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-0.577315) q[2];
rz(2.632085) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(-0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1381056) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(-1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(-0.12621005) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70996767) q[0];
sx q[0];
rz(-0.21490782) q[0];
sx q[0];
rz(0.14847319) q[0];
rz(0.38527617) q[2];
sx q[2];
rz(-2.3284973) q[2];
sx q[2];
rz(-0.77020459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.63657657) q[1];
sx q[1];
rz(-0.86537433) q[1];
sx q[1];
rz(-2.8935863) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0522271) q[3];
sx q[3];
rz(-2.7675981) q[3];
sx q[3];
rz(2.3911589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(-0.77511707) q[2];
rz(-2.3136247) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-3.0431842) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(0.55955204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3768809) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(-2.5401831) q[0];
rz(1.6227116) q[2];
sx q[2];
rz(-1.8318818) q[2];
sx q[2];
rz(-1.3494929) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.050484867) q[1];
sx q[1];
rz(-0.15639601) q[1];
sx q[1];
rz(0.97538235) q[1];
x q[2];
rz(-0.63013245) q[3];
sx q[3];
rz(-1.4117985) q[3];
sx q[3];
rz(2.7255448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(2.7977978) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(-0.47376537) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.375181) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-0.73076105) q[0];
rz(0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(2.2699845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9974737) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(3.0661422) q[0];
rz(-pi) q[1];
rz(2.7996054) q[2];
sx q[2];
rz(-0.42765289) q[2];
sx q[2];
rz(2.4004186) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.826088) q[1];
sx q[1];
rz(-1.9313889) q[1];
sx q[1];
rz(0.70198595) q[1];
x q[2];
rz(3.001508) q[3];
sx q[3];
rz(-2.1685371) q[3];
sx q[3];
rz(2.7261581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(1.139572) q[2];
rz(1.6428044) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(-2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9386439) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(-1.9198445) q[0];
rz(0.16601673) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(-1.6171914) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0696408) q[0];
sx q[0];
rz(-1.6531634) q[0];
sx q[0];
rz(3.1125493) q[0];
rz(-pi) q[1];
rz(-1.7622856) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(-2.7188403) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2251687) q[1];
sx q[1];
rz(-1.7290219) q[1];
sx q[1];
rz(1.9041063) q[1];
rz(-pi) q[2];
rz(0.41096656) q[3];
sx q[3];
rz(-2.4604359) q[3];
sx q[3];
rz(2.4364803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.021585492) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(0.35153708) q[2];
rz(-2.0848138) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(-2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64365023) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.836401) q[0];
rz(-0.38048831) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(0.25340733) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8587592) q[0];
sx q[0];
rz(-1.6738322) q[0];
sx q[0];
rz(-1.9673002) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9741653) q[2];
sx q[2];
rz(-1.2360459) q[2];
sx q[2];
rz(1.6850922) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0512143) q[1];
sx q[1];
rz(-1.9131294) q[1];
sx q[1];
rz(1.315209) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9990342) q[3];
sx q[3];
rz(-2.1381452) q[3];
sx q[3];
rz(-0.67728562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(-2.2425966) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(-1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9713365) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(-1.3118369) q[2];
sx q[2];
rz(-1.2987518) q[2];
sx q[2];
rz(-2.1546298) q[2];
rz(2.2181702) q[3];
sx q[3];
rz(-2.2150726) q[3];
sx q[3];
rz(1.1005145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];