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
rz(-1.8104115) q[0];
sx q[0];
rz(-1.6348732) q[0];
sx q[0];
rz(2.9659086) q[0];
rz(-2.076258) q[1];
sx q[1];
rz(-1.330436) q[1];
sx q[1];
rz(-2.4863844) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5436919) q[0];
sx q[0];
rz(-0.59698623) q[0];
sx q[0];
rz(-2.1587121) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51061012) q[2];
sx q[2];
rz(-3.0282927) q[2];
sx q[2];
rz(-2.4093546) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6083994) q[1];
sx q[1];
rz(-0.13020588) q[1];
sx q[1];
rz(0.59334468) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8484905) q[3];
sx q[3];
rz(-2.1610073) q[3];
sx q[3];
rz(2.9410998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4483999) q[2];
sx q[2];
rz(-0.11036631) q[2];
sx q[2];
rz(-2.8006862) q[2];
rz(-2.36813) q[3];
sx q[3];
rz(-1.0186467) q[3];
sx q[3];
rz(3.036518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.676749) q[0];
sx q[0];
rz(-2.8571547) q[0];
sx q[0];
rz(0.25962096) q[0];
rz(2.3530841) q[1];
sx q[1];
rz(-0.85616833) q[1];
sx q[1];
rz(-1.8964918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9282824) q[0];
sx q[0];
rz(-1.6716372) q[0];
sx q[0];
rz(-0.95243246) q[0];
rz(-1.670426) q[2];
sx q[2];
rz(-1.8848231) q[2];
sx q[2];
rz(-0.043836029) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40971995) q[1];
sx q[1];
rz(-1.221368) q[1];
sx q[1];
rz(-2.3278589) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29446843) q[3];
sx q[3];
rz(-0.52669062) q[3];
sx q[3];
rz(1.1931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9022687) q[2];
sx q[2];
rz(-0.50834346) q[2];
sx q[2];
rz(0.90947214) q[2];
rz(2.4746312) q[3];
sx q[3];
rz(-1.704155) q[3];
sx q[3];
rz(-0.10233574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-0.42720902) q[0];
sx q[0];
rz(-0.41218555) q[0];
sx q[0];
rz(1.7538196) q[0];
rz(2.1448403) q[1];
sx q[1];
rz(-0.69147384) q[1];
sx q[1];
rz(-0.48666418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5202058) q[0];
sx q[0];
rz(-0.4356111) q[0];
sx q[0];
rz(-1.2026879) q[0];
x q[1];
rz(2.4805807) q[2];
sx q[2];
rz(-0.95139388) q[2];
sx q[2];
rz(2.7255779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6909161) q[1];
sx q[1];
rz(-1.449739) q[1];
sx q[1];
rz(0.52074213) q[1];
rz(2.9938042) q[3];
sx q[3];
rz(-0.94919616) q[3];
sx q[3];
rz(0.10336598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77728689) q[2];
sx q[2];
rz(-1.8824258) q[2];
sx q[2];
rz(-2.6510009) q[2];
rz(-2.3250438) q[3];
sx q[3];
rz(-2.3994763) q[3];
sx q[3];
rz(-2.1013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8051324) q[0];
sx q[0];
rz(-1.8959624) q[0];
sx q[0];
rz(2.1601987) q[0];
rz(2.7384752) q[1];
sx q[1];
rz(-2.6336481) q[1];
sx q[1];
rz(0.53781646) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4889398) q[0];
sx q[0];
rz(-1.8171628) q[0];
sx q[0];
rz(-0.15423488) q[0];
rz(-pi) q[1];
rz(0.85015202) q[2];
sx q[2];
rz(-1.4701029) q[2];
sx q[2];
rz(-3.0792011) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0723426) q[1];
sx q[1];
rz(-1.1872655) q[1];
sx q[1];
rz(-2.093653) q[1];
rz(0.89610164) q[3];
sx q[3];
rz(-1.8180013) q[3];
sx q[3];
rz(-0.02442115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2526523) q[2];
sx q[2];
rz(-2.6344968) q[2];
sx q[2];
rz(-1.0520774) q[2];
rz(2.5893411) q[3];
sx q[3];
rz(-1.4058607) q[3];
sx q[3];
rz(-1.2205869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030516142) q[0];
sx q[0];
rz(-1.3330326) q[0];
sx q[0];
rz(3.1040177) q[0];
rz(-2.3887718) q[1];
sx q[1];
rz(-1.564874) q[1];
sx q[1];
rz(-0.084223025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1128453) q[0];
sx q[0];
rz(-1.0431093) q[0];
sx q[0];
rz(-0.28808388) q[0];
rz(2.7780121) q[2];
sx q[2];
rz(-1.2941133) q[2];
sx q[2];
rz(0.83438736) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9749637) q[1];
sx q[1];
rz(-1.5293603) q[1];
sx q[1];
rz(-1.1467481) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2920078) q[3];
sx q[3];
rz(-0.1670579) q[3];
sx q[3];
rz(-2.5045349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61421824) q[2];
sx q[2];
rz(-2.1671593) q[2];
sx q[2];
rz(-2.7091889) q[2];
rz(1.6480986) q[3];
sx q[3];
rz(-1.7688388) q[3];
sx q[3];
rz(0.73895946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768141) q[0];
sx q[0];
rz(-1.6141163) q[0];
sx q[0];
rz(-0.76171184) q[0];
rz(-1.0234045) q[1];
sx q[1];
rz(-0.49097148) q[1];
sx q[1];
rz(0.34058079) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6047594) q[0];
sx q[0];
rz(-2.3488461) q[0];
sx q[0];
rz(-1.6954846) q[0];
rz(-pi) q[1];
rz(-1.4144142) q[2];
sx q[2];
rz(-1.0151273) q[2];
sx q[2];
rz(1.2909691) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5300776) q[1];
sx q[1];
rz(-2.7002091) q[1];
sx q[1];
rz(-2.2468021) q[1];
rz(-pi) q[2];
rz(-0.75041308) q[3];
sx q[3];
rz(-1.6831796) q[3];
sx q[3];
rz(-3.0284798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9684101) q[2];
sx q[2];
rz(-1.1334271) q[2];
sx q[2];
rz(0.75562149) q[2];
rz(-0.38718811) q[3];
sx q[3];
rz(-1.3914934) q[3];
sx q[3];
rz(2.3473327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8887535) q[0];
sx q[0];
rz(-0.08789739) q[0];
sx q[0];
rz(1.9320236) q[0];
rz(1.6190489) q[1];
sx q[1];
rz(-2.3698898) q[1];
sx q[1];
rz(1.7542155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98859859) q[0];
sx q[0];
rz(-1.261088) q[0];
sx q[0];
rz(1.214908) q[0];
rz(-pi) q[1];
rz(0.31529324) q[2];
sx q[2];
rz(-1.0206501) q[2];
sx q[2];
rz(-1.8796876) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.55303326) q[1];
sx q[1];
rz(-0.84399763) q[1];
sx q[1];
rz(-0.56808205) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58297662) q[3];
sx q[3];
rz(-0.54123961) q[3];
sx q[3];
rz(-1.4547294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1986971) q[2];
sx q[2];
rz(-0.2936475) q[2];
sx q[2];
rz(-0.4100619) q[2];
rz(-0.53747082) q[3];
sx q[3];
rz(-2.0987857) q[3];
sx q[3];
rz(2.0736096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10802565) q[0];
sx q[0];
rz(-1.2355868) q[0];
sx q[0];
rz(-2.1349452) q[0];
rz(1.135723) q[1];
sx q[1];
rz(-2.6508811) q[1];
sx q[1];
rz(-0.27166414) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1170809) q[0];
sx q[0];
rz(-0.70444626) q[0];
sx q[0];
rz(-2.4570877) q[0];
x q[1];
rz(-0.67341398) q[2];
sx q[2];
rz(-1.0343164) q[2];
sx q[2];
rz(-1.3957889) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.063273009) q[1];
sx q[1];
rz(-0.98737475) q[1];
sx q[1];
rz(-0.88084817) q[1];
x q[2];
rz(2.397698) q[3];
sx q[3];
rz(-1.4049604) q[3];
sx q[3];
rz(-1.3506571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9378308) q[2];
sx q[2];
rz(-2.2088642) q[2];
sx q[2];
rz(1.1057828) q[2];
rz(-1.743099) q[3];
sx q[3];
rz(-2.1574557) q[3];
sx q[3];
rz(2.2891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6451013) q[0];
sx q[0];
rz(-2.190525) q[0];
sx q[0];
rz(2.7880461) q[0];
rz(-1.4191267) q[1];
sx q[1];
rz(-1.5357693) q[1];
sx q[1];
rz(3.0790192) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4434011) q[0];
sx q[0];
rz(-1.0834604) q[0];
sx q[0];
rz(-0.62701012) q[0];
rz(-pi) q[1];
x q[1];
rz(2.563971) q[2];
sx q[2];
rz(-2.2934539) q[2];
sx q[2];
rz(-1.3179145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8282916) q[1];
sx q[1];
rz(-0.57692674) q[1];
sx q[1];
rz(2.5531656) q[1];
x q[2];
rz(-1.4776286) q[3];
sx q[3];
rz(-0.9103295) q[3];
sx q[3];
rz(-2.4371393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8211557) q[2];
sx q[2];
rz(-1.7675567) q[2];
sx q[2];
rz(-2.721411) q[2];
rz(-0.34475103) q[3];
sx q[3];
rz(-0.77778608) q[3];
sx q[3];
rz(0.93728089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47278136) q[0];
sx q[0];
rz(-0.16952276) q[0];
sx q[0];
rz(-1.8602759) q[0];
rz(1.5877113) q[1];
sx q[1];
rz(-0.80931598) q[1];
sx q[1];
rz(0.70714998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65359914) q[0];
sx q[0];
rz(-0.63051134) q[0];
sx q[0];
rz(1.5911908) q[0];
rz(-0.94790876) q[2];
sx q[2];
rz(-2.4595408) q[2];
sx q[2];
rz(2.1932604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8445804) q[1];
sx q[1];
rz(-1.880901) q[1];
sx q[1];
rz(1.3363786) q[1];
x q[2];
rz(-0.94601241) q[3];
sx q[3];
rz(-2.1008228) q[3];
sx q[3];
rz(1.8998787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3215434) q[2];
sx q[2];
rz(-2.0115435) q[2];
sx q[2];
rz(0.496544) q[2];
rz(-1.8494362) q[3];
sx q[3];
rz(-0.29785952) q[3];
sx q[3];
rz(2.7636102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15688607) q[0];
sx q[0];
rz(-0.26837415) q[0];
sx q[0];
rz(-1.0347086) q[0];
rz(-0.59255076) q[1];
sx q[1];
rz(-0.68065803) q[1];
sx q[1];
rz(-1.5417644) q[1];
rz(1.7405657) q[2];
sx q[2];
rz(-0.23525722) q[2];
sx q[2];
rz(2.7973882) q[2];
rz(0.60382384) q[3];
sx q[3];
rz(-1.7573988) q[3];
sx q[3];
rz(0.16755541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
