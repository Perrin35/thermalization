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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9561477) q[0];
sx q[0];
rz(-1.7831793) q[0];
sx q[0];
rz(0.74622112) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5904434) q[2];
sx q[2];
rz(-0.95704776) q[2];
sx q[2];
rz(2.8710499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2394489) q[1];
sx q[1];
rz(-1.7001171) q[1];
sx q[1];
rz(-0.37954482) q[1];
rz(0.17957844) q[3];
sx q[3];
rz(-0.60185963) q[3];
sx q[3];
rz(1.7367401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52790102) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(-1.9159296) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801385) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(-2.8161312) q[0];
rz(1.7851967) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(1.9869841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5886473) q[0];
sx q[0];
rz(-1.5879022) q[0];
sx q[0];
rz(0.018789142) q[0];
rz(-pi) q[1];
rz(0.39312675) q[2];
sx q[2];
rz(-0.98193411) q[2];
sx q[2];
rz(-1.7413505) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6007538) q[1];
sx q[1];
rz(-0.80889091) q[1];
sx q[1];
rz(1.6879338) q[1];
x q[2];
rz(2.1784337) q[3];
sx q[3];
rz(-2.8249486) q[3];
sx q[3];
rz(0.39099993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4521728) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(-0.47131053) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.5154243) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(1.0916969) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9804304) q[0];
sx q[0];
rz(-2.3083901) q[0];
sx q[0];
rz(-0.79361332) q[0];
x q[1];
rz(-2.7269084) q[2];
sx q[2];
rz(-2.1846002) q[2];
sx q[2];
rz(-2.2874122) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.162902) q[1];
sx q[1];
rz(-2.2802417) q[1];
sx q[1];
rz(2.6497926) q[1];
x q[2];
rz(1.5063498) q[3];
sx q[3];
rz(-1.4430874) q[3];
sx q[3];
rz(-2.8116022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.320257) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(-2.7048892) q[0];
rz(-0.23315915) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(-0.31035796) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4219907) q[0];
sx q[0];
rz(-1.4388226) q[0];
sx q[0];
rz(-2.7601526) q[0];
rz(-pi) q[1];
rz(-0.68508673) q[2];
sx q[2];
rz(-1.473046) q[2];
sx q[2];
rz(0.098066559) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0536641) q[1];
sx q[1];
rz(-1.2128608) q[1];
sx q[1];
rz(3.1217561) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.849732) q[3];
sx q[3];
rz(-3.012987) q[3];
sx q[3];
rz(2.0898553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(-1.0774353) q[2];
rz(-3.0854026) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524277) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(0.24965832) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(0.87019428) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4207626) q[0];
sx q[0];
rz(-1.8846858) q[0];
sx q[0];
rz(-1.785196) q[0];
x q[1];
rz(-1.382803) q[2];
sx q[2];
rz(-0.88289875) q[2];
sx q[2];
rz(-0.57304136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.89228499) q[1];
sx q[1];
rz(-1.882949) q[1];
sx q[1];
rz(0.73242964) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16964511) q[3];
sx q[3];
rz(-2.1262453) q[3];
sx q[3];
rz(2.5559705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(0.30361787) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(1.3195066) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-0.2579903) q[0];
rz(2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.649883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0598037) q[0];
sx q[0];
rz(-2.7002618) q[0];
sx q[0];
rz(-0.23131891) q[0];
x q[1];
rz(2.4620373) q[2];
sx q[2];
rz(-1.9050042) q[2];
sx q[2];
rz(0.92781767) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7266453) q[1];
sx q[1];
rz(-0.84268314) q[1];
sx q[1];
rz(-2.3689518) q[1];
rz(-pi) q[2];
rz(-1.2061938) q[3];
sx q[3];
rz(-1.4758665) q[3];
sx q[3];
rz(1.7942384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(-2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.82350746) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(2.7303625) q[0];
rz(0.86589083) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(3.1076028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8695553) q[0];
sx q[0];
rz(-0.79622686) q[0];
sx q[0];
rz(-1.3433334) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54079536) q[2];
sx q[2];
rz(-2.135709) q[2];
sx q[2];
rz(-2.5164547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.781144) q[1];
sx q[1];
rz(-0.89372674) q[1];
sx q[1];
rz(-3.0767246) q[1];
rz(-pi) q[2];
rz(0.61492413) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(-1.8254335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(0.87654385) q[2];
rz(-2.792568) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(-0.36755964) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.4454909) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5515585) q[0];
sx q[0];
rz(-0.43338767) q[0];
sx q[0];
rz(-1.5831468) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82463512) q[2];
sx q[2];
rz(-2.3687009) q[2];
sx q[2];
rz(0.075721272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3280914) q[1];
sx q[1];
rz(-0.5760759) q[1];
sx q[1];
rz(0.2920132) q[1];
rz(-2.3950855) q[3];
sx q[3];
rz(-0.81614796) q[3];
sx q[3];
rz(1.5796803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8470856) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(0.40714804) q[2];
rz(1.5173222) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(2.8318185) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944086) q[0];
sx q[0];
rz(-1.697288) q[0];
sx q[0];
rz(-2.6148318) q[0];
x q[1];
rz(1.7477112) q[2];
sx q[2];
rz(-1.5801016) q[2];
sx q[2];
rz(-0.44405802) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8621091) q[1];
sx q[1];
rz(-2.9276507) q[1];
sx q[1];
rz(-1.2941542) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8272607) q[3];
sx q[3];
rz(-2.1520352) q[3];
sx q[3];
rz(-0.47282156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70242515) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(-0.65565482) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050215125) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.9357095) q[0];
rz(2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.6419798) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382032) q[0];
sx q[0];
rz(-1.2801542) q[0];
sx q[0];
rz(3.035726) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.070812289) q[2];
sx q[2];
rz(-0.76528463) q[2];
sx q[2];
rz(-0.27809696) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8780898) q[1];
sx q[1];
rz(-2.6260758) q[1];
sx q[1];
rz(1.131119) q[1];
x q[2];
rz(1.7435944) q[3];
sx q[3];
rz(-1.5683335) q[3];
sx q[3];
rz(2.5415004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.88400921) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(0.94669) q[2];
rz(2.7729014) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(3.070667) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(-2.6772066) q[2];
sx q[2];
rz(-1.1190363) q[2];
sx q[2];
rz(-2.6418532) q[2];
rz(-1.0740888) q[3];
sx q[3];
rz(-0.91377331) q[3];
sx q[3];
rz(2.5627315) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
