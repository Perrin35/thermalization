OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(2.5685413) q[0];
sx q[0];
rz(11.723784) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(8.3254568) q[1];
sx q[1];
rz(7.96666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18544491) q[0];
sx q[0];
rz(-1.3584134) q[0];
sx q[0];
rz(-2.3953715) q[0];
rz(-pi) q[1];
rz(-0.61383944) q[2];
sx q[2];
rz(-1.5868574) q[2];
sx q[2];
rz(1.8526555) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2394489) q[1];
sx q[1];
rz(-1.7001171) q[1];
sx q[1];
rz(0.37954482) q[1];
rz(-pi) q[2];
rz(2.5472766) q[3];
sx q[3];
rz(-1.4694957) q[3];
sx q[3];
rz(3.1241824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(-0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(0.32546145) q[0];
rz(-1.356396) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(1.9869841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1240631) q[0];
sx q[0];
rz(-1.5895827) q[0];
sx q[0];
rz(-1.5879052) q[0];
rz(-pi) q[1];
rz(-2.7484659) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(1.7413505) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3719912) q[1];
sx q[1];
rz(-0.76906119) q[1];
sx q[1];
rz(0.12188697) q[1];
x q[2];
rz(-2.9566544) q[3];
sx q[3];
rz(-1.8293081) q[3];
sx q[3];
rz(2.1188494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(0.47131053) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31323355) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(-2.5405163) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(1.0916969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1611623) q[0];
sx q[0];
rz(-2.3083901) q[0];
sx q[0];
rz(-0.79361332) q[0];
x q[1];
rz(2.0902363) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(1.5067593) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8852383) q[1];
sx q[1];
rz(-1.2043722) q[1];
sx q[1];
rz(-2.3430235) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6767119) q[3];
sx q[3];
rz(-2.998623) q[3];
sx q[3];
rz(0.79899341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.320257) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(-0.88095218) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110733) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(0.23315915) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(2.8312347) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096075637) q[0];
sx q[0];
rz(-1.9487511) q[0];
sx q[0];
rz(1.7128574) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15375806) q[2];
sx q[2];
rz(-0.69090828) q[2];
sx q[2];
rz(-1.549987) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.031361) q[1];
sx q[1];
rz(-0.35846113) q[1];
sx q[1];
rz(-1.5178174) q[1];
rz(1.4471099) q[3];
sx q[3];
rz(-1.5354772) q[3];
sx q[3];
rz(2.8992821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(1.0774353) q[2];
rz(3.0854026) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.189165) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(-2.8919343) q[0];
rz(1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(-0.87019428) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2244959) q[0];
sx q[0];
rz(-1.774569) q[0];
sx q[0];
rz(-0.32075551) q[0];
x q[1];
rz(-2.444961) q[2];
sx q[2];
rz(-1.7156892) q[2];
sx q[2];
rz(-2.2640413) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2493077) q[1];
sx q[1];
rz(-1.882949) q[1];
sx q[1];
rz(-0.73242964) q[1];
rz(-1.3051885) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(0.8997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(-1.3195066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(2.8836024) q[0];
rz(-0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(1.649883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0598037) q[0];
sx q[0];
rz(-2.7002618) q[0];
sx q[0];
rz(-2.9102737) q[0];
rz(-pi) q[1];
rz(2.4620373) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(-0.92781767) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7266453) q[1];
sx q[1];
rz(-2.2989095) q[1];
sx q[1];
rz(0.77264087) q[1];
x q[2];
rz(-1.8317354) q[3];
sx q[3];
rz(-2.76537) q[3];
sx q[3];
rz(0.46686831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(-2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(-0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82350746) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(-0.86589083) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-0.033989865) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720374) q[0];
sx q[0];
rz(-2.3453658) q[0];
sx q[0];
rz(1.3433334) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6007973) q[2];
sx q[2];
rz(-2.135709) q[2];
sx q[2];
rz(0.62513798) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36044869) q[1];
sx q[1];
rz(-2.2478659) q[1];
sx q[1];
rz(3.0767246) q[1];
rz(-pi) q[2];
rz(-0.20766081) q[3];
sx q[3];
rz(-0.62519473) q[3];
sx q[3];
rz(0.085426424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8975163) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.6961018) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5515585) q[0];
sx q[0];
rz(-2.708205) q[0];
sx q[0];
rz(-1.5831468) q[0];
rz(-pi) q[1];
rz(-0.94930737) q[2];
sx q[2];
rz(-2.0645803) q[2];
sx q[2];
rz(-0.91044237) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3280914) q[1];
sx q[1];
rz(-2.5655167) q[1];
sx q[1];
rz(0.2920132) q[1];
rz(-0.74650713) q[3];
sx q[3];
rz(-0.81614796) q[3];
sx q[3];
rz(-1.5796803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2945071) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-2.7344446) q[2];
rz(-1.5173222) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.80609926) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(-0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-0.30977419) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99684925) q[0];
sx q[0];
rz(-1.0486756) q[0];
sx q[0];
rz(-1.4247308) q[0];
rz(-1.7477112) q[2];
sx q[2];
rz(-1.5801016) q[2];
sx q[2];
rz(-2.6975346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.138315) q[1];
sx q[1];
rz(-1.365108) q[1];
sx q[1];
rz(-3.0823207) q[1];
x q[2];
rz(-1.8272607) q[3];
sx q[3];
rz(-0.98955742) q[3];
sx q[3];
rz(2.6687711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(1.6419798) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(3.0707804) q[2];
sx q[2];
rz(-2.376308) q[2];
sx q[2];
rz(-2.8634957) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.91883509) q[1];
sx q[1];
rz(-1.3593874) q[1];
sx q[1];
rz(2.0445776) q[1];
rz(-pi) q[2];
rz(-1.5851192) q[3];
sx q[3];
rz(-0.17281547) q[3];
sx q[3];
rz(2.1567791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-2.1949027) q[2];
rz(0.36869129) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(-2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(-0.82540853) q[2];
sx q[2];
rz(-2.5054629) q[2];
sx q[2];
rz(2.7873743) q[2];
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