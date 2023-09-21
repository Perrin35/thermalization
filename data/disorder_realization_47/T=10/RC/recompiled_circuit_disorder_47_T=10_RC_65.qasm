OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(-0.83431017) q[0];
sx q[0];
rz(3.0732529) q[0];
rz(-0.66032687) q[1];
sx q[1];
rz(-0.84815174) q[1];
sx q[1];
rz(0.037820427) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8605182) q[0];
sx q[0];
rz(-1.7639419) q[0];
sx q[0];
rz(-0.079695745) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7224563) q[2];
sx q[2];
rz(-2.4745387) q[2];
sx q[2];
rz(0.59282138) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5303858) q[1];
sx q[1];
rz(-1.5198277) q[1];
sx q[1];
rz(0.74688046) q[1];
rz(-pi) q[2];
rz(0.1511729) q[3];
sx q[3];
rz(-1.9133948) q[3];
sx q[3];
rz(1.3959988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7906856) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(1.8581871) q[2];
rz(3.0154199) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(3.0509907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.693817) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(-0.59666657) q[0];
rz(1.5555351) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(-1.7780875) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1648646) q[0];
sx q[0];
rz(-0.70978998) q[0];
sx q[0];
rz(-2.4564132) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5559276) q[2];
sx q[2];
rz(-1.613942) q[2];
sx q[2];
rz(-0.76669979) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1494257) q[1];
sx q[1];
rz(-0.70632315) q[1];
sx q[1];
rz(2.7041433) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76916839) q[3];
sx q[3];
rz(-2.2393919) q[3];
sx q[3];
rz(-1.234496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3229225) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(2.4070516) q[2];
rz(-2.5143886) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(-0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7221786) q[0];
sx q[0];
rz(-2.3086771) q[0];
sx q[0];
rz(2.869379) q[0];
rz(2.294337) q[1];
sx q[1];
rz(-1.8224742) q[1];
sx q[1];
rz(2.8289657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4011742) q[0];
sx q[0];
rz(-2.9156988) q[0];
sx q[0];
rz(0.24143879) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1001415) q[2];
sx q[2];
rz(-2.0690284) q[2];
sx q[2];
rz(0.91980308) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15309139) q[1];
sx q[1];
rz(-1.3107436) q[1];
sx q[1];
rz(3.0696763) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7496495) q[3];
sx q[3];
rz(-0.80229811) q[3];
sx q[3];
rz(2.9585569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0212038) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(0.81494251) q[2];
rz(0.96873823) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(-2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3880436) q[0];
sx q[0];
rz(-0.20064813) q[0];
sx q[0];
rz(0.62336212) q[0];
rz(0.81758824) q[1];
sx q[1];
rz(-0.31660429) q[1];
sx q[1];
rz(3.0923016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.277963) q[0];
sx q[0];
rz(-2.3290714) q[0];
sx q[0];
rz(0.066594007) q[0];
rz(-pi) q[1];
rz(-1.9360147) q[2];
sx q[2];
rz(-2.3962002) q[2];
sx q[2];
rz(-3.0129907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.019776736) q[1];
sx q[1];
rz(-0.41408086) q[1];
sx q[1];
rz(-2.2112234) q[1];
rz(-0.87818273) q[3];
sx q[3];
rz(-0.64681029) q[3];
sx q[3];
rz(1.9424903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4108882) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(2.1566186) q[2];
rz(2.2385521) q[3];
sx q[3];
rz(-1.9786381) q[3];
sx q[3];
rz(-2.8682017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34489283) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(3.0539736) q[0];
rz(-2.9852729) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(-2.2706251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0098457355) q[0];
sx q[0];
rz(-1.5994659) q[0];
sx q[0];
rz(0.28733758) q[0];
rz(-pi) q[1];
rz(-1.2147374) q[2];
sx q[2];
rz(-2.7241754) q[2];
sx q[2];
rz(2.2603214) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7333784) q[1];
sx q[1];
rz(-0.92792643) q[1];
sx q[1];
rz(-1.2182359) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17467588) q[3];
sx q[3];
rz(-0.81411618) q[3];
sx q[3];
rz(-1.0604309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.093420371) q[2];
sx q[2];
rz(-0.8478567) q[2];
sx q[2];
rz(-2.7887153) q[2];
rz(-0.86587632) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(0.29233366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6784994) q[0];
sx q[0];
rz(-1.1266288) q[0];
sx q[0];
rz(-0.10678664) q[0];
rz(-1.9550025) q[1];
sx q[1];
rz(-1.0266961) q[1];
sx q[1];
rz(-2.1170763) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2811919) q[0];
sx q[0];
rz(-1.2828579) q[0];
sx q[0];
rz(2.6698551) q[0];
rz(2.5249135) q[2];
sx q[2];
rz(-2.4384535) q[2];
sx q[2];
rz(2.6577735) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4752794) q[1];
sx q[1];
rz(-0.76994714) q[1];
sx q[1];
rz(-0.85789263) q[1];
x q[2];
rz(-2.6496274) q[3];
sx q[3];
rz(-1.2428189) q[3];
sx q[3];
rz(-2.5043951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(0.8708896) q[2];
rz(-1.332256) q[3];
sx q[3];
rz(-2.199316) q[3];
sx q[3];
rz(-1.4107305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9695327) q[0];
sx q[0];
rz(-1.7344069) q[0];
sx q[0];
rz(-0.55091888) q[0];
rz(2.6761966) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(-0.23682061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65695545) q[0];
sx q[0];
rz(-1.4716363) q[0];
sx q[0];
rz(0.037845503) q[0];
rz(-pi) q[1];
rz(0.15010712) q[2];
sx q[2];
rz(-1.4347715) q[2];
sx q[2];
rz(-1.0667691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6432861) q[1];
sx q[1];
rz(-0.28043881) q[1];
sx q[1];
rz(1.8449057) q[1];
rz(-pi) q[2];
rz(-2.3051466) q[3];
sx q[3];
rz(-1.8006547) q[3];
sx q[3];
rz(0.072007192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44234309) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(0.44011763) q[2];
rz(1.0951428) q[3];
sx q[3];
rz(-1.6710072) q[3];
sx q[3];
rz(1.346689) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8286164) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(-0.029504689) q[0];
rz(0.94738952) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(3.0227919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5598181) q[0];
sx q[0];
rz(-1.2604598) q[0];
sx q[0];
rz(-1.6161726) q[0];
x q[1];
rz(1.0256836) q[2];
sx q[2];
rz(-1.1720177) q[2];
sx q[2];
rz(1.1023956) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6606632) q[1];
sx q[1];
rz(-2.5813563) q[1];
sx q[1];
rz(-2.5653097) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0766255) q[3];
sx q[3];
rz(-1.4255187) q[3];
sx q[3];
rz(-2.6165917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3999195) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(1.5734394) q[2];
rz(-1.167477) q[3];
sx q[3];
rz(-1.7220595) q[3];
sx q[3];
rz(1.2742111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90299273) q[0];
sx q[0];
rz(-0.82505834) q[0];
sx q[0];
rz(2.9883244) q[0];
rz(2.080147) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(-1.9877888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35837367) q[0];
sx q[0];
rz(-1.8405387) q[0];
sx q[0];
rz(-2.7620035) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2527163) q[2];
sx q[2];
rz(-2.0098915) q[2];
sx q[2];
rz(0.097749226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3378539) q[1];
sx q[1];
rz(-1.6263464) q[1];
sx q[1];
rz(-1.3615723) q[1];
rz(-pi) q[2];
rz(0.12113573) q[3];
sx q[3];
rz(-2.2545358) q[3];
sx q[3];
rz(-0.63431206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8447421) q[2];
sx q[2];
rz(-1.1721609) q[2];
sx q[2];
rz(1.5273757) q[2];
rz(-2.1447003) q[3];
sx q[3];
rz(-1.4930054) q[3];
sx q[3];
rz(0.87866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24755724) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(0.19009185) q[0];
rz(-0.67063531) q[1];
sx q[1];
rz(-1.1963444) q[1];
sx q[1];
rz(1.6533096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3538441) q[0];
sx q[0];
rz(-1.4453381) q[0];
sx q[0];
rz(-1.6018484) q[0];
rz(-pi) q[1];
rz(-1.8411631) q[2];
sx q[2];
rz(-1.9989982) q[2];
sx q[2];
rz(2.2141475) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8390159) q[1];
sx q[1];
rz(-0.95514983) q[1];
sx q[1];
rz(-0.93977309) q[1];
rz(-pi) q[2];
rz(-0.90115746) q[3];
sx q[3];
rz(-2.3271051) q[3];
sx q[3];
rz(2.8709473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1378479) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(0.89912644) q[2];
rz(-0.51504618) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(-2.9639444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0638194) q[0];
sx q[0];
rz(-1.780092) q[0];
sx q[0];
rz(-0.5583981) q[0];
rz(-2.8521815) q[1];
sx q[1];
rz(-0.90129539) q[1];
sx q[1];
rz(1.7064066) q[1];
rz(-1.2895073) q[2];
sx q[2];
rz(-2.2608902) q[2];
sx q[2];
rz(2.8534941) q[2];
rz(1.9990986) q[3];
sx q[3];
rz(-2.8430568) q[3];
sx q[3];
rz(-1.6381016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];