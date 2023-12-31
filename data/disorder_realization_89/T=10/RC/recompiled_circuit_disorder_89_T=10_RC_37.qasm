OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(2.6309738) q[1];
sx q[1];
rz(9.0030158) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3573787) q[0];
sx q[0];
rz(-1.5970018) q[0];
sx q[0];
rz(-1.635301) q[0];
rz(-1.5316891) q[2];
sx q[2];
rz(-2.2570059) q[2];
sx q[2];
rz(-1.4456911) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64695839) q[1];
sx q[1];
rz(-1.1263493) q[1];
sx q[1];
rz(1.7723945) q[1];
rz(-2.5476417) q[3];
sx q[3];
rz(-1.6868601) q[3];
sx q[3];
rz(2.5719197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(1.0100693) q[2];
rz(-1.646237) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0533957) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(0.37047085) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7709748) q[0];
sx q[0];
rz(-0.3807225) q[0];
sx q[0];
rz(2.5230663) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6071762) q[2];
sx q[2];
rz(-2.538531) q[2];
sx q[2];
rz(0.39819983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49711455) q[1];
sx q[1];
rz(-1.2116355) q[1];
sx q[1];
rz(0.70980806) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8584077) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(-2.9420497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7039965) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(0.55830467) q[2];
rz(-2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(-2.5487652) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-2.9779789) q[0];
rz(2.4257461) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4931902) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(-0.52669749) q[0];
rz(-pi) q[1];
rz(-0.30544124) q[2];
sx q[2];
rz(-0.88103308) q[2];
sx q[2];
rz(-1.5058215) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.19792412) q[1];
sx q[1];
rz(-2.1197882) q[1];
sx q[1];
rz(2.3180069) q[1];
rz(-0.22294873) q[3];
sx q[3];
rz(-1.795459) q[3];
sx q[3];
rz(-0.49062452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(2.1219357) q[2];
rz(-2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(-2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402128) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2903039) q[0];
sx q[0];
rz(-1.5820869) q[0];
sx q[0];
rz(-1.4615061) q[0];
rz(-pi) q[1];
rz(-1.9539815) q[2];
sx q[2];
rz(-1.2387347) q[2];
sx q[2];
rz(3.1021169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1000881) q[1];
sx q[1];
rz(-1.3106292) q[1];
sx q[1];
rz(-0.86160223) q[1];
rz(-pi) q[2];
rz(2.652076) q[3];
sx q[3];
rz(-2.094141) q[3];
sx q[3];
rz(0.19260064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.443976) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(-0.90732968) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089371) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.6089815) q[0];
rz(0.73348796) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.3132494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8719296) q[0];
sx q[0];
rz(-2.4047244) q[0];
sx q[0];
rz(-2.065573) q[0];
x q[1];
rz(-2.1336742) q[2];
sx q[2];
rz(-3.0745227) q[2];
sx q[2];
rz(1.1941393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.65391958) q[1];
sx q[1];
rz(-2.6012585) q[1];
sx q[1];
rz(0.78507702) q[1];
rz(2.5562708) q[3];
sx q[3];
rz(-2.1748741) q[3];
sx q[3];
rz(-2.4620591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1398754) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(2.4772947) q[2];
rz(-0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(-0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649081) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(-0.054811906) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(0.48318133) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7532363) q[0];
sx q[0];
rz(-1.6325132) q[0];
sx q[0];
rz(-1.389099) q[0];
x q[1];
rz(-1.1081084) q[2];
sx q[2];
rz(-1.5414642) q[2];
sx q[2];
rz(-2.6372452) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9890081) q[1];
sx q[1];
rz(-2.2629988) q[1];
sx q[1];
rz(-0.89905507) q[1];
rz(-pi) q[2];
rz(0.11404927) q[3];
sx q[3];
rz(-1.2403135) q[3];
sx q[3];
rz(1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(0.085263578) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(0.95463395) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(-0.25209299) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4567586) q[0];
sx q[0];
rz(-1.7157468) q[0];
sx q[0];
rz(-2.7177939) q[0];
x q[1];
rz(-0.83129518) q[2];
sx q[2];
rz(-2.202946) q[2];
sx q[2];
rz(2.0668154) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5735057) q[1];
sx q[1];
rz(-1.9316257) q[1];
sx q[1];
rz(2.2682857) q[1];
x q[2];
rz(1.1105359) q[3];
sx q[3];
rz(-1.9122951) q[3];
sx q[3];
rz(-0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38348848) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(-2.2371116) q[2];
rz(-2.1379743) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43338183) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(1.7277539) q[0];
rz(2.4811603) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(-1.7699014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577514) q[0];
sx q[0];
rz(-0.62566602) q[0];
sx q[0];
rz(-1.7218504) q[0];
rz(-0.76099446) q[2];
sx q[2];
rz(-1.4174263) q[2];
sx q[2];
rz(0.93190565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.614537) q[1];
sx q[1];
rz(-1.0053047) q[1];
sx q[1];
rz(1.2177699) q[1];
rz(-pi) q[2];
rz(2.780704) q[3];
sx q[3];
rz(-1.4286255) q[3];
sx q[3];
rz(-1.7796381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(-0.39789847) q[2];
rz(-0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.21815498) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(2.9072705) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(2.871002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2345679) q[0];
sx q[0];
rz(-1.7922635) q[0];
sx q[0];
rz(1.7937167) q[0];
x q[1];
rz(-1.2241237) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(-0.62435645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8447664) q[1];
sx q[1];
rz(-1.4335732) q[1];
sx q[1];
rz(1.574135) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3769334) q[3];
sx q[3];
rz(-0.48741515) q[3];
sx q[3];
rz(2.2262115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33971912) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(-1.489893) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40437317) q[0];
sx q[0];
rz(-1.568856) q[0];
sx q[0];
rz(0.15429601) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-2.0097201) q[1];
sx q[1];
rz(-2.399209) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267747) q[0];
sx q[0];
rz(-0.97133884) q[0];
sx q[0];
rz(-0.74032289) q[0];
rz(-pi) q[1];
rz(-2.7102091) q[2];
sx q[2];
rz(-2.1200392) q[2];
sx q[2];
rz(-2.0768349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4780477) q[1];
sx q[1];
rz(-1.5596584) q[1];
sx q[1];
rz(-2.7958109) q[1];
x q[2];
rz(2.9418482) q[3];
sx q[3];
rz(-2.9962857) q[3];
sx q[3];
rz(0.47606459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8538889) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(1.5596191) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29006526) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(-1.3416946) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(2.3201597) q[2];
sx q[2];
rz(-1.3315622) q[2];
sx q[2];
rz(-0.17377725) q[2];
rz(0.7704173) q[3];
sx q[3];
rz(-0.19376783) q[3];
sx q[3];
rz(2.7289058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
