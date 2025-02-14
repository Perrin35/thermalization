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
rz(-0.35204044) q[0];
sx q[0];
rz(-0.8249324) q[0];
sx q[0];
rz(-2.6068249) q[0];
rz(0.98400247) q[1];
sx q[1];
rz(-2.6360631) q[1];
sx q[1];
rz(-1.8796138) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0875435) q[0];
sx q[0];
rz(-0.66960483) q[0];
sx q[0];
rz(-0.43144193) q[0];
x q[1];
rz(2.6017351) q[2];
sx q[2];
rz(-2.1261907) q[2];
sx q[2];
rz(2.843631) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83916503) q[1];
sx q[1];
rz(-1.5189956) q[1];
sx q[1];
rz(2.5940345) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46686904) q[3];
sx q[3];
rz(-0.43064603) q[3];
sx q[3];
rz(2.1809585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8925605) q[2];
sx q[2];
rz(-2.6292215) q[2];
sx q[2];
rz(-1.6585635) q[2];
rz(2.856971) q[3];
sx q[3];
rz(-2.5183545) q[3];
sx q[3];
rz(-2.0774138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8677218) q[0];
sx q[0];
rz(-0.72672788) q[0];
sx q[0];
rz(0.04100767) q[0];
rz(1.2076591) q[1];
sx q[1];
rz(-0.227808) q[1];
sx q[1];
rz(1.4606732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6448633) q[0];
sx q[0];
rz(-1.71642) q[0];
sx q[0];
rz(-1.54099) q[0];
rz(0.76717214) q[2];
sx q[2];
rz(-1.1572654) q[2];
sx q[2];
rz(-0.77595475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8833958) q[1];
sx q[1];
rz(-0.25942311) q[1];
sx q[1];
rz(-1.4197465) q[1];
x q[2];
rz(2.9462141) q[3];
sx q[3];
rz(-2.2535705) q[3];
sx q[3];
rz(1.186779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4760806) q[2];
sx q[2];
rz(-1.6802695) q[2];
sx q[2];
rz(1.3903138) q[2];
rz(0.0811854) q[3];
sx q[3];
rz(-0.66892162) q[3];
sx q[3];
rz(-1.4359052) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7471033) q[0];
sx q[0];
rz(-3.1293226) q[0];
sx q[0];
rz(-2.5315206) q[0];
rz(-1.3129781) q[1];
sx q[1];
rz(-2.4459631) q[1];
sx q[1];
rz(-0.55061603) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074215502) q[0];
sx q[0];
rz(-0.91892892) q[0];
sx q[0];
rz(-1.1135191) q[0];
rz(-pi) q[1];
rz(2.4464452) q[2];
sx q[2];
rz(-0.41054976) q[2];
sx q[2];
rz(-1.6352959) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0455835) q[1];
sx q[1];
rz(-2.7852045) q[1];
sx q[1];
rz(2.136904) q[1];
x q[2];
rz(-2.315963) q[3];
sx q[3];
rz(-1.6761002) q[3];
sx q[3];
rz(1.8838175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6637471) q[2];
sx q[2];
rz(-0.17243324) q[2];
sx q[2];
rz(-2.0752068) q[2];
rz(1.9829228) q[3];
sx q[3];
rz(-1.7585157) q[3];
sx q[3];
rz(-3.0995479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31024194) q[0];
sx q[0];
rz(-2.5284335) q[0];
sx q[0];
rz(-2.2245275) q[0];
rz(0.65545583) q[1];
sx q[1];
rz(-0.90924811) q[1];
sx q[1];
rz(-2.7520032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661267) q[0];
sx q[0];
rz(-1.4891169) q[0];
sx q[0];
rz(1.7577175) q[0];
rz(-pi) q[1];
rz(2.2245313) q[2];
sx q[2];
rz(-2.3417763) q[2];
sx q[2];
rz(2.4748477) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0181737) q[1];
sx q[1];
rz(-1.5416073) q[1];
sx q[1];
rz(-1.1694714) q[1];
x q[2];
rz(-2.2993777) q[3];
sx q[3];
rz(-0.73916736) q[3];
sx q[3];
rz(-0.30876866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83170825) q[2];
sx q[2];
rz(-0.94253057) q[2];
sx q[2];
rz(1.3523098) q[2];
rz(0.74812198) q[3];
sx q[3];
rz(-1.3040521) q[3];
sx q[3];
rz(-1.5266533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8378976) q[0];
sx q[0];
rz(-2.5717773) q[0];
sx q[0];
rz(1.8001528) q[0];
rz(-2.1603284) q[1];
sx q[1];
rz(-1.2813247) q[1];
sx q[1];
rz(2.431638) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2316596) q[0];
sx q[0];
rz(-1.2450448) q[0];
sx q[0];
rz(1.5909821) q[0];
rz(-pi) q[1];
rz(-0.23858647) q[2];
sx q[2];
rz(-2.1338531) q[2];
sx q[2];
rz(-0.67367879) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4932351) q[1];
sx q[1];
rz(-1.4005737) q[1];
sx q[1];
rz(-1.3551718) q[1];
rz(-pi) q[2];
rz(3.0026765) q[3];
sx q[3];
rz(-0.60231042) q[3];
sx q[3];
rz(1.7642782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.183341) q[2];
sx q[2];
rz(-0.7603344) q[2];
sx q[2];
rz(2.23488) q[2];
rz(-2.9592311) q[3];
sx q[3];
rz(-1.4401108) q[3];
sx q[3];
rz(-0.85842925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.3458503) q[0];
sx q[0];
rz(-2.1868571) q[0];
sx q[0];
rz(-0.20172754) q[0];
rz(1.3564823) q[1];
sx q[1];
rz(-2.4701665) q[1];
sx q[1];
rz(1.1511525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042554458) q[0];
sx q[0];
rz(-1.7793613) q[0];
sx q[0];
rz(2.0046356) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0315122) q[2];
sx q[2];
rz(-2.1270129) q[2];
sx q[2];
rz(-2.1002901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6342625) q[1];
sx q[1];
rz(-1.9729904) q[1];
sx q[1];
rz(-2.5172154) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3428976) q[3];
sx q[3];
rz(-1.0712475) q[3];
sx q[3];
rz(1.4394906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1543033) q[2];
sx q[2];
rz(-0.58649784) q[2];
sx q[2];
rz(-2.7941373) q[2];
rz(0.82184982) q[3];
sx q[3];
rz(-1.3653267) q[3];
sx q[3];
rz(-1.52389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8007941) q[0];
sx q[0];
rz(-2.4947385) q[0];
sx q[0];
rz(-2.0183753) q[0];
rz(0.42876354) q[1];
sx q[1];
rz(-0.88974297) q[1];
sx q[1];
rz(-0.14895983) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70388795) q[0];
sx q[0];
rz(-2.0876679) q[0];
sx q[0];
rz(-0.025512841) q[0];
rz(-2.8455714) q[2];
sx q[2];
rz(-1.3775577) q[2];
sx q[2];
rz(-3.0770965) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.625714) q[1];
sx q[1];
rz(-2.3209954) q[1];
sx q[1];
rz(-1.7734115) q[1];
rz(1.6520874) q[3];
sx q[3];
rz(-2.7346161) q[3];
sx q[3];
rz(-1.5291884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.011977) q[2];
sx q[2];
rz(-2.7836383) q[2];
sx q[2];
rz(2.820106) q[2];
rz(1.1164411) q[3];
sx q[3];
rz(-2.2827086) q[3];
sx q[3];
rz(1.4139676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1534934) q[0];
sx q[0];
rz(-1.7970002) q[0];
sx q[0];
rz(3.1186812) q[0];
rz(1.5494391) q[1];
sx q[1];
rz(-2.0982845) q[1];
sx q[1];
rz(-1.2124088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7717944) q[0];
sx q[0];
rz(-1.8317458) q[0];
sx q[0];
rz(-2.1290995) q[0];
rz(-pi) q[1];
rz(1.171807) q[2];
sx q[2];
rz(-0.771847) q[2];
sx q[2];
rz(-0.10153025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.072619) q[1];
sx q[1];
rz(-0.5104465) q[1];
sx q[1];
rz(-2.4859395) q[1];
x q[2];
rz(1.4597662) q[3];
sx q[3];
rz(-0.71655203) q[3];
sx q[3];
rz(-0.46886231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93078485) q[2];
sx q[2];
rz(-0.19889861) q[2];
sx q[2];
rz(2.1442294) q[2];
rz(1.2584244) q[3];
sx q[3];
rz(-1.1910028) q[3];
sx q[3];
rz(-1.9112126) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1579943) q[0];
sx q[0];
rz(-0.65584922) q[0];
sx q[0];
rz(2.6672145) q[0];
rz(-0.55167088) q[1];
sx q[1];
rz(-2.3730979) q[1];
sx q[1];
rz(-0.65753585) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5112572) q[0];
sx q[0];
rz(-1.5310107) q[0];
sx q[0];
rz(1.5944832) q[0];
rz(-pi) q[1];
rz(2.9716597) q[2];
sx q[2];
rz(-1.1701772) q[2];
sx q[2];
rz(-1.1378433) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8847981) q[1];
sx q[1];
rz(-1.9391283) q[1];
sx q[1];
rz(0.017466768) q[1];
rz(-1.4262496) q[3];
sx q[3];
rz(-1.4667505) q[3];
sx q[3];
rz(2.7964829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.23003301) q[2];
sx q[2];
rz(-2.7342789) q[2];
sx q[2];
rz(-1.3113021) q[2];
rz(2.4274872) q[3];
sx q[3];
rz(-1.2823391) q[3];
sx q[3];
rz(-2.5717521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1353726) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(2.5332992) q[0];
rz(-0.81781203) q[1];
sx q[1];
rz(-1.1759956) q[1];
sx q[1];
rz(2.3086595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.333182) q[0];
sx q[0];
rz(-1.6224553) q[0];
sx q[0];
rz(2.7201061) q[0];
x q[1];
rz(-2.8274546) q[2];
sx q[2];
rz(-1.8153662) q[2];
sx q[2];
rz(-0.42942522) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2257833) q[1];
sx q[1];
rz(-1.7034917) q[1];
sx q[1];
rz(-2.5578294) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5394443) q[3];
sx q[3];
rz(-2.1636181) q[3];
sx q[3];
rz(-1.624799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1828764) q[2];
sx q[2];
rz(-2.7887838) q[2];
sx q[2];
rz(0.58615169) q[2];
rz(2.2385249) q[3];
sx q[3];
rz(-2.51913) q[3];
sx q[3];
rz(-0.085845145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2387954) q[0];
sx q[0];
rz(-1.1302523) q[0];
sx q[0];
rz(-0.94373066) q[0];
rz(1.0974274) q[1];
sx q[1];
rz(-0.25249093) q[1];
sx q[1];
rz(2.9846334) q[1];
rz(-0.13931724) q[2];
sx q[2];
rz(-1.2323772) q[2];
sx q[2];
rz(2.7215794) q[2];
rz(-2.2678866) q[3];
sx q[3];
rz(-1.9601964) q[3];
sx q[3];
rz(-1.8662069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
