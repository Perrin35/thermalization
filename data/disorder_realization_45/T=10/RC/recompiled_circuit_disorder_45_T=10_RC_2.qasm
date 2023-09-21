OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(2.4581576) q[0];
sx q[0];
rz(12.0876) q[0];
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(-0.64136139) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7862579) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(-1.8281561) q[0];
rz(-1.1429206) q[2];
sx q[2];
rz(-1.2245721) q[2];
sx q[2];
rz(1.8510173) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8909059) q[1];
sx q[1];
rz(-0.85039925) q[1];
sx q[1];
rz(0.62178639) q[1];
rz(-0.62946837) q[3];
sx q[3];
rz(-1.2347722) q[3];
sx q[3];
rz(-1.1050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.550094) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(-0.47810289) q[2];
rz(1.6889307) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(-2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(-2.1226728) q[0];
rz(-1.6628751) q[1];
sx q[1];
rz(-0.61518413) q[1];
sx q[1];
rz(-2.5085124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3490152) q[0];
sx q[0];
rz(-1.4903869) q[0];
sx q[0];
rz(-0.76064674) q[0];
rz(-0.622153) q[2];
sx q[2];
rz(-2.2894147) q[2];
sx q[2];
rz(-2.6785786) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2411023) q[1];
sx q[1];
rz(-1.9111269) q[1];
sx q[1];
rz(-1.1452922) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88926104) q[3];
sx q[3];
rz(-1.325843) q[3];
sx q[3];
rz(-2.3219061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7818266) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(-0.11745545) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(-1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85787073) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(-0.088407956) q[0];
rz(1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(-2.450768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2937677) q[0];
sx q[0];
rz(-2.9710037) q[0];
sx q[0];
rz(-0.45844309) q[0];
x q[1];
rz(2.3766741) q[2];
sx q[2];
rz(-1.3996482) q[2];
sx q[2];
rz(2.314832) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9835811) q[1];
sx q[1];
rz(-2.9576655) q[1];
sx q[1];
rz(1.7137101) q[1];
rz(-pi) q[2];
rz(1.2191804) q[3];
sx q[3];
rz(-1.9295441) q[3];
sx q[3];
rz(-2.2946332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2174125) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(-3.0351191) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(0.52247125) q[0];
rz(-0.32896313) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(2.5879588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5579917) q[0];
sx q[0];
rz(-0.86692536) q[0];
sx q[0];
rz(2.5224199) q[0];
rz(-0.76769464) q[2];
sx q[2];
rz(-2.1622217) q[2];
sx q[2];
rz(-0.64085811) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5428634) q[1];
sx q[1];
rz(-0.79303629) q[1];
sx q[1];
rz(1.383177) q[1];
rz(-pi) q[2];
rz(0.028285154) q[3];
sx q[3];
rz(-0.79704282) q[3];
sx q[3];
rz(2.459211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.557495) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(0.49475691) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(-1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6506127) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(-2.0565128) q[0];
rz(-0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(-0.18403149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83443123) q[0];
sx q[0];
rz(-0.30858985) q[0];
sx q[0];
rz(0.35085268) q[0];
rz(2.7765772) q[2];
sx q[2];
rz(-1.3114309) q[2];
sx q[2];
rz(0.58155453) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3369493) q[1];
sx q[1];
rz(-0.78362432) q[1];
sx q[1];
rz(-2.4451392) q[1];
rz(-pi) q[2];
rz(2.9156978) q[3];
sx q[3];
rz(-2.4028006) q[3];
sx q[3];
rz(-1.6177288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(1.1408172) q[2];
rz(-0.42282894) q[3];
sx q[3];
rz(-1.4774277) q[3];
sx q[3];
rz(1.1269349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7891156) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(0.88678962) q[0];
rz(-2.9011762) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(2.9930847) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7667023) q[0];
sx q[0];
rz(-2.3935211) q[0];
sx q[0];
rz(-0.23776777) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6988585) q[2];
sx q[2];
rz(-0.25207439) q[2];
sx q[2];
rz(2.4661951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8588184) q[1];
sx q[1];
rz(-0.28206477) q[1];
sx q[1];
rz(1.1689405) q[1];
rz(-1.1735736) q[3];
sx q[3];
rz(-1.1363315) q[3];
sx q[3];
rz(-1.7064427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.285816) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(2.6532069) q[2];
rz(-2.6521902) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(1.874812) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4483036) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(2.3235902) q[0];
rz(-1.2524293) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(1.12524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97911191) q[0];
sx q[0];
rz(-1.1971803) q[0];
sx q[0];
rz(0.11066779) q[0];
rz(-0.15525012) q[2];
sx q[2];
rz(-1.8660188) q[2];
sx q[2];
rz(-1.0709907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.352467) q[1];
sx q[1];
rz(-1.2965634) q[1];
sx q[1];
rz(-1.2624361) q[1];
x q[2];
rz(-2.1358228) q[3];
sx q[3];
rz(-2.1566448) q[3];
sx q[3];
rz(-0.98423959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0730878) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(-0.64458624) q[2];
rz(1.4792431) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1709764) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(-1.9203141) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(-2.1309526) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0897652) q[0];
sx q[0];
rz(-0.83501378) q[0];
sx q[0];
rz(-2.9249973) q[0];
x q[1];
rz(-1.0534361) q[2];
sx q[2];
rz(-2.5854955) q[2];
sx q[2];
rz(2.9964787) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6992221) q[1];
sx q[1];
rz(-1.3500299) q[1];
sx q[1];
rz(-1.3352331) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6628077) q[3];
sx q[3];
rz(-1.9698471) q[3];
sx q[3];
rz(1.240318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6179787) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(-0.69407216) q[2];
rz(-0.56898919) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(0.49474299) q[0];
rz(2.5231979) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(3.0659952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15007524) q[0];
sx q[0];
rz(-2.4872094) q[0];
sx q[0];
rz(-1.0432711) q[0];
rz(-pi) q[1];
rz(1.3952903) q[2];
sx q[2];
rz(-1.5117466) q[2];
sx q[2];
rz(-2.4621071) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1739993) q[1];
sx q[1];
rz(-2.446397) q[1];
sx q[1];
rz(-2.980152) q[1];
x q[2];
rz(-1.3887614) q[3];
sx q[3];
rz(-0.958003) q[3];
sx q[3];
rz(0.36422563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33346924) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(1.3405651) q[2];
rz(-2.8373485) q[3];
sx q[3];
rz(-2.2777568) q[3];
sx q[3];
rz(-0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2952404) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(-0.17679581) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-2.5071564) q[1];
sx q[1];
rz(2.7005844) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84687418) q[0];
sx q[0];
rz(-0.28781578) q[0];
sx q[0];
rz(-0.78886445) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6454234) q[2];
sx q[2];
rz(-1.7241038) q[2];
sx q[2];
rz(-0.15664936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1383789) q[1];
sx q[1];
rz(-2.5795476) q[1];
sx q[1];
rz(-2.2467062) q[1];
rz(-pi) q[2];
rz(-0.62426626) q[3];
sx q[3];
rz(-1.8161976) q[3];
sx q[3];
rz(-1.0898958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56197721) q[2];
sx q[2];
rz(-0.56869555) q[2];
sx q[2];
rz(0.30187541) q[2];
rz(2.2484696) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(-0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72398913) q[0];
sx q[0];
rz(-1.8483193) q[0];
sx q[0];
rz(1.666477) q[0];
rz(0.026731116) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(-2.2970207) q[2];
sx q[2];
rz(-1.9324586) q[2];
sx q[2];
rz(2.4697138) q[2];
rz(-1.0711014) q[3];
sx q[3];
rz(-1.0715967) q[3];
sx q[3];
rz(-1.788492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];