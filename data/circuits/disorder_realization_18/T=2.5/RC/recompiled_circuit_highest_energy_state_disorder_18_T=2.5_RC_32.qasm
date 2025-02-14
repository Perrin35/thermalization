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
rz(-1.9189605) q[0];
sx q[0];
rz(-0.55813342) q[0];
sx q[0];
rz(0.65879917) q[0];
rz(0.88762033) q[1];
sx q[1];
rz(-0.69488156) q[1];
sx q[1];
rz(-3.053009) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026042875) q[0];
sx q[0];
rz(-0.47240531) q[0];
sx q[0];
rz(-1.1360854) q[0];
rz(-pi) q[1];
rz(-0.93841432) q[2];
sx q[2];
rz(-2.4749196) q[2];
sx q[2];
rz(-3.0423903) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7329053) q[1];
sx q[1];
rz(-0.63893203) q[1];
sx q[1];
rz(1.8270135) q[1];
rz(-2.1351571) q[3];
sx q[3];
rz(-0.86216037) q[3];
sx q[3];
rz(2.9624727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52458557) q[2];
sx q[2];
rz(-1.969939) q[2];
sx q[2];
rz(2.6535772) q[2];
rz(-0.46767849) q[3];
sx q[3];
rz(-0.52507639) q[3];
sx q[3];
rz(2.974143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2117598) q[0];
sx q[0];
rz(-0.81120315) q[0];
sx q[0];
rz(-1.8875341) q[0];
rz(-2.5414741) q[1];
sx q[1];
rz(-1.3328726) q[1];
sx q[1];
rz(2.7235203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1195898) q[0];
sx q[0];
rz(-1.3922577) q[0];
sx q[0];
rz(-0.54212944) q[0];
x q[1];
rz(2.977191) q[2];
sx q[2];
rz(-0.65509701) q[2];
sx q[2];
rz(3.0384105) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0954803) q[1];
sx q[1];
rz(-0.32278341) q[1];
sx q[1];
rz(1.7782636) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90466212) q[3];
sx q[3];
rz(-2.2719529) q[3];
sx q[3];
rz(-1.2010434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7347001) q[2];
sx q[2];
rz(-1.1310581) q[2];
sx q[2];
rz(0.1184173) q[2];
rz(-2.7393869) q[3];
sx q[3];
rz(-1.6536568) q[3];
sx q[3];
rz(1.3291298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6953832) q[0];
sx q[0];
rz(-0.24676794) q[0];
sx q[0];
rz(1.9870019) q[0];
rz(1.062695) q[1];
sx q[1];
rz(-2.1176391) q[1];
sx q[1];
rz(1.9564995) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1405615) q[0];
sx q[0];
rz(-2.2375282) q[0];
sx q[0];
rz(-0.35463766) q[0];
rz(-pi) q[1];
rz(0.22151557) q[2];
sx q[2];
rz(-1.0660401) q[2];
sx q[2];
rz(-1.5919478) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9251057) q[1];
sx q[1];
rz(-0.47678927) q[1];
sx q[1];
rz(0.9245704) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5873379) q[3];
sx q[3];
rz(-0.46559428) q[3];
sx q[3];
rz(1.3265644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2350754) q[2];
sx q[2];
rz(-2.2791028) q[2];
sx q[2];
rz(1.6974576) q[2];
rz(2.2200572) q[3];
sx q[3];
rz(-2.5369365) q[3];
sx q[3];
rz(0.055189341) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4398572) q[0];
sx q[0];
rz(-3.0210962) q[0];
sx q[0];
rz(-0.078027092) q[0];
rz(-1.1174508) q[1];
sx q[1];
rz(-1.8945339) q[1];
sx q[1];
rz(1.6695581) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32737752) q[0];
sx q[0];
rz(-2.3348885) q[0];
sx q[0];
rz(-2.9931158) q[0];
rz(0.83547893) q[2];
sx q[2];
rz(-0.72670055) q[2];
sx q[2];
rz(-1.6146605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8824655) q[1];
sx q[1];
rz(-1.5992016) q[1];
sx q[1];
rz(1.2106497) q[1];
rz(0.72301282) q[3];
sx q[3];
rz(-1.8545056) q[3];
sx q[3];
rz(-2.8151388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3052519) q[2];
sx q[2];
rz(-0.64771104) q[2];
sx q[2];
rz(2.9913537) q[2];
rz(0.59208208) q[3];
sx q[3];
rz(-1.838107) q[3];
sx q[3];
rz(-1.1881812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2109461) q[0];
sx q[0];
rz(-0.35363126) q[0];
sx q[0];
rz(0.76552248) q[0];
rz(-0.80134478) q[1];
sx q[1];
rz(-1.067602) q[1];
sx q[1];
rz(-1.1917005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.649851) q[0];
sx q[0];
rz(-2.5132127) q[0];
sx q[0];
rz(-2.6091486) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52764738) q[2];
sx q[2];
rz(-1.9867267) q[2];
sx q[2];
rz(-1.6017583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6995865) q[1];
sx q[1];
rz(-2.4558459) q[1];
sx q[1];
rz(1.6722635) q[1];
rz(0.7804773) q[3];
sx q[3];
rz(-0.37364081) q[3];
sx q[3];
rz(-2.3821156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2661065) q[2];
sx q[2];
rz(-2.2849639) q[2];
sx q[2];
rz(0.53492707) q[2];
rz(-2.5180425) q[3];
sx q[3];
rz(-1.1112735) q[3];
sx q[3];
rz(0.88304869) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542338) q[0];
sx q[0];
rz(-0.41340241) q[0];
sx q[0];
rz(2.2075388) q[0];
rz(1.8961228) q[1];
sx q[1];
rz(-1.586069) q[1];
sx q[1];
rz(2.5159871) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4419785) q[0];
sx q[0];
rz(-1.7480339) q[0];
sx q[0];
rz(-0.81798792) q[0];
x q[1];
rz(-2.2410237) q[2];
sx q[2];
rz(-1.4975542) q[2];
sx q[2];
rz(-1.9894969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8865247) q[1];
sx q[1];
rz(-2.2400948) q[1];
sx q[1];
rz(1.1978576) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17600893) q[3];
sx q[3];
rz(-1.7932682) q[3];
sx q[3];
rz(-2.5493718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.049456747) q[2];
sx q[2];
rz(-0.9149887) q[2];
sx q[2];
rz(-3.0165239) q[2];
rz(2.4140029) q[3];
sx q[3];
rz(-1.5743419) q[3];
sx q[3];
rz(-2.6653813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8759988) q[0];
sx q[0];
rz(-0.47530526) q[0];
sx q[0];
rz(-1.884961) q[0];
rz(1.1988634) q[1];
sx q[1];
rz(-1.4224982) q[1];
sx q[1];
rz(1.2748324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0555818) q[0];
sx q[0];
rz(-1.1390098) q[0];
sx q[0];
rz(2.5837269) q[0];
rz(2.7706128) q[2];
sx q[2];
rz(-1.0643963) q[2];
sx q[2];
rz(-1.9547403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5526455) q[1];
sx q[1];
rz(-1.1324404) q[1];
sx q[1];
rz(2.752315) q[1];
x q[2];
rz(-1.9406464) q[3];
sx q[3];
rz(-2.8027034) q[3];
sx q[3];
rz(1.3047993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4510497) q[2];
sx q[2];
rz(-2.1284911) q[2];
sx q[2];
rz(2.1742382) q[2];
rz(-1.5830154) q[3];
sx q[3];
rz(-1.6396921) q[3];
sx q[3];
rz(-0.30174747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3485182) q[0];
sx q[0];
rz(-2.5418042) q[0];
sx q[0];
rz(3.0317958) q[0];
rz(-1.173136) q[1];
sx q[1];
rz(-2.1497612) q[1];
sx q[1];
rz(2.0593624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17324461) q[0];
sx q[0];
rz(-0.81615198) q[0];
sx q[0];
rz(-1.9108755) q[0];
rz(-pi) q[1];
rz(-0.10854461) q[2];
sx q[2];
rz(-1.0883779) q[2];
sx q[2];
rz(0.71297836) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6761192) q[1];
sx q[1];
rz(-0.88485588) q[1];
sx q[1];
rz(-0.53829389) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5030062) q[3];
sx q[3];
rz(-1.2632086) q[3];
sx q[3];
rz(1.547055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3829696) q[2];
sx q[2];
rz(-1.0926282) q[2];
sx q[2];
rz(0.36337241) q[2];
rz(2.162497) q[3];
sx q[3];
rz(-2.4978814) q[3];
sx q[3];
rz(-2.6399829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30524224) q[0];
sx q[0];
rz(-0.79638052) q[0];
sx q[0];
rz(0.35261944) q[0];
rz(-1.7800219) q[1];
sx q[1];
rz(-1.9510599) q[1];
sx q[1];
rz(-0.1415267) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3985719) q[0];
sx q[0];
rz(-1.5367265) q[0];
sx q[0];
rz(2.9719246) q[0];
x q[1];
rz(-2.5472801) q[2];
sx q[2];
rz(-1.1804575) q[2];
sx q[2];
rz(0.28536404) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5474816) q[1];
sx q[1];
rz(-1.2073646) q[1];
sx q[1];
rz(0.44977447) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4961595) q[3];
sx q[3];
rz(-1.107599) q[3];
sx q[3];
rz(-2.2768159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8286459) q[2];
sx q[2];
rz(-1.6673648) q[2];
sx q[2];
rz(2.055577) q[2];
rz(-0.18276754) q[3];
sx q[3];
rz(-1.026231) q[3];
sx q[3];
rz(2.1695547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6163841) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(-0.64424789) q[0];
rz(0.066990189) q[1];
sx q[1];
rz(-1.7552152) q[1];
sx q[1];
rz(-3.0279874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841235) q[0];
sx q[0];
rz(-2.9175287) q[0];
sx q[0];
rz(-1.5518918) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.058931) q[2];
sx q[2];
rz(-1.6726521) q[2];
sx q[2];
rz(-1.5422623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1190824) q[1];
sx q[1];
rz(-1.4003229) q[1];
sx q[1];
rz(1.4524231) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4573649) q[3];
sx q[3];
rz(-0.3598752) q[3];
sx q[3];
rz(-2.7168093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1429448) q[2];
sx q[2];
rz(-1.8871658) q[2];
sx q[2];
rz(2.3764853) q[2];
rz(0.1782002) q[3];
sx q[3];
rz(-1.5604115) q[3];
sx q[3];
rz(0.83711886) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9255623) q[0];
sx q[0];
rz(-0.75556527) q[0];
sx q[0];
rz(2.8788993) q[0];
rz(0.33949159) q[1];
sx q[1];
rz(-1.9120293) q[1];
sx q[1];
rz(1.0614352) q[1];
rz(-0.4364261) q[2];
sx q[2];
rz(-2.438475) q[2];
sx q[2];
rz(-1.9154127) q[2];
rz(2.8398561) q[3];
sx q[3];
rz(-1.269125) q[3];
sx q[3];
rz(2.453809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
