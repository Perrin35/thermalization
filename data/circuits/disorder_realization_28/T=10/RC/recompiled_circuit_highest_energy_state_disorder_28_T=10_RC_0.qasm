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
rz(-2.7201535) q[0];
sx q[0];
rz(-1.8869737) q[0];
sx q[0];
rz(-0.52809554) q[0];
rz(-0.34673196) q[1];
sx q[1];
rz(4.4693153) q[1];
sx q[1];
rz(9.1413924) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91862667) q[0];
sx q[0];
rz(-2.1518927) q[0];
sx q[0];
rz(-1.6870935) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0714226) q[2];
sx q[2];
rz(-0.81489043) q[2];
sx q[2];
rz(-1.7931149) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6470601) q[1];
sx q[1];
rz(-1.2732098) q[1];
sx q[1];
rz(-1.8370166) q[1];
x q[2];
rz(-1.3259726) q[3];
sx q[3];
rz(-2.7283165) q[3];
sx q[3];
rz(-0.63140376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9458719) q[2];
sx q[2];
rz(-1.4737782) q[2];
sx q[2];
rz(-0.26738581) q[2];
rz(0.17313677) q[3];
sx q[3];
rz(-1.1370398) q[3];
sx q[3];
rz(0.6828298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46083573) q[0];
sx q[0];
rz(-0.84592485) q[0];
sx q[0];
rz(0.2980921) q[0];
rz(-2.8744892) q[1];
sx q[1];
rz(-2.2593081) q[1];
sx q[1];
rz(0.90046901) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6783645) q[0];
sx q[0];
rz(-1.3229538) q[0];
sx q[0];
rz(-0.24893399) q[0];
rz(-1.4203594) q[2];
sx q[2];
rz(-1.4263881) q[2];
sx q[2];
rz(0.23061801) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.63486491) q[1];
sx q[1];
rz(-1.7657451) q[1];
sx q[1];
rz(-0.69201236) q[1];
rz(1.2163085) q[3];
sx q[3];
rz(-2.1976314) q[3];
sx q[3];
rz(1.1594866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4233826) q[2];
sx q[2];
rz(-2.155828) q[2];
sx q[2];
rz(0.59744376) q[2];
rz(-1.0159703) q[3];
sx q[3];
rz(-1.6705325) q[3];
sx q[3];
rz(-1.5684675) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628767) q[0];
sx q[0];
rz(-0.91575423) q[0];
sx q[0];
rz(0.065091982) q[0];
rz(2.0386631) q[1];
sx q[1];
rz(-1.256559) q[1];
sx q[1];
rz(2.7535313) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267305) q[0];
sx q[0];
rz(-0.81860733) q[0];
sx q[0];
rz(1.2486876) q[0];
rz(0.36696649) q[2];
sx q[2];
rz(-1.5005558) q[2];
sx q[2];
rz(0.63634576) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0124828) q[1];
sx q[1];
rz(-1.4133682) q[1];
sx q[1];
rz(-1.1557259) q[1];
rz(-pi) q[2];
rz(-0.23406844) q[3];
sx q[3];
rz(-1.2897629) q[3];
sx q[3];
rz(0.57408702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1661561) q[2];
sx q[2];
rz(-2.1076951) q[2];
sx q[2];
rz(-0.088509716) q[2];
rz(-2.4837808) q[3];
sx q[3];
rz(-0.43840539) q[3];
sx q[3];
rz(1.2737761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6863962) q[0];
sx q[0];
rz(-2.1786067) q[0];
sx q[0];
rz(-2.2121867) q[0];
rz(1.9724253) q[1];
sx q[1];
rz(-0.87375748) q[1];
sx q[1];
rz(-3.015231) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9360775) q[0];
sx q[0];
rz(-1.9761855) q[0];
sx q[0];
rz(0.074294093) q[0];
x q[1];
rz(-2.5917327) q[2];
sx q[2];
rz(-2.7603995) q[2];
sx q[2];
rz(1.3709244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3060386) q[1];
sx q[1];
rz(-1.1321196) q[1];
sx q[1];
rz(-1.2293596) q[1];
x q[2];
rz(0.56512079) q[3];
sx q[3];
rz(-2.6905746) q[3];
sx q[3];
rz(-2.7489611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29752842) q[2];
sx q[2];
rz(-2.0487831) q[2];
sx q[2];
rz(-1.4275985) q[2];
rz(2.0370238) q[3];
sx q[3];
rz(-1.5449056) q[3];
sx q[3];
rz(2.4692718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70668689) q[0];
sx q[0];
rz(-0.66449419) q[0];
sx q[0];
rz(1.8988761) q[0];
rz(-0.13936123) q[1];
sx q[1];
rz(-0.31529537) q[1];
sx q[1];
rz(-0.79162663) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46027943) q[0];
sx q[0];
rz(-0.15817197) q[0];
sx q[0];
rz(1.3307443) q[0];
x q[1];
rz(-1.3940195) q[2];
sx q[2];
rz(-1.0097754) q[2];
sx q[2];
rz(-1.068312) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0549889) q[1];
sx q[1];
rz(-1.5468925) q[1];
sx q[1];
rz(1.2633282) q[1];
rz(-pi) q[2];
x q[2];
rz(1.743312) q[3];
sx q[3];
rz(-0.60997395) q[3];
sx q[3];
rz(2.0016746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41607729) q[2];
sx q[2];
rz(-1.6107586) q[2];
sx q[2];
rz(1.7804954) q[2];
rz(-2.1040037) q[3];
sx q[3];
rz(-1.4525843) q[3];
sx q[3];
rz(3.1328372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0233651) q[0];
sx q[0];
rz(-0.26015493) q[0];
sx q[0];
rz(0.17876974) q[0];
rz(-0.44890064) q[1];
sx q[1];
rz(-1.9763549) q[1];
sx q[1];
rz(-1.1753722) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3942053) q[0];
sx q[0];
rz(-0.76524599) q[0];
sx q[0];
rz(-0.09455643) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8959013) q[2];
sx q[2];
rz(-0.81907228) q[2];
sx q[2];
rz(2.9931835) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6028831) q[1];
sx q[1];
rz(-2.3228163) q[1];
sx q[1];
rz(-2.0497889) q[1];
rz(0.5171295) q[3];
sx q[3];
rz(-1.7742566) q[3];
sx q[3];
rz(0.77282016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41345227) q[2];
sx q[2];
rz(-2.7588625) q[2];
sx q[2];
rz(-2.9276796) q[2];
rz(-1.7119857) q[3];
sx q[3];
rz(-1.471328) q[3];
sx q[3];
rz(-1.3235929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0578617) q[0];
sx q[0];
rz(-1.8566751) q[0];
sx q[0];
rz(-0.65823746) q[0];
rz(-1.4312875) q[1];
sx q[1];
rz(-0.88600102) q[1];
sx q[1];
rz(1.5828945) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3806136) q[0];
sx q[0];
rz(-1.7581994) q[0];
sx q[0];
rz(-0.44664573) q[0];
x q[1];
rz(-1.8020242) q[2];
sx q[2];
rz(-1.8699081) q[2];
sx q[2];
rz(2.0264152) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66872193) q[1];
sx q[1];
rz(-1.2518365) q[1];
sx q[1];
rz(1.215056) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56230265) q[3];
sx q[3];
rz(-0.65861693) q[3];
sx q[3];
rz(-0.035619481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7610641) q[2];
sx q[2];
rz(-1.5917799) q[2];
sx q[2];
rz(1.9350249) q[2];
rz(-1.0208463) q[3];
sx q[3];
rz(-2.6267093) q[3];
sx q[3];
rz(-0.90768874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9970488) q[0];
sx q[0];
rz(-2.9289065) q[0];
sx q[0];
rz(-2.8295243) q[0];
rz(-2.8128305) q[1];
sx q[1];
rz(-1.6203251) q[1];
sx q[1];
rz(-2.388248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3088209) q[0];
sx q[0];
rz(-2.1276908) q[0];
sx q[0];
rz(-0.69016074) q[0];
x q[1];
rz(1.5320832) q[2];
sx q[2];
rz(-2.6305894) q[2];
sx q[2];
rz(-1.3965565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7406941) q[1];
sx q[1];
rz(-2.3339865) q[1];
sx q[1];
rz(-0.66011564) q[1];
x q[2];
rz(-2.6238717) q[3];
sx q[3];
rz(-2.0690527) q[3];
sx q[3];
rz(-0.39420989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83851492) q[2];
sx q[2];
rz(-1.7531771) q[2];
sx q[2];
rz(2.370749) q[2];
rz(0.54909697) q[3];
sx q[3];
rz(-1.8531468) q[3];
sx q[3];
rz(-2.2652266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86641208) q[0];
sx q[0];
rz(-0.65264767) q[0];
sx q[0];
rz(1.3394248) q[0];
rz(-1.132698) q[1];
sx q[1];
rz(-2.7481672) q[1];
sx q[1];
rz(3.0963617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6409174) q[0];
sx q[0];
rz(-1.7637327) q[0];
sx q[0];
rz(-1.1290324) q[0];
rz(3.0196683) q[2];
sx q[2];
rz(-1.2164494) q[2];
sx q[2];
rz(-1.6533802) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69581) q[1];
sx q[1];
rz(-0.28430609) q[1];
sx q[1];
rz(-0.48281702) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24465235) q[3];
sx q[3];
rz(-2.1095697) q[3];
sx q[3];
rz(-0.95726171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.905978) q[2];
sx q[2];
rz(-2.2048042) q[2];
sx q[2];
rz(2.3255685) q[2];
rz(1.3927381) q[3];
sx q[3];
rz(-1.1002898) q[3];
sx q[3];
rz(-1.6284774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8428962) q[0];
sx q[0];
rz(-0.5492292) q[0];
sx q[0];
rz(1.5811051) q[0];
rz(2.9807978) q[1];
sx q[1];
rz(-1.1409047) q[1];
sx q[1];
rz(2.2198832) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35730001) q[0];
sx q[0];
rz(-2.9521715) q[0];
sx q[0];
rz(-1.5206536) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3733487) q[2];
sx q[2];
rz(-0.74809725) q[2];
sx q[2];
rz(2.04271) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79419151) q[1];
sx q[1];
rz(-0.32116613) q[1];
sx q[1];
rz(-2.0570517) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7245737) q[3];
sx q[3];
rz(-1.3594846) q[3];
sx q[3];
rz(2.7696424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2738652) q[2];
sx q[2];
rz(-2.2160857) q[2];
sx q[2];
rz(-2.0694536) q[2];
rz(-0.22465651) q[3];
sx q[3];
rz(-1.6499237) q[3];
sx q[3];
rz(0.91607654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.12162019) q[0];
sx q[0];
rz(-0.92588035) q[0];
sx q[0];
rz(0.31582381) q[0];
rz(-0.074180457) q[1];
sx q[1];
rz(-2.0769495) q[1];
sx q[1];
rz(3.1202797) q[1];
rz(0.027728524) q[2];
sx q[2];
rz(-2.3578845) q[2];
sx q[2];
rz(3.1393928) q[2];
rz(1.9613135) q[3];
sx q[3];
rz(-2.2979679) q[3];
sx q[3];
rz(-1.2818492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
