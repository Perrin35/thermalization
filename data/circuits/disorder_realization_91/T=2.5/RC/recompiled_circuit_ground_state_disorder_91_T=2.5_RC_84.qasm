OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.33136) q[0];
sx q[0];
rz(-1.2547837) q[0];
sx q[0];
rz(0.64594185) q[0];
rz(-1.4589925) q[1];
sx q[1];
rz(-0.12812935) q[1];
sx q[1];
rz(2.473414) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7433919) q[0];
sx q[0];
rz(-1.8132134) q[0];
sx q[0];
rz(-1.4546118) q[0];
rz(-1.1485841) q[2];
sx q[2];
rz(-0.6919043) q[2];
sx q[2];
rz(0.9949323) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9757802) q[1];
sx q[1];
rz(-2.4620612) q[1];
sx q[1];
rz(-0.91109101) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1124437) q[3];
sx q[3];
rz(-2.3437269) q[3];
sx q[3];
rz(-3.0276379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51840034) q[2];
sx q[2];
rz(-0.65541583) q[2];
sx q[2];
rz(2.0584959) q[2];
rz(2.8862503) q[3];
sx q[3];
rz(-0.89699236) q[3];
sx q[3];
rz(1.159509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.979368) q[0];
sx q[0];
rz(-3.0001682) q[0];
sx q[0];
rz(-0.033578385) q[0];
rz(2.0032739) q[1];
sx q[1];
rz(-2.1470224) q[1];
sx q[1];
rz(-0.23695645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7342252) q[0];
sx q[0];
rz(-1.7091284) q[0];
sx q[0];
rz(-2.9107735) q[0];
rz(-pi) q[1];
rz(1.9230546) q[2];
sx q[2];
rz(-0.50261231) q[2];
sx q[2];
rz(-1.9723122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0396871) q[1];
sx q[1];
rz(-0.19409212) q[1];
sx q[1];
rz(1.7029352) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2036425) q[3];
sx q[3];
rz(-0.68162912) q[3];
sx q[3];
rz(1.5860032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.61243764) q[2];
sx q[2];
rz(-1.0470904) q[2];
sx q[2];
rz(-3.0212413) q[2];
rz(2.555661) q[3];
sx q[3];
rz(-1.3804932) q[3];
sx q[3];
rz(-1.8732635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2633857) q[0];
sx q[0];
rz(-0.96579856) q[0];
sx q[0];
rz(2.1534488) q[0];
rz(1.6506317) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(1.5536701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8107306) q[0];
sx q[0];
rz(-1.5548883) q[0];
sx q[0];
rz(-1.5718979) q[0];
rz(-pi) q[1];
rz(-0.13440172) q[2];
sx q[2];
rz(-1.5024868) q[2];
sx q[2];
rz(1.7640863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37582477) q[1];
sx q[1];
rz(-1.8319523) q[1];
sx q[1];
rz(0.029726083) q[1];
rz(-1.9918898) q[3];
sx q[3];
rz(-2.097887) q[3];
sx q[3];
rz(0.89638174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3042018) q[2];
sx q[2];
rz(-1.7981671) q[2];
sx q[2];
rz(0.070579441) q[2];
rz(-2.6523759) q[3];
sx q[3];
rz(-0.990812) q[3];
sx q[3];
rz(-0.14744559) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1864784) q[0];
sx q[0];
rz(-2.4461353) q[0];
sx q[0];
rz(1.4008993) q[0];
rz(-2.9715624) q[1];
sx q[1];
rz(-1.4100807) q[1];
sx q[1];
rz(1.2331351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2182541) q[0];
sx q[0];
rz(-1.4460576) q[0];
sx q[0];
rz(-1.1472923) q[0];
rz(2.4242086) q[2];
sx q[2];
rz(-1.1108494) q[2];
sx q[2];
rz(1.1954952) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68076949) q[1];
sx q[1];
rz(-1.5857547) q[1];
sx q[1];
rz(1.2911002) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6771116) q[3];
sx q[3];
rz(-1.5600403) q[3];
sx q[3];
rz(-1.8821723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3418545) q[2];
sx q[2];
rz(-1.4350472) q[2];
sx q[2];
rz(-2.5119761) q[2];
rz(3.0047505) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(-1.7262044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28302309) q[0];
sx q[0];
rz(-2.6333599) q[0];
sx q[0];
rz(3.0392905) q[0];
rz(1.4981859) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(0.35710517) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5095561) q[0];
sx q[0];
rz(-0.3118383) q[0];
sx q[0];
rz(-2.3760892) q[0];
rz(-pi) q[1];
rz(-1.5979366) q[2];
sx q[2];
rz(-1.6939591) q[2];
sx q[2];
rz(-2.9371098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6372924) q[1];
sx q[1];
rz(-1.2724814) q[1];
sx q[1];
rz(0.53108414) q[1];
x q[2];
rz(1.0856241) q[3];
sx q[3];
rz(-1.1729122) q[3];
sx q[3];
rz(-2.0069938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47306481) q[2];
sx q[2];
rz(-1.4184971) q[2];
sx q[2];
rz(1.8190039) q[2];
rz(-1.4332917) q[3];
sx q[3];
rz(-1.8247484) q[3];
sx q[3];
rz(-2.9776261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54356164) q[0];
sx q[0];
rz(-0.99995166) q[0];
sx q[0];
rz(1.438197) q[0];
rz(-1.9012798) q[1];
sx q[1];
rz(-2.4867609) q[1];
sx q[1];
rz(1.3528489) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677296) q[0];
sx q[0];
rz(-2.7244748) q[0];
sx q[0];
rz(1.3365082) q[0];
x q[1];
rz(-0.2850432) q[2];
sx q[2];
rz(-0.54423287) q[2];
sx q[2];
rz(-2.2533992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6042418) q[1];
sx q[1];
rz(-2.4927995) q[1];
sx q[1];
rz(0.93334812) q[1];
rz(-pi) q[2];
rz(-0.39522533) q[3];
sx q[3];
rz(-1.6754284) q[3];
sx q[3];
rz(-2.6992309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9271586) q[2];
sx q[2];
rz(-0.36403251) q[2];
sx q[2];
rz(-2.7511609) q[2];
rz(1.3123784) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(0.81982476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44535962) q[0];
sx q[0];
rz(-0.92388988) q[0];
sx q[0];
rz(-1.7286638) q[0];
rz(-1.6784809) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(-0.63953343) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4832925) q[0];
sx q[0];
rz(-2.4118703) q[0];
sx q[0];
rz(-2.150282) q[0];
rz(-0.15070559) q[2];
sx q[2];
rz(-1.4960714) q[2];
sx q[2];
rz(0.13241235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6653135) q[1];
sx q[1];
rz(-1.8355882) q[1];
sx q[1];
rz(2.4864343) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9652548) q[3];
sx q[3];
rz(-1.4603851) q[3];
sx q[3];
rz(-2.4019456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21706906) q[2];
sx q[2];
rz(-1.3288493) q[2];
sx q[2];
rz(-2.6893943) q[2];
rz(-0.2229812) q[3];
sx q[3];
rz(-0.28135869) q[3];
sx q[3];
rz(0.15534672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.13101354) q[0];
sx q[0];
rz(-0.92229811) q[0];
sx q[0];
rz(-0.16192326) q[0];
rz(0.34795347) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(0.23439342) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1944626) q[0];
sx q[0];
rz(-0.85887733) q[0];
sx q[0];
rz(0.80545896) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4215464) q[2];
sx q[2];
rz(-1.8077501) q[2];
sx q[2];
rz(2.4958378) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1474851) q[1];
sx q[1];
rz(-0.6692769) q[1];
sx q[1];
rz(-1.2662751) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27674562) q[3];
sx q[3];
rz(-1.6349155) q[3];
sx q[3];
rz(-1.8217877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28635412) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(2.3983119) q[3];
sx q[3];
rz(-2.1478839) q[3];
sx q[3];
rz(0.68964094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8648935) q[0];
sx q[0];
rz(-1.445329) q[0];
sx q[0];
rz(1.0869166) q[0];
rz(-1.9089606) q[1];
sx q[1];
rz(-1.107736) q[1];
sx q[1];
rz(-1.3023652) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3051599) q[0];
sx q[0];
rz(-1.4576685) q[0];
sx q[0];
rz(-0.21048429) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97424284) q[2];
sx q[2];
rz(-1.8678209) q[2];
sx q[2];
rz(-3.0926306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86587807) q[1];
sx q[1];
rz(-0.68531407) q[1];
sx q[1];
rz(-0.65071836) q[1];
rz(-pi) q[2];
rz(1.5518673) q[3];
sx q[3];
rz(-2.193748) q[3];
sx q[3];
rz(1.2894693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46633259) q[2];
sx q[2];
rz(-1.2559428) q[2];
sx q[2];
rz(-0.41435286) q[2];
rz(0.77053344) q[3];
sx q[3];
rz(-1.4221752) q[3];
sx q[3];
rz(0.93366247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3122124) q[0];
sx q[0];
rz(-0.82013622) q[0];
sx q[0];
rz(0.46863753) q[0];
rz(-1.8679484) q[1];
sx q[1];
rz(-1.2761152) q[1];
sx q[1];
rz(-2.830107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92665206) q[0];
sx q[0];
rz(-0.72065852) q[0];
sx q[0];
rz(-1.65833) q[0];
rz(0.63705541) q[2];
sx q[2];
rz(-0.14547507) q[2];
sx q[2];
rz(-1.5381952) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35986082) q[1];
sx q[1];
rz(-1.5241429) q[1];
sx q[1];
rz(-1.0830888) q[1];
rz(-pi) q[2];
rz(2.4021637) q[3];
sx q[3];
rz(-0.66777705) q[3];
sx q[3];
rz(1.8831933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1062539) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(1.7485471) q[2];
rz(-0.92750183) q[3];
sx q[3];
rz(-1.5774957) q[3];
sx q[3];
rz(-0.86021304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8504234) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(-2.9604079) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(1.3109315) q[2];
sx q[2];
rz(-0.54878546) q[2];
sx q[2];
rz(2.6981163) q[2];
rz(-0.90418935) q[3];
sx q[3];
rz(-1.0720357) q[3];
sx q[3];
rz(-2.8513249) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
