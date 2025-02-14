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
rz(0.42143917) q[0];
sx q[0];
rz(-1.2546189) q[0];
sx q[0];
rz(0.52809554) q[0];
rz(2.7948607) q[1];
sx q[1];
rz(-1.3277227) q[1];
sx q[1];
rz(-2.858207) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58812773) q[0];
sx q[0];
rz(-1.4736543) q[0];
sx q[0];
rz(2.5573822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47094496) q[2];
sx q[2];
rz(-2.2631553) q[2];
sx q[2];
rz(0.67519855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8979075) q[1];
sx q[1];
rz(-0.3966316) q[1];
sx q[1];
rz(-0.70901386) q[1];
rz(-pi) q[2];
rz(-1.9730522) q[3];
sx q[3];
rz(-1.4732971) q[3];
sx q[3];
rz(1.1643226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1957207) q[2];
sx q[2];
rz(-1.4737782) q[2];
sx q[2];
rz(0.26738581) q[2];
rz(0.17313677) q[3];
sx q[3];
rz(-2.0045529) q[3];
sx q[3];
rz(-0.6828298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6807569) q[0];
sx q[0];
rz(-0.84592485) q[0];
sx q[0];
rz(2.8435006) q[0];
rz(0.2671034) q[1];
sx q[1];
rz(-2.2593081) q[1];
sx q[1];
rz(0.90046901) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6783645) q[0];
sx q[0];
rz(-1.8186388) q[0];
sx q[0];
rz(-0.24893399) q[0];
x q[1];
rz(1.4203594) q[2];
sx q[2];
rz(-1.7152046) q[2];
sx q[2];
rz(0.23061801) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9760321) q[1];
sx q[1];
rz(-2.4270279) q[1];
sx q[1];
rz(0.30010414) q[1];
x q[2];
rz(-2.4839738) q[3];
sx q[3];
rz(-1.8557576) q[3];
sx q[3];
rz(-0.62510014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4233826) q[2];
sx q[2];
rz(-2.155828) q[2];
sx q[2];
rz(2.5441489) q[2];
rz(-1.0159703) q[3];
sx q[3];
rz(-1.4710602) q[3];
sx q[3];
rz(-1.5731251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.51282561) q[0];
sx q[0];
rz(-0.91575423) q[0];
sx q[0];
rz(3.0765007) q[0];
rz(2.0386631) q[1];
sx q[1];
rz(-1.8850336) q[1];
sx q[1];
rz(0.38806134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93174879) q[0];
sx q[0];
rz(-1.3375306) q[0];
sx q[0];
rz(-2.363028) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36696649) q[2];
sx q[2];
rz(-1.6410368) q[2];
sx q[2];
rz(-0.63634576) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.90012031) q[1];
sx q[1];
rz(-2.6992976) q[1];
sx q[1];
rz(1.9458179) q[1];
x q[2];
rz(-2.2476419) q[3];
sx q[3];
rz(-2.7778447) q[3];
sx q[3];
rz(-1.8574024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.97543657) q[2];
sx q[2];
rz(-2.1076951) q[2];
sx q[2];
rz(-0.088509716) q[2];
rz(-2.4837808) q[3];
sx q[3];
rz(-0.43840539) q[3];
sx q[3];
rz(-1.8678166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(0.45519644) q[0];
sx q[0];
rz(-2.1786067) q[0];
sx q[0];
rz(-0.92940593) q[0];
rz(1.1691673) q[1];
sx q[1];
rz(-0.87375748) q[1];
sx q[1];
rz(3.015231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1226144) q[0];
sx q[0];
rz(-0.41176957) q[0];
sx q[0];
rz(1.7420578) q[0];
rz(2.8123145) q[2];
sx q[2];
rz(-1.3751404) q[2];
sx q[2];
rz(-2.8242127) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1388701) q[1];
sx q[1];
rz(-2.5926081) q[1];
sx q[1];
rz(0.61985888) q[1];
rz(-2.5764719) q[3];
sx q[3];
rz(-2.6905746) q[3];
sx q[3];
rz(-2.7489611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.29752842) q[2];
sx q[2];
rz(-1.0928096) q[2];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70668689) q[0];
sx q[0];
rz(-2.4770985) q[0];
sx q[0];
rz(1.8988761) q[0];
rz(0.13936123) q[1];
sx q[1];
rz(-2.8262973) q[1];
sx q[1];
rz(2.349966) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2682429) q[0];
sx q[0];
rz(-1.6082544) q[0];
sx q[0];
rz(1.7245049) q[0];
rz(-pi) q[1];
rz(0.27288066) q[2];
sx q[2];
rz(-2.5562393) q[2];
sx q[2];
rz(1.7493474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5593223) q[1];
sx q[1];
rz(-0.30836654) q[1];
sx q[1];
rz(1.6496303) q[1];
rz(-0.96782622) q[3];
sx q[3];
rz(-1.669291) q[3];
sx q[3];
rz(-0.28901327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7255154) q[2];
sx q[2];
rz(-1.5308341) q[2];
sx q[2];
rz(-1.7804954) q[2];
rz(2.1040037) q[3];
sx q[3];
rz(-1.6890084) q[3];
sx q[3];
rz(3.1328372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0233651) q[0];
sx q[0];
rz(-0.26015493) q[0];
sx q[0];
rz(0.17876974) q[0];
rz(0.44890064) q[1];
sx q[1];
rz(-1.1652378) q[1];
sx q[1];
rz(-1.1753722) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7473874) q[0];
sx q[0];
rz(-2.3763467) q[0];
sx q[0];
rz(0.09455643) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.825338) q[2];
sx q[2];
rz(-2.3582728) q[2];
sx q[2];
rz(-2.6412727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6028831) q[1];
sx q[1];
rz(-2.3228163) q[1];
sx q[1];
rz(2.0497889) q[1];
rz(1.8038347) q[3];
sx q[3];
rz(-2.0762328) q[3];
sx q[3];
rz(-2.4580372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7281404) q[2];
sx q[2];
rz(-0.38273013) q[2];
sx q[2];
rz(-0.21391301) q[2];
rz(-1.4296069) q[3];
sx q[3];
rz(-1.6702646) q[3];
sx q[3];
rz(-1.3235929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0578617) q[0];
sx q[0];
rz(-1.2849176) q[0];
sx q[0];
rz(-0.65823746) q[0];
rz(-1.4312875) q[1];
sx q[1];
rz(-2.2555916) q[1];
sx q[1];
rz(-1.5828945) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4207672) q[0];
sx q[0];
rz(-1.1325192) q[0];
sx q[0];
rz(-1.363561) q[0];
rz(0.30679254) q[2];
sx q[2];
rz(-1.791583) q[2];
sx q[2];
rz(0.38635269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78608785) q[1];
sx q[1];
rz(-1.2337323) q[1];
sx q[1];
rz(-2.8028767) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1794847) q[3];
sx q[3];
rz(-1.0265304) q[3];
sx q[3];
rz(0.63718347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7610641) q[2];
sx q[2];
rz(-1.5498127) q[2];
sx q[2];
rz(1.2065678) q[2];
rz(-1.0208463) q[3];
sx q[3];
rz(-2.6267093) q[3];
sx q[3];
rz(-0.90768874) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9970488) q[0];
sx q[0];
rz(-2.9289065) q[0];
sx q[0];
rz(2.8295243) q[0];
rz(2.8128305) q[1];
sx q[1];
rz(-1.6203251) q[1];
sx q[1];
rz(-0.75334466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67344224) q[0];
sx q[0];
rz(-0.99981013) q[0];
sx q[0];
rz(-2.2500413) q[0];
rz(-0.021696731) q[2];
sx q[2];
rz(-1.0602131) q[2];
sx q[2];
rz(1.7006602) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6623913) q[1];
sx q[1];
rz(-2.0298784) q[1];
sx q[1];
rz(2.4513112) q[1];
rz(-pi) q[2];
rz(0.83266074) q[3];
sx q[3];
rz(-0.7023905) q[3];
sx q[3];
rz(-0.47846068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83851492) q[2];
sx q[2];
rz(-1.3884156) q[2];
sx q[2];
rz(-2.370749) q[2];
rz(-0.54909697) q[3];
sx q[3];
rz(-1.2884459) q[3];
sx q[3];
rz(-2.2652266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86641208) q[0];
sx q[0];
rz(-2.488945) q[0];
sx q[0];
rz(1.8021679) q[0];
rz(2.0088947) q[1];
sx q[1];
rz(-2.7481672) q[1];
sx q[1];
rz(-0.045230953) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6409174) q[0];
sx q[0];
rz(-1.7637327) q[0];
sx q[0];
rz(1.1290324) q[0];
rz(3.0196683) q[2];
sx q[2];
rz(-1.2164494) q[2];
sx q[2];
rz(-1.6533802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4457827) q[1];
sx q[1];
rz(-2.8572866) q[1];
sx q[1];
rz(-2.6587756) q[1];
rz(-pi) q[2];
rz(-1.0186152) q[3];
sx q[3];
rz(-1.3613627) q[3];
sx q[3];
rz(2.6554573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.905978) q[2];
sx q[2];
rz(-2.2048042) q[2];
sx q[2];
rz(2.3255685) q[2];
rz(1.7488545) q[3];
sx q[3];
rz(-1.1002898) q[3];
sx q[3];
rz(-1.5131153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.8428962) q[0];
sx q[0];
rz(-2.5923634) q[0];
sx q[0];
rz(-1.5604875) q[0];
rz(-0.16079482) q[1];
sx q[1];
rz(-2.000688) q[1];
sx q[1];
rz(-2.2198832) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9773437) q[0];
sx q[0];
rz(-1.5613587) q[0];
sx q[0];
rz(1.759985) q[0];
rz(-pi) q[1];
rz(2.3733487) q[2];
sx q[2];
rz(-2.3934954) q[2];
sx q[2];
rz(1.0988826) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.900093) q[1];
sx q[1];
rz(-1.422736) q[1];
sx q[1];
rz(-1.2847406) q[1];
rz(1.4170189) q[3];
sx q[3];
rz(-1.3594846) q[3];
sx q[3];
rz(-2.7696424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.86772743) q[2];
sx q[2];
rz(-2.2160857) q[2];
sx q[2];
rz(2.0694536) q[2];
rz(0.22465651) q[3];
sx q[3];
rz(-1.6499237) q[3];
sx q[3];
rz(2.2255161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0199725) q[0];
sx q[0];
rz(-0.92588035) q[0];
sx q[0];
rz(0.31582381) q[0];
rz(0.074180457) q[1];
sx q[1];
rz(-1.0646432) q[1];
sx q[1];
rz(-0.021312996) q[1];
rz(-2.3580768) q[2];
sx q[2];
rz(-1.5903689) q[2];
sx q[2];
rz(-1.5533536) q[2];
rz(-0.40423468) q[3];
sx q[3];
rz(-2.3334097) q[3];
sx q[3];
rz(2.4142053) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
