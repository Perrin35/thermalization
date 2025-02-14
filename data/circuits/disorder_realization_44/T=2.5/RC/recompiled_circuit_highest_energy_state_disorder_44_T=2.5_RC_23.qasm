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
rz(-0.62701464) q[0];
sx q[0];
rz(-2.505317) q[0];
sx q[0];
rz(1.0778435) q[0];
rz(-2.0650504) q[1];
sx q[1];
rz(-2.512518) q[1];
sx q[1];
rz(-2.10973) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54993486) q[0];
sx q[0];
rz(-3.1374133) q[0];
sx q[0];
rz(3.0463534) q[0];
x q[1];
rz(2.5450942) q[2];
sx q[2];
rz(-0.41587999) q[2];
sx q[2];
rz(-2.2592827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.061432211) q[1];
sx q[1];
rz(-2.5974063) q[1];
sx q[1];
rz(-0.70744343) q[1];
rz(-2.9239465) q[3];
sx q[3];
rz(-2.7983694) q[3];
sx q[3];
rz(1.9453334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2304113) q[2];
sx q[2];
rz(-1.958467) q[2];
sx q[2];
rz(-2.6739056) q[2];
rz(-2.6905401) q[3];
sx q[3];
rz(-2.7428198) q[3];
sx q[3];
rz(0.27802813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69538799) q[0];
sx q[0];
rz(-0.22286335) q[0];
sx q[0];
rz(0.15431246) q[0];
rz(0.34126869) q[1];
sx q[1];
rz(-2.7879265) q[1];
sx q[1];
rz(-2.7214859) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.882928) q[0];
sx q[0];
rz(-2.1999584) q[0];
sx q[0];
rz(1.8101519) q[0];
x q[1];
rz(-1.5163563) q[2];
sx q[2];
rz(-2.4307361) q[2];
sx q[2];
rz(2.9756851) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5737916) q[1];
sx q[1];
rz(-2.670799) q[1];
sx q[1];
rz(3.1385697) q[1];
x q[2];
rz(0.16815925) q[3];
sx q[3];
rz(-2.3871867) q[3];
sx q[3];
rz(-1.8969632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61098617) q[2];
sx q[2];
rz(-0.89908081) q[2];
sx q[2];
rz(2.3685624) q[2];
rz(0.84345877) q[3];
sx q[3];
rz(-1.5638899) q[3];
sx q[3];
rz(-2.5696866) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.968349) q[0];
sx q[0];
rz(-1.0272212) q[0];
sx q[0];
rz(-0.20208836) q[0];
rz(-2.659722) q[1];
sx q[1];
rz(-0.20129573) q[1];
sx q[1];
rz(-0.906382) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2807686) q[0];
sx q[0];
rz(-1.1204136) q[0];
sx q[0];
rz(0.76669873) q[0];
rz(1.1163122) q[2];
sx q[2];
rz(-1.2589728) q[2];
sx q[2];
rz(-0.21909595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14531586) q[1];
sx q[1];
rz(-1.7309124) q[1];
sx q[1];
rz(-2.5969863) q[1];
rz(-pi) q[2];
rz(-0.1679095) q[3];
sx q[3];
rz(-1.2793667) q[3];
sx q[3];
rz(1.7999906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0704982) q[2];
sx q[2];
rz(-2.0281894) q[2];
sx q[2];
rz(-0.58007288) q[2];
rz(-1.5602559) q[3];
sx q[3];
rz(-1.7507078) q[3];
sx q[3];
rz(-1.4238547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4105014) q[0];
sx q[0];
rz(-3.0803362) q[0];
sx q[0];
rz(-3.0938003) q[0];
rz(-1.2740678) q[1];
sx q[1];
rz(-1.2893226) q[1];
sx q[1];
rz(0.045225708) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6134098) q[0];
sx q[0];
rz(-0.65887132) q[0];
sx q[0];
rz(-1.8152366) q[0];
rz(-pi) q[1];
rz(2.8130262) q[2];
sx q[2];
rz(-0.96146482) q[2];
sx q[2];
rz(2.6549465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.2312551) q[1];
sx q[1];
rz(-1.5738166) q[1];
sx q[1];
rz(1.5652204) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0257231) q[3];
sx q[3];
rz(-1.9506321) q[3];
sx q[3];
rz(-1.4295417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4101326) q[2];
sx q[2];
rz(-1.8803909) q[2];
sx q[2];
rz(-2.251808) q[2];
rz(0.54640031) q[3];
sx q[3];
rz(-1.0011287) q[3];
sx q[3];
rz(2.8304097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6379717) q[0];
sx q[0];
rz(-1.4542955) q[0];
sx q[0];
rz(1.7244435) q[0];
rz(-2.0218938) q[1];
sx q[1];
rz(-1.7599301) q[1];
sx q[1];
rz(-1.959257) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9989026) q[0];
sx q[0];
rz(-1.9574528) q[0];
sx q[0];
rz(-1.645734) q[0];
x q[1];
rz(-1.117302) q[2];
sx q[2];
rz(-1.1536479) q[2];
sx q[2];
rz(0.81062775) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2675543) q[1];
sx q[1];
rz(-0.8505162) q[1];
sx q[1];
rz(-1.505386) q[1];
x q[2];
rz(-2.0145217) q[3];
sx q[3];
rz(-2.2699353) q[3];
sx q[3];
rz(0.8389896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.917439) q[2];
sx q[2];
rz(-0.76442337) q[2];
sx q[2];
rz(-0.42382851) q[2];
rz(-1.1118927) q[3];
sx q[3];
rz(-2.6424776) q[3];
sx q[3];
rz(-2.4901938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79913419) q[0];
sx q[0];
rz(-3.0693711) q[0];
sx q[0];
rz(0.95241958) q[0];
rz(0.77596387) q[1];
sx q[1];
rz(-2.2064078) q[1];
sx q[1];
rz(1.5752569) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0512038) q[0];
sx q[0];
rz(-3.0984833) q[0];
sx q[0];
rz(-0.48844047) q[0];
rz(3.1029019) q[2];
sx q[2];
rz(-0.96444791) q[2];
sx q[2];
rz(-2.4232466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4500931) q[1];
sx q[1];
rz(-2.1686692) q[1];
sx q[1];
rz(-1.3045425) q[1];
rz(-pi) q[2];
rz(1.2346047) q[3];
sx q[3];
rz(-0.51516525) q[3];
sx q[3];
rz(-0.060843918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39885193) q[2];
sx q[2];
rz(-2.3336918) q[2];
sx q[2];
rz(-0.74637949) q[2];
rz(1.0910723) q[3];
sx q[3];
rz(-1.4336136) q[3];
sx q[3];
rz(0.55238849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2090476) q[0];
sx q[0];
rz(-3.0952251) q[0];
sx q[0];
rz(2.5982502) q[0];
rz(-2.6047193) q[1];
sx q[1];
rz(-2.8165292) q[1];
sx q[1];
rz(3.0184025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7167599) q[0];
sx q[0];
rz(-0.93532978) q[0];
sx q[0];
rz(-3.1034971) q[0];
rz(-pi) q[1];
rz(-2.4692612) q[2];
sx q[2];
rz(-2.5200082) q[2];
sx q[2];
rz(-2.7171135) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1063042) q[1];
sx q[1];
rz(-1.0266487) q[1];
sx q[1];
rz(-0.72334163) q[1];
rz(-pi) q[2];
rz(-0.80135342) q[3];
sx q[3];
rz(-2.1983302) q[3];
sx q[3];
rz(2.827044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8366375) q[2];
sx q[2];
rz(-1.387549) q[2];
sx q[2];
rz(-0.48509625) q[2];
rz(1.1525611) q[3];
sx q[3];
rz(-1.7771114) q[3];
sx q[3];
rz(0.047957234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7570067) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(0.49569976) q[0];
rz(3.0777625) q[1];
sx q[1];
rz(-0.7494691) q[1];
sx q[1];
rz(0.071050342) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9793689) q[0];
sx q[0];
rz(-2.6155229) q[0];
sx q[0];
rz(-2.6499676) q[0];
rz(0.2500791) q[2];
sx q[2];
rz(-0.61998487) q[2];
sx q[2];
rz(2.300547) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7531705) q[1];
sx q[1];
rz(-2.8533397) q[1];
sx q[1];
rz(0.99797319) q[1];
rz(-pi) q[2];
rz(-0.44720113) q[3];
sx q[3];
rz(-1.4942823) q[3];
sx q[3];
rz(0.90257989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12585982) q[2];
sx q[2];
rz(-1.4733529) q[2];
sx q[2];
rz(-0.86137548) q[2];
rz(-1.8371948) q[3];
sx q[3];
rz(-2.5671037) q[3];
sx q[3];
rz(1.5931574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9651589) q[0];
sx q[0];
rz(-2.6506944) q[0];
sx q[0];
rz(1.4392256) q[0];
rz(3.0610415) q[1];
sx q[1];
rz(-1.6702024) q[1];
sx q[1];
rz(-1.7505987) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1192525) q[0];
sx q[0];
rz(-1.0588405) q[0];
sx q[0];
rz(0.48959022) q[0];
rz(-pi) q[1];
rz(1.6356231) q[2];
sx q[2];
rz(-2.8362084) q[2];
sx q[2];
rz(2.7750654) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94165588) q[1];
sx q[1];
rz(-0.98683954) q[1];
sx q[1];
rz(2.5879509) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7709522) q[3];
sx q[3];
rz(-1.2145059) q[3];
sx q[3];
rz(0.72539893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.78802687) q[2];
sx q[2];
rz(-1.2949508) q[2];
sx q[2];
rz(0.12727748) q[2];
rz(-0.92653972) q[3];
sx q[3];
rz(-0.24181952) q[3];
sx q[3];
rz(1.658879) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7552898) q[0];
sx q[0];
rz(-0.64978623) q[0];
sx q[0];
rz(-1.6408828) q[0];
rz(-0.81037784) q[1];
sx q[1];
rz(-0.85226285) q[1];
sx q[1];
rz(-2.3694029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9861575) q[0];
sx q[0];
rz(-2.2732908) q[0];
sx q[0];
rz(0.8849573) q[0];
rz(-2.5224546) q[2];
sx q[2];
rz(-0.35751128) q[2];
sx q[2];
rz(-0.50019568) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3898728) q[1];
sx q[1];
rz(-1.0214318) q[1];
sx q[1];
rz(3.1279039) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5072892) q[3];
sx q[3];
rz(-2.4606703) q[3];
sx q[3];
rz(1.6977967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3810252) q[2];
sx q[2];
rz(-0.80485359) q[2];
sx q[2];
rz(0.62925657) q[2];
rz(-0.73090807) q[3];
sx q[3];
rz(-2.2980502) q[3];
sx q[3];
rz(1.4430911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44169852) q[0];
sx q[0];
rz(-0.48342539) q[0];
sx q[0];
rz(-0.40406686) q[0];
rz(-0.40456698) q[1];
sx q[1];
rz(-1.7351983) q[1];
sx q[1];
rz(2.4529967) q[1];
rz(2.8597833) q[2];
sx q[2];
rz(-1.5055613) q[2];
sx q[2];
rz(0.67508634) q[2];
rz(1.2734781) q[3];
sx q[3];
rz(-2.086629) q[3];
sx q[3];
rz(1.119864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
