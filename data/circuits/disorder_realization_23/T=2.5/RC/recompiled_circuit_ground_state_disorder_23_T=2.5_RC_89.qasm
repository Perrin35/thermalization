OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.58148122) q[0];
sx q[0];
rz(-1.886263) q[0];
sx q[0];
rz(2.6502967) q[0];
rz(-0.30839857) q[1];
sx q[1];
rz(-1.8412794) q[1];
sx q[1];
rz(-0.058606776) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9973544) q[0];
sx q[0];
rz(-0.73737007) q[0];
sx q[0];
rz(2.2585346) q[0];
rz(2.357448) q[2];
sx q[2];
rz(-1.9040739) q[2];
sx q[2];
rz(-0.12034697) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8303918) q[1];
sx q[1];
rz(-1.7313384) q[1];
sx q[1];
rz(2.83648) q[1];
rz(1.8080072) q[3];
sx q[3];
rz(-0.32792842) q[3];
sx q[3];
rz(1.6746317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0356174) q[2];
sx q[2];
rz(-2.3507037) q[2];
sx q[2];
rz(-2.7604575) q[2];
rz(-3.1086339) q[3];
sx q[3];
rz(-1.4573174) q[3];
sx q[3];
rz(-1.0491252) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9829262) q[0];
sx q[0];
rz(-0.52678147) q[0];
sx q[0];
rz(0.43989936) q[0];
rz(-0.061821763) q[1];
sx q[1];
rz(-1.5301751) q[1];
sx q[1];
rz(0.27438146) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3456164) q[0];
sx q[0];
rz(-1.1152949) q[0];
sx q[0];
rz(-2.2296395) q[0];
rz(-2.9129162) q[2];
sx q[2];
rz(-1.9878721) q[2];
sx q[2];
rz(-1.4859391) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.741863) q[1];
sx q[1];
rz(-1.4951839) q[1];
sx q[1];
rz(1.1221328) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8058978) q[3];
sx q[3];
rz(-1.0367852) q[3];
sx q[3];
rz(-1.4015568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8476094) q[2];
sx q[2];
rz(-0.01394883) q[2];
sx q[2];
rz(0.070276109) q[2];
rz(-2.516732) q[3];
sx q[3];
rz(-1.8127541) q[3];
sx q[3];
rz(3.0666053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.4570419) q[0];
sx q[0];
rz(-1.6383189) q[0];
sx q[0];
rz(-2.7497838) q[0];
rz(2.1458972) q[1];
sx q[1];
rz(-2.3052146) q[1];
sx q[1];
rz(3.0973184) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0842154) q[0];
sx q[0];
rz(-2.4829933) q[0];
sx q[0];
rz(1.0530217) q[0];
rz(2.2703553) q[2];
sx q[2];
rz(-0.92513621) q[2];
sx q[2];
rz(-1.439656) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4391195) q[1];
sx q[1];
rz(-2.8608192) q[1];
sx q[1];
rz(2.4149815) q[1];
x q[2];
rz(0.1416154) q[3];
sx q[3];
rz(-0.60185233) q[3];
sx q[3];
rz(2.600649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5135045) q[2];
sx q[2];
rz(-0.91614437) q[2];
sx q[2];
rz(-1.2543031) q[2];
rz(3.0301376) q[3];
sx q[3];
rz(-2.1696551) q[3];
sx q[3];
rz(-0.85533065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8348389) q[0];
sx q[0];
rz(-2.4618537) q[0];
sx q[0];
rz(0.87598959) q[0];
rz(2.7442979) q[1];
sx q[1];
rz(-1.5715503) q[1];
sx q[1];
rz(1.3321336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0615342) q[0];
sx q[0];
rz(-1.3803015) q[0];
sx q[0];
rz(-0.47866042) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0202254) q[2];
sx q[2];
rz(-0.75006333) q[2];
sx q[2];
rz(-1.5062065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5046103) q[1];
sx q[1];
rz(-2.9688997) q[1];
sx q[1];
rz(-0.47422282) q[1];
rz(-1.0334084) q[3];
sx q[3];
rz(-1.3674481) q[3];
sx q[3];
rz(2.5739079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1490271) q[2];
sx q[2];
rz(-1.488648) q[2];
sx q[2];
rz(-1.7472902) q[2];
rz(0.99327883) q[3];
sx q[3];
rz(-1.1634049) q[3];
sx q[3];
rz(1.2629898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8472327) q[0];
sx q[0];
rz(-2.9672186) q[0];
sx q[0];
rz(-2.2827523) q[0];
rz(-0.42293388) q[1];
sx q[1];
rz(-0.49086389) q[1];
sx q[1];
rz(1.6805958) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92991023) q[0];
sx q[0];
rz(-2.0693464) q[0];
sx q[0];
rz(1.431026) q[0];
x q[1];
rz(-3.1390343) q[2];
sx q[2];
rz(-2.0160437) q[2];
sx q[2];
rz(1.3819577) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7778421) q[1];
sx q[1];
rz(-1.7763174) q[1];
sx q[1];
rz(-1.3745893) q[1];
rz(-pi) q[2];
rz(2.9988838) q[3];
sx q[3];
rz(-2.0721216) q[3];
sx q[3];
rz(3.0317948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6964263) q[2];
sx q[2];
rz(-2.7327765) q[2];
sx q[2];
rz(-1.6430631) q[2];
rz(2.0131352) q[3];
sx q[3];
rz(-1.107629) q[3];
sx q[3];
rz(0.64277738) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6978825) q[0];
sx q[0];
rz(-1.8118129) q[0];
sx q[0];
rz(-2.4238996) q[0];
rz(0.36092654) q[1];
sx q[1];
rz(-0.91054994) q[1];
sx q[1];
rz(2.8275729) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0180394) q[0];
sx q[0];
rz(-1.6121056) q[0];
sx q[0];
rz(-0.49561758) q[0];
rz(-pi) q[1];
rz(-1.6144361) q[2];
sx q[2];
rz(-0.72481643) q[2];
sx q[2];
rz(-2.5284655) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18923387) q[1];
sx q[1];
rz(-1.2554607) q[1];
sx q[1];
rz(2.0604412) q[1];
x q[2];
rz(-0.18743214) q[3];
sx q[3];
rz(-0.77875528) q[3];
sx q[3];
rz(-1.4608778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6641984) q[2];
sx q[2];
rz(-1.3496642) q[2];
sx q[2];
rz(-0.2142218) q[2];
rz(-2.2231806) q[3];
sx q[3];
rz(-2.1954506) q[3];
sx q[3];
rz(1.4001747) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4825831) q[0];
sx q[0];
rz(-0.56203401) q[0];
sx q[0];
rz(-2.5157978) q[0];
rz(-1.2311426) q[1];
sx q[1];
rz(-1.8742671) q[1];
sx q[1];
rz(-2.321718) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5441262) q[0];
sx q[0];
rz(-1.0080999) q[0];
sx q[0];
rz(1.94348) q[0];
rz(-2.4307239) q[2];
sx q[2];
rz(-1.9108624) q[2];
sx q[2];
rz(-2.1785871) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7149799) q[1];
sx q[1];
rz(-2.8016114) q[1];
sx q[1];
rz(-1.3424385) q[1];
rz(-pi) q[2];
rz(3.047057) q[3];
sx q[3];
rz(-1.5526287) q[3];
sx q[3];
rz(1.8985192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95255533) q[2];
sx q[2];
rz(-1.8841212) q[2];
sx q[2];
rz(0.38743585) q[2];
rz(0.45346013) q[3];
sx q[3];
rz(-0.87415868) q[3];
sx q[3];
rz(-0.15997729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979935) q[0];
sx q[0];
rz(-1.1389808) q[0];
sx q[0];
rz(-1.1238264) q[0];
rz(-0.75346142) q[1];
sx q[1];
rz(-1.1898142) q[1];
sx q[1];
rz(-1.9728647) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92896398) q[0];
sx q[0];
rz(-0.88200404) q[0];
sx q[0];
rz(-0.11879158) q[0];
x q[1];
rz(1.5828791) q[2];
sx q[2];
rz(-0.45736499) q[2];
sx q[2];
rz(-1.7810019) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1271229) q[1];
sx q[1];
rz(-1.8934087) q[1];
sx q[1];
rz(-2.2029976) q[1];
rz(-pi) q[2];
rz(2.1451352) q[3];
sx q[3];
rz(-0.74241168) q[3];
sx q[3];
rz(-1.3848828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0231102) q[2];
sx q[2];
rz(-1.0175984) q[2];
sx q[2];
rz(-2.6712096) q[2];
rz(1.6810301) q[3];
sx q[3];
rz(-1.2602256) q[3];
sx q[3];
rz(-0.98106658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39101446) q[0];
sx q[0];
rz(-1.7711201) q[0];
sx q[0];
rz(-1.4960666) q[0];
rz(1.1356614) q[1];
sx q[1];
rz(-1.2874425) q[1];
sx q[1];
rz(-0.48887238) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35621168) q[0];
sx q[0];
rz(-1.4309034) q[0];
sx q[0];
rz(-3.1091305) q[0];
rz(-pi) q[1];
x q[1];
rz(0.096285419) q[2];
sx q[2];
rz(-2.0803148) q[2];
sx q[2];
rz(-2.4063039) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5597361) q[1];
sx q[1];
rz(-1.5213517) q[1];
sx q[1];
rz(1.54848) q[1];
rz(2.5219157) q[3];
sx q[3];
rz(-1.5159199) q[3];
sx q[3];
rz(0.60407274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6966072) q[2];
sx q[2];
rz(-2.5312436) q[2];
sx q[2];
rz(-1.102591) q[2];
rz(-1.6164814) q[3];
sx q[3];
rz(-2.3307255) q[3];
sx q[3];
rz(-1.6703687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1333604) q[0];
sx q[0];
rz(-2.6120549) q[0];
sx q[0];
rz(-1.3326921) q[0];
rz(-0.38707271) q[1];
sx q[1];
rz(-1.2553071) q[1];
sx q[1];
rz(2.0852087) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41829506) q[0];
sx q[0];
rz(-1.3271062) q[0];
sx q[0];
rz(0.14362986) q[0];
x q[1];
rz(-1.6956639) q[2];
sx q[2];
rz(-1.1556632) q[2];
sx q[2];
rz(1.268569) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0936959) q[1];
sx q[1];
rz(-2.2655089) q[1];
sx q[1];
rz(1.0662778) q[1];
rz(-pi) q[2];
rz(-0.9659395) q[3];
sx q[3];
rz(-0.98839983) q[3];
sx q[3];
rz(2.8386311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.010783823) q[2];
sx q[2];
rz(-0.84785145) q[2];
sx q[2];
rz(1.3991785) q[2];
rz(-2.0684659) q[3];
sx q[3];
rz(-1.224204) q[3];
sx q[3];
rz(0.46260241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.241093) q[0];
sx q[0];
rz(-1.6443962) q[0];
sx q[0];
rz(-1.6045438) q[0];
rz(2.4201139) q[1];
sx q[1];
rz(-2.571048) q[1];
sx q[1];
rz(-1.2068988) q[1];
rz(-0.73312326) q[2];
sx q[2];
rz(-0.55126581) q[2];
sx q[2];
rz(0.049964213) q[2];
rz(-1.6695475) q[3];
sx q[3];
rz(-0.80174123) q[3];
sx q[3];
rz(1.5855736) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
