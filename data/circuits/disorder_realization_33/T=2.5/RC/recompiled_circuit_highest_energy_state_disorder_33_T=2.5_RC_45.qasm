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
rz(-2.9774732) q[0];
sx q[0];
rz(-1.0241221) q[0];
sx q[0];
rz(2.6574988) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(-1.5634544) q[1];
sx q[1];
rz(-0.054952316) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9021716) q[0];
sx q[0];
rz(-0.67089426) q[0];
sx q[0];
rz(-1.5205199) q[0];
x q[1];
rz(-1.820716) q[2];
sx q[2];
rz(-1.5139765) q[2];
sx q[2];
rz(-2.394258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7701227) q[1];
sx q[1];
rz(-0.16487637) q[1];
sx q[1];
rz(-2.5046964) q[1];
rz(-0.94934978) q[3];
sx q[3];
rz(-1.2934522) q[3];
sx q[3];
rz(0.30888939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83076465) q[2];
sx q[2];
rz(-0.97042933) q[2];
sx q[2];
rz(2.8669299) q[2];
rz(-0.87485391) q[3];
sx q[3];
rz(-1.091205) q[3];
sx q[3];
rz(-0.27354512) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4817568) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(0.39875317) q[0];
rz(2.8975471) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(1.1057378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.669848) q[0];
sx q[0];
rz(-2.6378184) q[0];
sx q[0];
rz(-0.48724799) q[0];
rz(-pi) q[1];
rz(-0.87200882) q[2];
sx q[2];
rz(-2.3559743) q[2];
sx q[2];
rz(-0.18839041) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2011021) q[1];
sx q[1];
rz(-3.1374212) q[1];
sx q[1];
rz(-0.68912403) q[1];
rz(-1.8375754) q[3];
sx q[3];
rz(-1.4543797) q[3];
sx q[3];
rz(2.4698225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4411321) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(2.0665118) q[2];
rz(0.69617802) q[3];
sx q[3];
rz(-1.2735561) q[3];
sx q[3];
rz(2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8657846) q[0];
sx q[0];
rz(-1.8901261) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(0.61947852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3724646) q[0];
sx q[0];
rz(-1.541168) q[0];
sx q[0];
rz(1.6086786) q[0];
x q[1];
rz(1.6939075) q[2];
sx q[2];
rz(-1.5242212) q[2];
sx q[2];
rz(-1.87093) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9668321) q[1];
sx q[1];
rz(-1.7836387) q[1];
sx q[1];
rz(1.879746) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1501174) q[3];
sx q[3];
rz(-1.0414755) q[3];
sx q[3];
rz(-2.4059699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2565903) q[2];
sx q[2];
rz(-1.8424748) q[2];
sx q[2];
rz(-0.090506434) q[2];
rz(0.45804405) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(1.1080144) q[3];
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
rz(-1.3312382) q[0];
sx q[0];
rz(-0.88307035) q[0];
sx q[0];
rz(2.2099387) q[0];
rz(-2.9725507) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(-2.1021252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0213892) q[0];
sx q[0];
rz(-1.3303555) q[0];
sx q[0];
rz(-0.73657764) q[0];
rz(0.2257963) q[2];
sx q[2];
rz(-1.3232627) q[2];
sx q[2];
rz(-0.087866656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66450602) q[1];
sx q[1];
rz(-1.6305822) q[1];
sx q[1];
rz(1.3633481) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63757293) q[3];
sx q[3];
rz(-2.7927783) q[3];
sx q[3];
rz(1.82774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52183759) q[2];
sx q[2];
rz(-2.0401185) q[2];
sx q[2];
rz(-0.49631611) q[2];
rz(2.3450092) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(-1.0930141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26688823) q[0];
sx q[0];
rz(-2.5570091) q[0];
sx q[0];
rz(-2.1697178) q[0];
rz(-3.019849) q[1];
sx q[1];
rz(-1.5309445) q[1];
sx q[1];
rz(-1.3294539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15518256) q[0];
sx q[0];
rz(-1.2647795) q[0];
sx q[0];
rz(-0.06432342) q[0];
rz(2.2586063) q[2];
sx q[2];
rz(-2.5341883) q[2];
sx q[2];
rz(-2.4508053) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.54534528) q[1];
sx q[1];
rz(-1.9831428) q[1];
sx q[1];
rz(-2.418787) q[1];
rz(1.6869808) q[3];
sx q[3];
rz(-3.0619762) q[3];
sx q[3];
rz(1.5594359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10151265) q[2];
sx q[2];
rz(-1.9209346) q[2];
sx q[2];
rz(0.15138781) q[2];
rz(2.3769412) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5449272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695892) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(-2.0701011) q[0];
rz(-0.53120652) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(-0.62613097) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3817978) q[0];
sx q[0];
rz(-3.1239428) q[0];
sx q[0];
rz(1.8438898) q[0];
rz(0.35289571) q[2];
sx q[2];
rz(-2.5145686) q[2];
sx q[2];
rz(-0.74909808) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28865675) q[1];
sx q[1];
rz(-0.96785883) q[1];
sx q[1];
rz(-1.9441685) q[1];
x q[2];
rz(0.71685426) q[3];
sx q[3];
rz(-0.62590137) q[3];
sx q[3];
rz(-0.38960534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35817394) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(-2.1273071) q[2];
rz(-1.8861534) q[3];
sx q[3];
rz(-1.7858601) q[3];
sx q[3];
rz(1.0534508) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96100539) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(-2.2210806) q[0];
rz(3.0275184) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(2.0106409) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81984542) q[0];
sx q[0];
rz(-1.5764569) q[0];
sx q[0];
rz(-0.59771718) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3384096) q[2];
sx q[2];
rz(-2.4121662) q[2];
sx q[2];
rz(-2.8620697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29207001) q[1];
sx q[1];
rz(-0.91150586) q[1];
sx q[1];
rz(2.9738722) q[1];
rz(3.00131) q[3];
sx q[3];
rz(-0.52061235) q[3];
sx q[3];
rz(-1.3365819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7941234) q[2];
sx q[2];
rz(-0.64212126) q[2];
sx q[2];
rz(-0.40880173) q[2];
rz(2.9583904) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(-2.0028152) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278397) q[0];
sx q[0];
rz(-0.04700679) q[0];
sx q[0];
rz(2.8507932) q[0];
rz(-0.94789061) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(-1.9416169) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3481846) q[0];
sx q[0];
rz(-2.4538592) q[0];
sx q[0];
rz(-2.4853112) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3340204) q[2];
sx q[2];
rz(-2.117273) q[2];
sx q[2];
rz(-1.228491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34741286) q[1];
sx q[1];
rz(-2.844226) q[1];
sx q[1];
rz(1.7296687) q[1];
rz(-pi) q[2];
rz(1.7539976) q[3];
sx q[3];
rz(-1.4716513) q[3];
sx q[3];
rz(0.20618901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7802508) q[2];
sx q[2];
rz(-1.6929408) q[2];
sx q[2];
rz(-1.7928436) q[2];
rz(1.5032984) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(-0.1851113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1143484) q[0];
sx q[0];
rz(-0.37186563) q[0];
sx q[0];
rz(2.0674904) q[0];
rz(0.4298003) q[1];
sx q[1];
rz(-2.0975515) q[1];
sx q[1];
rz(0.46357402) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060424711) q[0];
sx q[0];
rz(-1.3827818) q[0];
sx q[0];
rz(-0.096676143) q[0];
rz(-pi) q[1];
rz(0.26185449) q[2];
sx q[2];
rz(-1.3758278) q[2];
sx q[2];
rz(-1.5972114) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9755755) q[1];
sx q[1];
rz(-2.3834627) q[1];
sx q[1];
rz(2.1339971) q[1];
rz(-pi) q[2];
rz(-2.998407) q[3];
sx q[3];
rz(-0.56326635) q[3];
sx q[3];
rz(-1.8031424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.27434906) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(0.80844936) q[2];
rz(0.48458734) q[3];
sx q[3];
rz(-0.37823585) q[3];
sx q[3];
rz(-2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7782068) q[0];
sx q[0];
rz(-0.45553842) q[0];
sx q[0];
rz(2.0090012) q[0];
rz(0.038381902) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(2.8448232) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39958056) q[0];
sx q[0];
rz(-2.098408) q[0];
sx q[0];
rz(-0.3216089) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77905853) q[2];
sx q[2];
rz(-2.546306) q[2];
sx q[2];
rz(2.1434458) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6925507) q[1];
sx q[1];
rz(-1.7925279) q[1];
sx q[1];
rz(1.4936844) q[1];
rz(-pi) q[2];
rz(-0.39408306) q[3];
sx q[3];
rz(-0.9235477) q[3];
sx q[3];
rz(2.955472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9639637) q[2];
sx q[2];
rz(-2.0202426) q[2];
sx q[2];
rz(0.56619823) q[2];
rz(-2.0171793) q[3];
sx q[3];
rz(-0.12523139) q[3];
sx q[3];
rz(-1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2531256) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(0.26168564) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(-1.546312) q[2];
sx q[2];
rz(-1.6527805) q[2];
sx q[2];
rz(2.3775227) q[2];
rz(-0.86033173) q[3];
sx q[3];
rz(-1.0515778) q[3];
sx q[3];
rz(-0.39099494) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
