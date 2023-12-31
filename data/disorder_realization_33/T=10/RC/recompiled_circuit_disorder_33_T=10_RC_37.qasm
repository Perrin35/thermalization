OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(-2.2059031) q[1];
sx q[1];
rz(1.5712665) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9309064) q[0];
sx q[0];
rz(-1.2331729) q[0];
sx q[0];
rz(-2.7772285) q[0];
rz(0.76531305) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(-3.090976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5174487) q[1];
sx q[1];
rz(-1.9074719) q[1];
sx q[1];
rz(1.3144073) q[1];
rz(-pi) q[2];
rz(1.0899815) q[3];
sx q[3];
rz(-2.7365723) q[3];
sx q[3];
rz(2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.87542614) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(-1.1323294) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(-1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(2.9557513) q[0];
rz(2.5813685) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(-2.9247608) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6235979) q[0];
sx q[0];
rz(-1.8580084) q[0];
sx q[0];
rz(-0.90201305) q[0];
rz(-pi) q[1];
rz(-0.88044135) q[2];
sx q[2];
rz(-0.67697064) q[2];
sx q[2];
rz(1.134269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9358881) q[1];
sx q[1];
rz(-1.3591896) q[1];
sx q[1];
rz(-2.2536224) q[1];
rz(-pi) q[2];
rz(2.5012245) q[3];
sx q[3];
rz(-0.98207563) q[3];
sx q[3];
rz(2.0224188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(-1.2878093) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(-1.8151981) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-0.93260971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9661449) q[0];
sx q[0];
rz(-1.517059) q[0];
sx q[0];
rz(-1.2530112) q[0];
x q[1];
rz(1.3435752) q[2];
sx q[2];
rz(-1.3933239) q[2];
sx q[2];
rz(2.9858659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0886791) q[1];
sx q[1];
rz(-2.1267849) q[1];
sx q[1];
rz(1.710379) q[1];
rz(-0.91498418) q[3];
sx q[3];
rz(-2.4768156) q[3];
sx q[3];
rz(-1.8613601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(-2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(-0.50022593) q[0];
rz(-2.3362828) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(1.6436228) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26101199) q[0];
sx q[0];
rz(-1.2769165) q[0];
sx q[0];
rz(-1.0512933) q[0];
rz(-1.8565606) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(-0.49516585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8343463) q[1];
sx q[1];
rz(-1.3386968) q[1];
sx q[1];
rz(0.26718617) q[1];
x q[2];
rz(0.76969947) q[3];
sx q[3];
rz(-0.61004988) q[3];
sx q[3];
rz(1.1806012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3952289) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(2.4397819) q[2];
rz(-0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(-0.81533122) q[0];
rz(-1.6197846) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(1.048208) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002776) q[0];
sx q[0];
rz(-1.8198697) q[0];
sx q[0];
rz(-1.8427909) q[0];
rz(1.5830718) q[2];
sx q[2];
rz(-2.2222812) q[2];
sx q[2];
rz(-2.2001681) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9427467) q[1];
sx q[1];
rz(-1.8837351) q[1];
sx q[1];
rz(1.8206157) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1225554) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(-0.90941959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6158225) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(2.0416416) q[2];
rz(2.3163017) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(-2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.0681756) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-0.90240479) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(-0.12983233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8311365) q[0];
sx q[0];
rz(-1.2015011) q[0];
sx q[0];
rz(2.5184758) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36826276) q[2];
sx q[2];
rz(-2.1149181) q[2];
sx q[2];
rz(-2.1899109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2643913) q[1];
sx q[1];
rz(-2.1346722) q[1];
sx q[1];
rz(-2.0045723) q[1];
rz(-pi) q[2];
rz(2.5927605) q[3];
sx q[3];
rz(-1.4026814) q[3];
sx q[3];
rz(2.7041534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(-2.9373346) q[2];
rz(1.2060818) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4181353) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(-1.6947421) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(0.68626219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0745569) q[0];
sx q[0];
rz(-2.0445163) q[0];
sx q[0];
rz(2.0635598) q[0];
rz(-pi) q[1];
rz(2.1967728) q[2];
sx q[2];
rz(-2.2959024) q[2];
sx q[2];
rz(3.0996029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2974907) q[1];
sx q[1];
rz(-2.9417848) q[1];
sx q[1];
rz(-1.954345) q[1];
x q[2];
rz(2.5043082) q[3];
sx q[3];
rz(-2.6316959) q[3];
sx q[3];
rz(-0.84264681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4454322) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(0.0017722842) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(-1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6034265) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(1.4461393) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(-1.5015645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1743463) q[0];
sx q[0];
rz(-1.5513199) q[0];
sx q[0];
rz(2.0635701) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.026560606) q[2];
sx q[2];
rz(-1.5090669) q[2];
sx q[2];
rz(-0.82690566) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1135243) q[1];
sx q[1];
rz(-1.6599732) q[1];
sx q[1];
rz(2.813617) q[1];
x q[2];
rz(-2.1498508) q[3];
sx q[3];
rz(-2.1175044) q[3];
sx q[3];
rz(-2.9255097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.8224576) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(0.31931988) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(0.38326344) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4708913) q[0];
sx q[0];
rz(-2.0235217) q[0];
sx q[0];
rz(0.41505138) q[0];
x q[1];
rz(1.8750538) q[2];
sx q[2];
rz(-2.2396302) q[2];
sx q[2];
rz(2.5077016) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0723567) q[1];
sx q[1];
rz(-1.623739) q[1];
sx q[1];
rz(0.80979053) q[1];
rz(-pi) q[2];
rz(2.6978108) q[3];
sx q[3];
rz(-0.1212596) q[3];
sx q[3];
rz(-1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7982771) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.2822255) q[2];
rz(-1.6451689) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(1.0378029) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(-2.1077572) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4170096) q[0];
sx q[0];
rz(-1.7010744) q[0];
sx q[0];
rz(0.91086046) q[0];
rz(-pi) q[1];
rz(0.86976544) q[2];
sx q[2];
rz(-1.1681721) q[2];
sx q[2];
rz(1.0588156) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2227877) q[1];
sx q[1];
rz(-0.58090392) q[1];
sx q[1];
rz(-1.3851628) q[1];
rz(-0.71318993) q[3];
sx q[3];
rz(-1.7367559) q[3];
sx q[3];
rz(-2.1470269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(1.5283782) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6939659) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(1.7059965) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(0.96799093) q[2];
sx q[2];
rz(-1.5445166) q[2];
sx q[2];
rz(-1.9894285) q[2];
rz(1.3676436) q[3];
sx q[3];
rz(-2.805134) q[3];
sx q[3];
rz(0.72611879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
