OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(-2.3119976) q[0];
sx q[0];
rz(-0.15396804) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(0.33831236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4727288) q[0];
sx q[0];
rz(-2.7358486) q[0];
sx q[0];
rz(2.2292482) q[0];
rz(1.1233653) q[2];
sx q[2];
rz(-0.72927232) q[2];
sx q[2];
rz(-2.4016618) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.073804341) q[1];
sx q[1];
rz(-1.5177625) q[1];
sx q[1];
rz(2.8552613) q[1];
rz(1.0584352) q[3];
sx q[3];
rz(-1.5570939) q[3];
sx q[3];
rz(2.0126359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.9677229) q[2];
rz(0.075803444) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(-3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(-0.064963438) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(-1.2423135) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9984098) q[0];
sx q[0];
rz(-1.8185203) q[0];
sx q[0];
rz(-1.4383016) q[0];
x q[1];
rz(-0.1588891) q[2];
sx q[2];
rz(-2.1974568) q[2];
sx q[2];
rz(-1.5140669) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9638483) q[1];
sx q[1];
rz(-1.4114393) q[1];
sx q[1];
rz(-2.6049155) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70988016) q[3];
sx q[3];
rz(-1.1871561) q[3];
sx q[3];
rz(2.3185454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80766455) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(0.57404533) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42049256) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(0.72845355) q[0];
rz(-1.6473673) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(-1.0167936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3166312) q[0];
sx q[0];
rz(-1.6065856) q[0];
sx q[0];
rz(-0.081598452) q[0];
rz(2.0605893) q[2];
sx q[2];
rz(-2.042633) q[2];
sx q[2];
rz(-2.8420198) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6779855) q[1];
sx q[1];
rz(-1.4336168) q[1];
sx q[1];
rz(2.3198818) q[1];
rz(1.8524283) q[3];
sx q[3];
rz(-1.5236119) q[3];
sx q[3];
rz(-0.55164528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(1.0495079) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(-1.1857741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8124354) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(-1.6695492) q[0];
rz(2.4064348) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(0.24681117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4287764) q[0];
sx q[0];
rz(-1.719559) q[0];
sx q[0];
rz(2.2674198) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39101379) q[2];
sx q[2];
rz(-2.6303929) q[2];
sx q[2];
rz(-0.15399394) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66407953) q[1];
sx q[1];
rz(-1.941136) q[1];
sx q[1];
rz(1.627229) q[1];
x q[2];
rz(2.187192) q[3];
sx q[3];
rz(-1.8064926) q[3];
sx q[3];
rz(2.1637722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4776769) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(-1.5412615) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(-2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108869) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(-1.8141618) q[0];
rz(-1.56303) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(0.24838233) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43863338) q[0];
sx q[0];
rz(-0.54134936) q[0];
sx q[0];
rz(-3.0754509) q[0];
rz(-2.9551198) q[2];
sx q[2];
rz(-1.4154134) q[2];
sx q[2];
rz(1.965167) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76368139) q[1];
sx q[1];
rz(-1.3994201) q[1];
sx q[1];
rz(-0.59314368) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31110839) q[3];
sx q[3];
rz(-0.6725544) q[3];
sx q[3];
rz(-0.5141408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43626943) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(-0.79745897) q[2];
rz(-0.38875368) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(-0.50271547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.1335063) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(-1.3457993) q[0];
rz(-1.0812409) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(-3.016901) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4437618) q[0];
sx q[0];
rz(-0.19653453) q[0];
sx q[0];
rz(2.6854808) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2776676) q[2];
sx q[2];
rz(-2.1595402) q[2];
sx q[2];
rz(-0.53520203) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2248968) q[1];
sx q[1];
rz(-1.5142913) q[1];
sx q[1];
rz(-2.4042261) q[1];
rz(-pi) q[2];
rz(-1.7894621) q[3];
sx q[3];
rz(-1.339774) q[3];
sx q[3];
rz(0.71393379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51320118) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(-0.61895269) q[2];
rz(-2.0882873) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(-2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597647) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(-1.0429617) q[0];
rz(-0.46328059) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(1.0707062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28744222) q[0];
sx q[0];
rz(-1.5713912) q[0];
sx q[0];
rz(2.6575412) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4378909) q[2];
sx q[2];
rz(-1.9153567) q[2];
sx q[2];
rz(0.29495707) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13094014) q[1];
sx q[1];
rz(-0.7796692) q[1];
sx q[1];
rz(0.14203771) q[1];
x q[2];
rz(0.33220746) q[3];
sx q[3];
rz(-1.5515965) q[3];
sx q[3];
rz(1.7393877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36859194) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(0.32361844) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(-1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535646) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.9352242) q[0];
rz(-1.9288829) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(-2.1941197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3377209) q[0];
sx q[0];
rz(-2.5944355) q[0];
sx q[0];
rz(-1.7659811) q[0];
rz(-0.73080365) q[2];
sx q[2];
rz(-2.9060504) q[2];
sx q[2];
rz(1.3256324) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35093388) q[1];
sx q[1];
rz(-0.62218636) q[1];
sx q[1];
rz(-1.15637) q[1];
rz(-2.3341228) q[3];
sx q[3];
rz(-1.315457) q[3];
sx q[3];
rz(1.0367928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5802713) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(2.6718111) q[2];
rz(-1.3011159) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(-0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(-1.3289733) q[0];
rz(2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-2.9387617) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.348939) q[0];
sx q[0];
rz(-1.5447504) q[0];
sx q[0];
rz(0.039020122) q[0];
x q[1];
rz(-1.3615666) q[2];
sx q[2];
rz(-1.3522569) q[2];
sx q[2];
rz(2.2184559) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1111787) q[1];
sx q[1];
rz(-1.6200388) q[1];
sx q[1];
rz(-1.4700252) q[1];
x q[2];
rz(1.7248475) q[3];
sx q[3];
rz(-2.027958) q[3];
sx q[3];
rz(-0.16897133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(-0.99651304) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(-0.80741185) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080169454) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(0.043047992) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(0.28082401) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94175324) q[0];
sx q[0];
rz(-2.0240677) q[0];
sx q[0];
rz(-2.0487294) q[0];
rz(-pi) q[1];
rz(2.4091987) q[2];
sx q[2];
rz(-2.8576982) q[2];
sx q[2];
rz(0.64901272) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3481969) q[1];
sx q[1];
rz(-1.6284202) q[1];
sx q[1];
rz(-3.1053931) q[1];
x q[2];
rz(1.2383575) q[3];
sx q[3];
rz(-1.3849764) q[3];
sx q[3];
rz(-0.18248724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3165555) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(-0.62310702) q[2];
rz(2.1394219) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(-2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(-0.80550823) q[2];
sx q[2];
rz(-1.1136354) q[2];
sx q[2];
rz(-1.8352933) q[2];
rz(-1.1602041) q[3];
sx q[3];
rz(-1.2413597) q[3];
sx q[3];
rz(1.6551457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
