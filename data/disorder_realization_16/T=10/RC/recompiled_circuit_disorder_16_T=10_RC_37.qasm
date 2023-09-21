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
rz(0.83377588) q[1];
sx q[1];
rz(-2.1492465) q[1];
sx q[1];
rz(-0.33831236) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1725537) q[0];
sx q[0];
rz(-1.8882897) q[0];
sx q[0];
rz(-2.88455) q[0];
x q[1];
rz(0.89262427) q[2];
sx q[2];
rz(-1.2783588) q[2];
sx q[2];
rz(-0.48722789) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.073804341) q[1];
sx q[1];
rz(-1.5177625) q[1];
sx q[1];
rz(-0.28633134) q[1];
x q[2];
rz(-2.0831574) q[3];
sx q[3];
rz(-1.5570939) q[3];
sx q[3];
rz(-1.1289568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.9677229) q[2];
rz(-3.0657892) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(3.0766292) q[0];
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
rz(0.64105469) q[0];
sx q[0];
rz(-0.28028742) q[0];
sx q[0];
rz(-2.6602402) q[0];
rz(-1.785948) q[2];
sx q[2];
rz(-0.64385507) q[2];
sx q[2];
rz(1.3607963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8426659) q[1];
sx q[1];
rz(-2.0999523) q[1];
sx q[1];
rz(1.3859205) q[1];
rz(0.70988016) q[3];
sx q[3];
rz(-1.1871561) q[3];
sx q[3];
rz(0.82304728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3339281) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(2.5675473) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(-3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7211001) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(-2.4131391) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(2.1247991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66660488) q[0];
sx q[0];
rz(-0.089086108) q[0];
sx q[0];
rz(-2.7276917) q[0];
rz(-pi) q[1];
rz(2.0605893) q[2];
sx q[2];
rz(-1.0989597) q[2];
sx q[2];
rz(2.8420198) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46360717) q[1];
sx q[1];
rz(-1.4336168) q[1];
sx q[1];
rz(0.82171085) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4025027) q[3];
sx q[3];
rz(-0.28545359) q[3];
sx q[3];
rz(-1.1807549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(-0.63878757) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(-1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291572) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(-1.6695492) q[0];
rz(-2.4064348) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(2.8947815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4287764) q[0];
sx q[0];
rz(-1.719559) q[0];
sx q[0];
rz(-2.2674198) q[0];
rz(-0.39101379) q[2];
sx q[2];
rz(-2.6303929) q[2];
sx q[2];
rz(2.9875987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88627316) q[1];
sx q[1];
rz(-1.5181932) q[1];
sx q[1];
rz(-2.7707151) q[1];
rz(-pi) q[2];
rz(-0.28624268) q[3];
sx q[3];
rz(-2.1677368) q[3];
sx q[3];
rz(0.7569353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4776769) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.5412615) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-2.0440846) q[3];
sx q[3];
rz(-0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8108869) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(-1.3274308) q[0];
rz(-1.56303) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(0.24838233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9527234) q[0];
sx q[0];
rz(-1.5367322) q[0];
sx q[0];
rz(2.6012095) q[0];
rz(-pi) q[1];
rz(-1.4127172) q[2];
sx q[2];
rz(-1.7549967) q[2];
sx q[2];
rz(0.36518156) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0866962) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(2.8413248) q[1];
rz(0.31110839) q[3];
sx q[3];
rz(-2.4690383) q[3];
sx q[3];
rz(0.5141408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7053232) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(2.3441337) q[2];
rz(-2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-0.50271547) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(-1.7957934) q[0];
rz(2.0603518) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(0.1246917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69783083) q[0];
sx q[0];
rz(-2.9450581) q[0];
sx q[0];
rz(-0.4561119) q[0];
rz(-pi) q[1];
rz(-1.8639251) q[2];
sx q[2];
rz(-0.98205245) q[2];
sx q[2];
rz(-2.6063906) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4253937) q[1];
sx q[1];
rz(-0.73912207) q[1];
sx q[1];
rz(-3.0576586) q[1];
rz(-0.74495875) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(1.6572286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6283915) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(2.0882873) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(-0.9296023) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597647) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(0.46328059) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(1.0707062) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844855) q[0];
sx q[0];
rz(-2.6575408) q[0];
sx q[0];
rz(0.0012782106) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7878739) q[2];
sx q[2];
rz(-2.7732447) q[2];
sx q[2];
rz(3.0596717) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0740944) q[1];
sx q[1];
rz(-0.80103445) q[1];
sx q[1];
rz(1.4317516) q[1];
rz(-3.0827818) q[3];
sx q[3];
rz(-2.8088514) q[3];
sx q[3];
rz(2.9174093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36859194) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(-2.1598024) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.6535646) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(1.9288829) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(-0.94747296) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5652126) q[0];
sx q[0];
rz(-1.0351666) q[0];
sx q[0];
rz(-0.11760786) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73080365) q[2];
sx q[2];
rz(-2.9060504) q[2];
sx q[2];
rz(-1.3256324) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35093388) q[1];
sx q[1];
rz(-0.62218636) q[1];
sx q[1];
rz(1.15637) q[1];
rz(-pi) q[2];
rz(-2.3341228) q[3];
sx q[3];
rz(-1.8261357) q[3];
sx q[3];
rz(-1.0367928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(2.6718111) q[2];
rz(-1.8404768) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(2.8619213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8895421) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(2.3503616) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(-0.20283094) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77712599) q[0];
sx q[0];
rz(-1.6098032) q[0];
sx q[0];
rz(-1.5447306) q[0];
rz(-pi) q[1];
rz(-2.9183396) q[2];
sx q[2];
rz(-1.3666144) q[2];
sx q[2];
rz(2.5399361) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.030414) q[1];
sx q[1];
rz(-1.5215538) q[1];
sx q[1];
rz(-1.4700252) q[1];
rz(-2.6796954) q[3];
sx q[3];
rz(-1.4326722) q[3];
sx q[3];
rz(-1.671333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7245076) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(2.1450796) q[2];
rz(0.35342446) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(-0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0614232) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(-2.9598575) q[0];
rz(3.0985447) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(0.28082401) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94175324) q[0];
sx q[0];
rz(-1.117525) q[0];
sx q[0];
rz(-1.0928632) q[0];
rz(-pi) q[1];
rz(-0.73239399) q[2];
sx q[2];
rz(-2.8576982) q[2];
sx q[2];
rz(-2.4925799) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3481969) q[1];
sx q[1];
rz(-1.6284202) q[1];
sx q[1];
rz(-3.1053931) q[1];
rz(-0.19631581) q[3];
sx q[3];
rz(-1.8972978) q[3];
sx q[3];
rz(-1.6895837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3165555) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(-2.5184856) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(-2.5785057) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6476718) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(-0.87396809) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-2.3360844) q[2];
sx q[2];
rz(-2.0279573) q[2];
sx q[2];
rz(1.3062994) q[2];
rz(0.35691805) q[3];
sx q[3];
rz(-1.9580943) q[3];
sx q[3];
rz(3.0860268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
