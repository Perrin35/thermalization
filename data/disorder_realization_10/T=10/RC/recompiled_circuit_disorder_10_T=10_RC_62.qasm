OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4410285) q[0];
sx q[0];
rz(-1.1428042) q[0];
sx q[0];
rz(-1.2115275) q[0];
rz(-0.22663528) q[1];
sx q[1];
rz(-1.5770788) q[1];
sx q[1];
rz(-2.8432863) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5601215) q[0];
sx q[0];
rz(-1.5157962) q[0];
sx q[0];
rz(2.4762857) q[0];
rz(-pi) q[1];
rz(1.9036129) q[2];
sx q[2];
rz(-1.4822072) q[2];
sx q[2];
rz(2.7804136) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71516192) q[1];
sx q[1];
rz(-0.75725812) q[1];
sx q[1];
rz(0.84233474) q[1];
rz(1.467642) q[3];
sx q[3];
rz(-1.6543596) q[3];
sx q[3];
rz(0.32675693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6478708) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(0.99386627) q[2];
rz(0.99938756) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(-0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(-3.1233741) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(-2.6699064) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6149711) q[0];
sx q[0];
rz(-1.7978298) q[0];
sx q[0];
rz(1.660166) q[0];
x q[1];
rz(1.291044) q[2];
sx q[2];
rz(-1.0824167) q[2];
sx q[2];
rz(2.5807057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3968518) q[1];
sx q[1];
rz(-0.62780118) q[1];
sx q[1];
rz(0.46698924) q[1];
rz(1.3470115) q[3];
sx q[3];
rz(-2.5219678) q[3];
sx q[3];
rz(0.97298813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54962426) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(-0.96898752) q[2];
rz(-2.5668868) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(2.95978) q[0];
rz(2.0388942) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(-1.6859432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239379) q[0];
sx q[0];
rz(-1.0929937) q[0];
sx q[0];
rz(-0.91207232) q[0];
x q[1];
rz(-1.1481029) q[2];
sx q[2];
rz(-2.2659677) q[2];
sx q[2];
rz(2.9425651) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0297444) q[1];
sx q[1];
rz(-1.6964456) q[1];
sx q[1];
rz(1.5602342) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6772179) q[3];
sx q[3];
rz(-2.1933746) q[3];
sx q[3];
rz(2.4604083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2640947) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(-2.1739615) q[2];
rz(2.4140221) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85686344) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(2.5033584) q[0];
rz(-1.1278641) q[1];
sx q[1];
rz(-2.3141839) q[1];
sx q[1];
rz(1.9086054) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5467984) q[0];
sx q[0];
rz(-1.6398755) q[0];
sx q[0];
rz(2.6462206) q[0];
x q[1];
rz(-0.58926438) q[2];
sx q[2];
rz(-1.1009842) q[2];
sx q[2];
rz(-1.9643009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9945558) q[1];
sx q[1];
rz(-1.6068659) q[1];
sx q[1];
rz(-1.7493164) q[1];
rz(1.2992371) q[3];
sx q[3];
rz(-2.4304667) q[3];
sx q[3];
rz(1.3815051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2234852) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(2.3275862) q[2];
rz(-1.043184) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2994613) q[0];
sx q[0];
rz(-1.786754) q[0];
sx q[0];
rz(0.89170757) q[0];
rz(-1.2437598) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(2.9290501) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1067057) q[0];
sx q[0];
rz(-1.5505152) q[0];
sx q[0];
rz(-3.1253392) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5151305) q[2];
sx q[2];
rz(-2.5362483) q[2];
sx q[2];
rz(-0.84673131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5383496) q[1];
sx q[1];
rz(-1.0298567) q[1];
sx q[1];
rz(-2.4451838) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46194525) q[3];
sx q[3];
rz(-0.89832234) q[3];
sx q[3];
rz(0.40601054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1084958) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(2.6203716) q[2];
rz(-1.3850348) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(1.2683755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8591156) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(2.916472) q[0];
rz(-1.7865932) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(-0.37757847) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8416653) q[0];
sx q[0];
rz(-2.8746434) q[0];
sx q[0];
rz(-1.4873234) q[0];
x q[1];
rz(-2.9794681) q[2];
sx q[2];
rz(-2.8505278) q[2];
sx q[2];
rz(-1.7578917) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1212335) q[1];
sx q[1];
rz(-1.4660144) q[1];
sx q[1];
rz(-2.0038414) q[1];
x q[2];
rz(1.5402921) q[3];
sx q[3];
rz(-0.77611938) q[3];
sx q[3];
rz(3.0200849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(-2.6605576) q[2];
rz(-0.40361079) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0525381) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(2.9470434) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(-2.887168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5146778) q[0];
sx q[0];
rz(-1.8906381) q[0];
sx q[0];
rz(-0.11492782) q[0];
rz(-pi) q[1];
rz(1.7888072) q[2];
sx q[2];
rz(-0.62354747) q[2];
sx q[2];
rz(2.4965198) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2563045) q[1];
sx q[1];
rz(-1.2997775) q[1];
sx q[1];
rz(-2.6726252) q[1];
x q[2];
rz(-2.6894327) q[3];
sx q[3];
rz(-1.9981355) q[3];
sx q[3];
rz(2.6219581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1640132) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(1.3158201) q[2];
rz(2.2655462) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2106237) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(0.82398206) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(1.5751858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63021916) q[0];
sx q[0];
rz(-0.89238088) q[0];
sx q[0];
rz(1.1648965) q[0];
rz(2.1112061) q[2];
sx q[2];
rz(-1.4166797) q[2];
sx q[2];
rz(0.87481462) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6187001) q[1];
sx q[1];
rz(-2.61781) q[1];
sx q[1];
rz(2.3431542) q[1];
rz(-pi) q[2];
rz(-0.50623399) q[3];
sx q[3];
rz(-1.9444379) q[3];
sx q[3];
rz(0.95844275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2660797) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(-0.27754647) q[2];
rz(-1.4510441) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16185109) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.6850527) q[0];
rz(2.5121571) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(1.1368407) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0757383) q[0];
sx q[0];
rz(-2.542001) q[0];
sx q[0];
rz(2.0980741) q[0];
x q[1];
rz(-1.9929664) q[2];
sx q[2];
rz(-1.2264226) q[2];
sx q[2];
rz(-3.0712155) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6362308) q[1];
sx q[1];
rz(-0.33126918) q[1];
sx q[1];
rz(-1.2760217) q[1];
rz(-pi) q[2];
x q[2];
rz(1.198248) q[3];
sx q[3];
rz(-0.87009831) q[3];
sx q[3];
rz(-0.31853279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4920766) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(2.1006404) q[2];
rz(-3.1395636) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(-2.8295529) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390389) q[0];
sx q[0];
rz(-0.22452393) q[0];
sx q[0];
rz(0.94605207) q[0];
rz(-2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(-2.5295703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5123972) q[0];
sx q[0];
rz(-2.0266268) q[0];
sx q[0];
rz(0.7645316) q[0];
x q[1];
rz(2.7962748) q[2];
sx q[2];
rz(-1.3441663) q[2];
sx q[2];
rz(-2.741284) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7358688) q[1];
sx q[1];
rz(-1.4475665) q[1];
sx q[1];
rz(2.481639) q[1];
rz(1.1935812) q[3];
sx q[3];
rz(-1.0672788) q[3];
sx q[3];
rz(-0.52770381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0570021) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(0.65336147) q[2];
rz(2.7907794) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(-0.70070926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54031298) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(2.3868949) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(2.4738612) q[2];
sx q[2];
rz(-2.3420391) q[2];
sx q[2];
rz(2.5426368) q[2];
rz(0.084372088) q[3];
sx q[3];
rz(-1.2481239) q[3];
sx q[3];
rz(0.39831755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
