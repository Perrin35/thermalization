OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37501332) q[0];
sx q[0];
rz(4.3809173) q[0];
sx q[0];
rz(8.3349174) q[0];
rz(-2.0774948) q[1];
sx q[1];
rz(-1.2531333) q[1];
sx q[1];
rz(-0.90322948) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922674) q[0];
sx q[0];
rz(-2.196133) q[0];
sx q[0];
rz(-2.3614285) q[0];
x q[1];
rz(3.1268812) q[2];
sx q[2];
rz(-2.2539469) q[2];
sx q[2];
rz(0.25757441) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.049202327) q[1];
sx q[1];
rz(-2.4396439) q[1];
sx q[1];
rz(-0.54767056) q[1];
rz(3.0211307) q[3];
sx q[3];
rz(-2.0673923) q[3];
sx q[3];
rz(0.42460006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4452867) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(-0.43329263) q[2];
rz(-2.8516234) q[3];
sx q[3];
rz(-2.2765997) q[3];
sx q[3];
rz(0.32052952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9669773) q[0];
sx q[0];
rz(-2.2791635) q[0];
sx q[0];
rz(-1.2906661) q[0];
rz(-2.4621452) q[1];
sx q[1];
rz(-2.3652855) q[1];
sx q[1];
rz(0.89250934) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36223012) q[0];
sx q[0];
rz(-2.9625872) q[0];
sx q[0];
rz(1.711861) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9694445) q[2];
sx q[2];
rz(-2.1855693) q[2];
sx q[2];
rz(-3.1190256) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.635517) q[1];
sx q[1];
rz(-0.92787023) q[1];
sx q[1];
rz(-1.9159813) q[1];
rz(2.8909182) q[3];
sx q[3];
rz(-1.2077577) q[3];
sx q[3];
rz(2.9533197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14900011) q[2];
sx q[2];
rz(-1.0789824) q[2];
sx q[2];
rz(1.1326257) q[2];
rz(0.66631404) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(-1.9168436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3742974) q[0];
sx q[0];
rz(-0.39300028) q[0];
sx q[0];
rz(-0.98180109) q[0];
rz(-0.94217316) q[1];
sx q[1];
rz(-1.3092382) q[1];
sx q[1];
rz(0.91032496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3142598) q[0];
sx q[0];
rz(-1.1480486) q[0];
sx q[0];
rz(-2.6092922) q[0];
x q[1];
rz(-1.4508171) q[2];
sx q[2];
rz(-1.2960805) q[2];
sx q[2];
rz(-0.071823013) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1012464) q[1];
sx q[1];
rz(-1.9034428) q[1];
sx q[1];
rz(-0.31942792) q[1];
rz(-pi) q[2];
rz(2.0098088) q[3];
sx q[3];
rz(-1.3512423) q[3];
sx q[3];
rz(-1.1846402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2008449) q[2];
sx q[2];
rz(-0.36483279) q[2];
sx q[2];
rz(0.18690898) q[2];
rz(2.9777891) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(0.12266172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0201037) q[0];
sx q[0];
rz(-1.4153471) q[0];
sx q[0];
rz(2.4439268) q[0];
rz(2.5949219) q[1];
sx q[1];
rz(-0.29269871) q[1];
sx q[1];
rz(-1.1950511) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23604017) q[0];
sx q[0];
rz(-2.4396606) q[0];
sx q[0];
rz(-2.1208956) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9505492) q[2];
sx q[2];
rz(-1.4378387) q[2];
sx q[2];
rz(1.3933448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4711485) q[1];
sx q[1];
rz(-0.19757825) q[1];
sx q[1];
rz(-1.0001282) q[1];
rz(-pi) q[2];
rz(1.2863808) q[3];
sx q[3];
rz(-1.5002325) q[3];
sx q[3];
rz(-2.0224935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4635072) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(2.9441693) q[2];
rz(0.10489634) q[3];
sx q[3];
rz(-2.0320804) q[3];
sx q[3];
rz(-3.0535898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5941493) q[0];
sx q[0];
rz(-0.53590411) q[0];
sx q[0];
rz(2.1323668) q[0];
rz(-0.028845305) q[1];
sx q[1];
rz(-1.4776769) q[1];
sx q[1];
rz(-0.65863329) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.665425) q[0];
sx q[0];
rz(-1.1575677) q[0];
sx q[0];
rz(-0.40979235) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0844857) q[2];
sx q[2];
rz(-1.3895814) q[2];
sx q[2];
rz(2.0985556) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.027623) q[1];
sx q[1];
rz(-1.7391035) q[1];
sx q[1];
rz(-1.6644003) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8010718) q[3];
sx q[3];
rz(-2.252823) q[3];
sx q[3];
rz(1.3547699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7625526) q[2];
sx q[2];
rz(-0.54658824) q[2];
sx q[2];
rz(-0.23400447) q[2];
rz(-1.0903357) q[3];
sx q[3];
rz(-1.5749616) q[3];
sx q[3];
rz(2.0719297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2284018) q[0];
sx q[0];
rz(-2.985432) q[0];
sx q[0];
rz(1.2983904) q[0];
rz(0.78650728) q[1];
sx q[1];
rz(-2.2857917) q[1];
sx q[1];
rz(1.7826805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0010064) q[0];
sx q[0];
rz(-1.6504262) q[0];
sx q[0];
rz(2.4796159) q[0];
rz(-pi) q[1];
rz(2.4491389) q[2];
sx q[2];
rz(-0.34231774) q[2];
sx q[2];
rz(-2.3794425) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1383776) q[1];
sx q[1];
rz(-1.5691461) q[1];
sx q[1];
rz(-1.5695851) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2511061) q[3];
sx q[3];
rz(-2.0123008) q[3];
sx q[3];
rz(2.6755345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7274373) q[2];
sx q[2];
rz(-1.6748019) q[2];
sx q[2];
rz(2.9690913) q[2];
rz(-2.6532459) q[3];
sx q[3];
rz(-1.9988873) q[3];
sx q[3];
rz(1.1138227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1682424) q[0];
sx q[0];
rz(-2.2689447) q[0];
sx q[0];
rz(-2.8016222) q[0];
rz(-2.8151457) q[1];
sx q[1];
rz(-2.4531334) q[1];
sx q[1];
rz(1.2350157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1025196) q[0];
sx q[0];
rz(-0.93329859) q[0];
sx q[0];
rz(-1.1028642) q[0];
x q[1];
rz(0.45845647) q[2];
sx q[2];
rz(-1.2453015) q[2];
sx q[2];
rz(-2.1874962) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3307747) q[1];
sx q[1];
rz(-2.0502809) q[1];
sx q[1];
rz(-1.564784) q[1];
x q[2];
rz(-2.6458287) q[3];
sx q[3];
rz(-0.44124441) q[3];
sx q[3];
rz(-0.10291162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1116703) q[2];
sx q[2];
rz(-2.3387574) q[2];
sx q[2];
rz(-2.5862582) q[2];
rz(-0.16767821) q[3];
sx q[3];
rz(-1.6043112) q[3];
sx q[3];
rz(2.663747) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70917201) q[0];
sx q[0];
rz(-2.2116311) q[0];
sx q[0];
rz(1.1093371) q[0];
rz(-1.1322016) q[1];
sx q[1];
rz(-0.39674509) q[1];
sx q[1];
rz(0.94165492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1718194) q[0];
sx q[0];
rz(-1.0189153) q[0];
sx q[0];
rz(0.18751796) q[0];
x q[1];
rz(1.0815762) q[2];
sx q[2];
rz(-2.8205546) q[2];
sx q[2];
rz(-0.58656091) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4577427) q[1];
sx q[1];
rz(-2.2096429) q[1];
sx q[1];
rz(-0.15775494) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38576084) q[3];
sx q[3];
rz(-1.3732135) q[3];
sx q[3];
rz(2.945874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7534916) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(1.6027742) q[2];
rz(-1.7704891) q[3];
sx q[3];
rz(-1.1857827) q[3];
sx q[3];
rz(-0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0592773) q[0];
sx q[0];
rz(-2.5264854) q[0];
sx q[0];
rz(-0.96784651) q[0];
rz(-1.2713426) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(-0.06180067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58690155) q[0];
sx q[0];
rz(-2.1036316) q[0];
sx q[0];
rz(-2.0237049) q[0];
rz(-pi) q[1];
rz(-0.65431739) q[2];
sx q[2];
rz(-2.8192602) q[2];
sx q[2];
rz(-2.2155025) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1872431) q[1];
sx q[1];
rz(-2.4496485) q[1];
sx q[1];
rz(2.991597) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4260208) q[3];
sx q[3];
rz(-2.712783) q[3];
sx q[3];
rz(2.9675067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9515848) q[2];
sx q[2];
rz(-2.9374359) q[2];
sx q[2];
rz(3.0657213) q[2];
rz(1.9742981) q[3];
sx q[3];
rz(-1.1018402) q[3];
sx q[3];
rz(-0.30437881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9645914) q[0];
sx q[0];
rz(-2.8622506) q[0];
sx q[0];
rz(0.28468537) q[0];
rz(3.0668861) q[1];
sx q[1];
rz(-1.2295405) q[1];
sx q[1];
rz(1.7237192) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95686326) q[0];
sx q[0];
rz(-1.7635582) q[0];
sx q[0];
rz(-0.75988976) q[0];
x q[1];
rz(2.1277587) q[2];
sx q[2];
rz(-1.6979453) q[2];
sx q[2];
rz(-0.47455041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1571341) q[1];
sx q[1];
rz(-0.95370953) q[1];
sx q[1];
rz(-1.0544712) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5261493) q[3];
sx q[3];
rz(-2.4952125) q[3];
sx q[3];
rz(0.97585362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3384) q[2];
sx q[2];
rz(-2.000838) q[2];
sx q[2];
rz(1.2791862) q[2];
rz(0.98215669) q[3];
sx q[3];
rz(-1.4582062) q[3];
sx q[3];
rz(-0.88472432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90047705) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(-2.1263532) q[1];
sx q[1];
rz(-1.992234) q[1];
sx q[1];
rz(2.0043859) q[1];
rz(-0.10812689) q[2];
sx q[2];
rz(-1.3095289) q[2];
sx q[2];
rz(1.1332569) q[2];
rz(-0.99526631) q[3];
sx q[3];
rz(-0.91777663) q[3];
sx q[3];
rz(-1.4025626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
