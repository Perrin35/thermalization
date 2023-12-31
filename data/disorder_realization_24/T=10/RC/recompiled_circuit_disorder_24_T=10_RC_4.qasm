OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(-1.0868602) q[0];
sx q[0];
rz(-1.342919) q[0];
rz(7.8328447) q[1];
sx q[1];
rz(3.0631493) q[1];
sx q[1];
rz(6.9258239) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064518236) q[0];
sx q[0];
rz(-1.0318349) q[0];
sx q[0];
rz(-1.5383188) q[0];
rz(-pi) q[1];
rz(-1.7625916) q[2];
sx q[2];
rz(-2.0711053) q[2];
sx q[2];
rz(1.2230011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3712284) q[1];
sx q[1];
rz(-0.1703573) q[1];
sx q[1];
rz(0.78070663) q[1];
rz(-pi) q[2];
rz(-2.9075165) q[3];
sx q[3];
rz(-1.5353068) q[3];
sx q[3];
rz(-1.8518098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47444433) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(-1.2791963) q[2];
rz(-0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(-0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779125) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(2.1610778) q[0];
rz(2.9837043) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(-0.78871361) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0448488) q[0];
sx q[0];
rz(-1.4474488) q[0];
sx q[0];
rz(-0.11476536) q[0];
rz(-pi) q[1];
rz(-1.6025701) q[2];
sx q[2];
rz(-1.8268158) q[2];
sx q[2];
rz(2.1928744) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84856725) q[1];
sx q[1];
rz(-1.2663406) q[1];
sx q[1];
rz(2.5679563) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7699645) q[3];
sx q[3];
rz(-1.0125481) q[3];
sx q[3];
rz(1.6127197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.366189) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(-0.186084) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9300951) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(-2.511456) q[0];
rz(-3.0139626) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(-2.4198467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8080224) q[0];
sx q[0];
rz(-2.434242) q[0];
sx q[0];
rz(2.7471514) q[0];
rz(1.1100936) q[2];
sx q[2];
rz(-2.2113872) q[2];
sx q[2];
rz(-5.0355807e-05) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7194781) q[1];
sx q[1];
rz(-2.0146703) q[1];
sx q[1];
rz(1.5281954) q[1];
x q[2];
rz(-2.1980522) q[3];
sx q[3];
rz(-1.2949847) q[3];
sx q[3];
rz(2.4081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3815986) q[2];
sx q[2];
rz(-0.95278946) q[2];
sx q[2];
rz(-0.90908137) q[2];
rz(2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58406126) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(-3.0990565) q[0];
rz(2.361239) q[1];
sx q[1];
rz(-0.50023729) q[1];
sx q[1];
rz(2.3775878) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9616868) q[0];
sx q[0];
rz(-0.67877239) q[0];
sx q[0];
rz(1.5508482) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4970384) q[2];
sx q[2];
rz(-1.4696215) q[2];
sx q[2];
rz(2.8855756) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9909775) q[1];
sx q[1];
rz(-2.5533278) q[1];
sx q[1];
rz(0.88612835) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4776811) q[3];
sx q[3];
rz(-1.9466562) q[3];
sx q[3];
rz(2.1223048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8445231) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(-0.53304535) q[2];
rz(-2.879203) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9168636) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(1.8977144) q[0];
rz(-0.22661701) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(-0.40333834) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6809083) q[0];
sx q[0];
rz(-1.7366689) q[0];
sx q[0];
rz(0.18780639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68571217) q[2];
sx q[2];
rz(-1.2843686) q[2];
sx q[2];
rz(-0.098766947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0163739) q[1];
sx q[1];
rz(-2.5902936) q[1];
sx q[1];
rz(-0.2972879) q[1];
rz(-pi) q[2];
rz(-0.67540695) q[3];
sx q[3];
rz(-1.1408148) q[3];
sx q[3];
rz(2.5436385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.32101813) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(-0.78249758) q[2];
rz(-2.0292422) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(-1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.430442) q[0];
sx q[0];
rz(-2.7395881) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(-0.09952155) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(-0.0064370357) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59505263) q[0];
sx q[0];
rz(-1.8227302) q[0];
sx q[0];
rz(-2.3810054) q[0];
rz(-pi) q[1];
rz(-1.1941031) q[2];
sx q[2];
rz(-1.6695392) q[2];
sx q[2];
rz(1.8284947) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41758075) q[1];
sx q[1];
rz(-1.9829826) q[1];
sx q[1];
rz(-2.8860303) q[1];
rz(-pi) q[2];
rz(-2.2925917) q[3];
sx q[3];
rz(-0.95522049) q[3];
sx q[3];
rz(-1.701223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36879888) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(-0.84645611) q[2];
rz(2.1438697) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0434175) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(-0.055710677) q[0];
rz(-0.78272351) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(-1.1605211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8502064) q[0];
sx q[0];
rz(-1.123748) q[0];
sx q[0];
rz(0.38711754) q[0];
rz(1.0497401) q[2];
sx q[2];
rz(-1.1546635) q[2];
sx q[2];
rz(0.001948826) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6351663) q[1];
sx q[1];
rz(-0.70235683) q[1];
sx q[1];
rz(-2.173645) q[1];
x q[2];
rz(-2.3485687) q[3];
sx q[3];
rz(-2.6245955) q[3];
sx q[3];
rz(1.5152064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(2.356142) q[2];
rz(0.75585946) q[3];
sx q[3];
rz(-1.2952341) q[3];
sx q[3];
rz(2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(-1.1428517) q[0];
rz(-1.8354592) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(2.725504) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04836719) q[0];
sx q[0];
rz(-2.0443633) q[0];
sx q[0];
rz(0.22305365) q[0];
x q[1];
rz(1.2942737) q[2];
sx q[2];
rz(-2.6575436) q[2];
sx q[2];
rz(2.7798142) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.63562288) q[1];
sx q[1];
rz(-1.1715679) q[1];
sx q[1];
rz(-1.3969621) q[1];
x q[2];
rz(1.3564524) q[3];
sx q[3];
rz(-0.48210258) q[3];
sx q[3];
rz(-1.6899504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1064421) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(0.32315928) q[2];
rz(2.9390826) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(0.37339661) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(1.1449822) q[0];
rz(1.1960944) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(2.6224565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0275092) q[0];
sx q[0];
rz(-2.7402096) q[0];
sx q[0];
rz(1.2252349) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22186188) q[2];
sx q[2];
rz(-1.8010745) q[2];
sx q[2];
rz(-0.78267539) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0874333) q[1];
sx q[1];
rz(-2.3490153) q[1];
sx q[1];
rz(1.4448318) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9037644) q[3];
sx q[3];
rz(-1.3999709) q[3];
sx q[3];
rz(1.2549787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3141994) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(3.1006151) q[2];
rz(-2.273902) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(-0.51122558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39559078) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(1.6145153) q[0];
rz(1.4279667) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(1.6428927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1536627) q[0];
sx q[0];
rz(-1.6462353) q[0];
sx q[0];
rz(-3.1213785) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0936436) q[2];
sx q[2];
rz(-1.7257294) q[2];
sx q[2];
rz(2.5028554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.26951075) q[1];
sx q[1];
rz(-1.2334358) q[1];
sx q[1];
rz(0.81706394) q[1];
x q[2];
rz(2.5999971) q[3];
sx q[3];
rz(-0.23532) q[3];
sx q[3];
rz(0.02863392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(2.995058) q[2];
rz(0.81418973) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9183337) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(-1.8756443) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(2.0064034) q[2];
sx q[2];
rz(-2.8786567) q[2];
sx q[2];
rz(-1.565956) q[2];
rz(2.2371348) q[3];
sx q[3];
rz(-0.96257985) q[3];
sx q[3];
rz(-1.0812159) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
