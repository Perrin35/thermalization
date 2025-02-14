OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.12282523) q[0];
sx q[0];
rz(-0.0039761877) q[0];
sx q[0];
rz(-3.0206326) q[0];
rz(-2.6236985) q[1];
sx q[1];
rz(-2.8105812) q[1];
sx q[1];
rz(-1.9537227) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8594583) q[0];
sx q[0];
rz(-1.4842502) q[0];
sx q[0];
rz(0.29973932) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14622525) q[2];
sx q[2];
rz(-1.0492965) q[2];
sx q[2];
rz(0.63298775) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8459699) q[1];
sx q[1];
rz(-0.69087183) q[1];
sx q[1];
rz(-0.9102896) q[1];
rz(-pi) q[2];
rz(0.066066001) q[3];
sx q[3];
rz(-2.5089536) q[3];
sx q[3];
rz(-3.0939915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0594222) q[2];
sx q[2];
rz(-2.1990081) q[2];
sx q[2];
rz(-1.4284632) q[2];
rz(1.7510471) q[3];
sx q[3];
rz(-2.3640552) q[3];
sx q[3];
rz(-0.20974222) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21935894) q[0];
sx q[0];
rz(-1.5809504) q[0];
sx q[0];
rz(-1.6844164) q[0];
rz(2.4600785) q[1];
sx q[1];
rz(-0.69029713) q[1];
sx q[1];
rz(0.95726475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5752561) q[0];
sx q[0];
rz(-2.0342376) q[0];
sx q[0];
rz(-1.8011679) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.886734) q[2];
sx q[2];
rz(-0.83410848) q[2];
sx q[2];
rz(-2.2535549) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87220472) q[1];
sx q[1];
rz(-1.1630668) q[1];
sx q[1];
rz(-3.0881626) q[1];
rz(-pi) q[2];
rz(-1.718347) q[3];
sx q[3];
rz(-3.0719764) q[3];
sx q[3];
rz(-0.92629647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7822781) q[2];
sx q[2];
rz(-1.0474397) q[2];
sx q[2];
rz(-2.2347343) q[2];
rz(2.4685229) q[3];
sx q[3];
rz(-0.38452092) q[3];
sx q[3];
rz(-1.8460021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68536711) q[0];
sx q[0];
rz(-2.4816368) q[0];
sx q[0];
rz(2.5947156) q[0];
rz(-1.7710255) q[1];
sx q[1];
rz(-1.1616881) q[1];
sx q[1];
rz(0.10017698) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91104394) q[0];
sx q[0];
rz(-1.3063138) q[0];
sx q[0];
rz(0.13767155) q[0];
rz(-pi) q[1];
rz(2.5008548) q[2];
sx q[2];
rz(-0.60534873) q[2];
sx q[2];
rz(-1.5126765) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0292058) q[1];
sx q[1];
rz(-0.76923623) q[1];
sx q[1];
rz(-0.33035128) q[1];
rz(1.1421575) q[3];
sx q[3];
rz(-1.996962) q[3];
sx q[3];
rz(-1.7249191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3935516) q[2];
sx q[2];
rz(-1.2057722) q[2];
sx q[2];
rz(1.0432517) q[2];
rz(0.68759632) q[3];
sx q[3];
rz(-2.6396773) q[3];
sx q[3];
rz(2.8858394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8461175) q[0];
sx q[0];
rz(-0.44545528) q[0];
sx q[0];
rz(-2.1050146) q[0];
rz(1.8702033) q[1];
sx q[1];
rz(-2.3174353) q[1];
sx q[1];
rz(-1.8545256) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66726738) q[0];
sx q[0];
rz(-2.0653353) q[0];
sx q[0];
rz(0.70566341) q[0];
x q[1];
rz(-2.4802568) q[2];
sx q[2];
rz(-1.2244389) q[2];
sx q[2];
rz(2.6544184) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8657664) q[1];
sx q[1];
rz(-0.681502) q[1];
sx q[1];
rz(-1.1033415) q[1];
rz(-pi) q[2];
rz(1.1274476) q[3];
sx q[3];
rz(-1.2928767) q[3];
sx q[3];
rz(-1.9390566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6048772) q[2];
sx q[2];
rz(-2.776919) q[2];
sx q[2];
rz(-0.019901179) q[2];
rz(-0.59602916) q[3];
sx q[3];
rz(-0.28087956) q[3];
sx q[3];
rz(-1.178044) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48536456) q[0];
sx q[0];
rz(-1.8940268) q[0];
sx q[0];
rz(1.1744936) q[0];
rz(2.8354722) q[1];
sx q[1];
rz(-1.0447634) q[1];
sx q[1];
rz(1.3694481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1533617) q[0];
sx q[0];
rz(-1.5643459) q[0];
sx q[0];
rz(3.1383443) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3197863) q[2];
sx q[2];
rz(-2.5282871) q[2];
sx q[2];
rz(-1.5410739) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28182733) q[1];
sx q[1];
rz(-1.9143811) q[1];
sx q[1];
rz(-2.9815893) q[1];
rz(-pi) q[2];
rz(2.4804875) q[3];
sx q[3];
rz(-2.5206158) q[3];
sx q[3];
rz(1.0266765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8451763) q[2];
sx q[2];
rz(-2.3974819) q[2];
sx q[2];
rz(-0.80294341) q[2];
rz(1.1154277) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(1.7723134) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95336103) q[0];
sx q[0];
rz(-0.53934923) q[0];
sx q[0];
rz(-3.0292335) q[0];
rz(0.27680963) q[1];
sx q[1];
rz(-1.1667629) q[1];
sx q[1];
rz(-2.4940431) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66747276) q[0];
sx q[0];
rz(-1.9074567) q[0];
sx q[0];
rz(-2.1483627) q[0];
rz(-pi) q[1];
rz(-0.66866691) q[2];
sx q[2];
rz(-1.200182) q[2];
sx q[2];
rz(-0.47793717) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4989483) q[1];
sx q[1];
rz(-1.1234087) q[1];
sx q[1];
rz(1.4397463) q[1];
rz(-pi) q[2];
rz(2.9660951) q[3];
sx q[3];
rz(-1.8030522) q[3];
sx q[3];
rz(-1.0009223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4623146) q[2];
sx q[2];
rz(-2.9913112) q[2];
sx q[2];
rz(-1.798299) q[2];
rz(0.51370931) q[3];
sx q[3];
rz(-1.6535583) q[3];
sx q[3];
rz(-2.546379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3419471) q[0];
sx q[0];
rz(-1.3222313) q[0];
sx q[0];
rz(-0.42022002) q[0];
rz(-1.7814024) q[1];
sx q[1];
rz(-1.9748634) q[1];
sx q[1];
rz(0.80418783) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59475454) q[0];
sx q[0];
rz(-1.9512984) q[0];
sx q[0];
rz(-2.1516933) q[0];
x q[1];
rz(1.2857513) q[2];
sx q[2];
rz(-0.80696304) q[2];
sx q[2];
rz(-0.80978823) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7090164) q[1];
sx q[1];
rz(-0.99919277) q[1];
sx q[1];
rz(2.9663621) q[1];
x q[2];
rz(3.0325417) q[3];
sx q[3];
rz(-1.4067459) q[3];
sx q[3];
rz(-1.4993639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5516025) q[2];
sx q[2];
rz(-1.0162105) q[2];
sx q[2];
rz(1.1691947) q[2];
rz(0.17360887) q[3];
sx q[3];
rz(-0.9089402) q[3];
sx q[3];
rz(1.199523) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6731897) q[0];
sx q[0];
rz(-2.115374) q[0];
sx q[0];
rz(1.9386559) q[0];
rz(2.9007593) q[1];
sx q[1];
rz(-0.92643654) q[1];
sx q[1];
rz(-0.9333207) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.653395) q[0];
sx q[0];
rz(-1.1943333) q[0];
sx q[0];
rz(-2.7791136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1740994) q[2];
sx q[2];
rz(-2.5317414) q[2];
sx q[2];
rz(1.4791453) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.62040239) q[1];
sx q[1];
rz(-2.7082084) q[1];
sx q[1];
rz(0.49232011) q[1];
rz(-pi) q[2];
rz(-1.0068137) q[3];
sx q[3];
rz(-1.4419334) q[3];
sx q[3];
rz(0.098524898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11747083) q[2];
sx q[2];
rz(-1.4535707) q[2];
sx q[2];
rz(3.0061099) q[2];
rz(-0.99902117) q[3];
sx q[3];
rz(-0.24805598) q[3];
sx q[3];
rz(-0.55618709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67021543) q[0];
sx q[0];
rz(-1.1286292) q[0];
sx q[0];
rz(1.7673329) q[0];
rz(-2.0595835) q[1];
sx q[1];
rz(-2.3550985) q[1];
sx q[1];
rz(-1.5793922) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.608273) q[0];
sx q[0];
rz(-2.1100304) q[0];
sx q[0];
rz(3.1031552) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57317971) q[2];
sx q[2];
rz(-0.64116353) q[2];
sx q[2];
rz(1.4631084) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97171451) q[1];
sx q[1];
rz(-2.5757667) q[1];
sx q[1];
rz(-0.91822718) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7855603) q[3];
sx q[3];
rz(-2.1904328) q[3];
sx q[3];
rz(2.5621264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9595327) q[2];
sx q[2];
rz(-1.5675194) q[2];
sx q[2];
rz(0.62425295) q[2];
rz(0.67638451) q[3];
sx q[3];
rz(-1.9703777) q[3];
sx q[3];
rz(-2.306365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71037173) q[0];
sx q[0];
rz(-1.8932764) q[0];
sx q[0];
rz(1.7927908) q[0];
rz(-0.30385941) q[1];
sx q[1];
rz(-1.6108578) q[1];
sx q[1];
rz(-2.2644728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5069339) q[0];
sx q[0];
rz(-1.7153694) q[0];
sx q[0];
rz(-2.4094482) q[0];
x q[1];
rz(1.2563324) q[2];
sx q[2];
rz(-0.34032492) q[2];
sx q[2];
rz(-0.2257502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.34906975) q[1];
sx q[1];
rz(-1.9510837) q[1];
sx q[1];
rz(-2.040634) q[1];
rz(-pi) q[2];
rz(0.56352625) q[3];
sx q[3];
rz(-1.0358264) q[3];
sx q[3];
rz(0.68171147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0554241) q[2];
sx q[2];
rz(-0.40579) q[2];
sx q[2];
rz(-2.1093192) q[2];
rz(-2.0530733) q[3];
sx q[3];
rz(-1.4884596) q[3];
sx q[3];
rz(0.77435875) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1879723) q[0];
sx q[0];
rz(-2.3379876) q[0];
sx q[0];
rz(0.048820989) q[0];
rz(1.8232518) q[1];
sx q[1];
rz(-1.7346458) q[1];
sx q[1];
rz(0.55191747) q[1];
rz(0.55808347) q[2];
sx q[2];
rz(-0.48102745) q[2];
sx q[2];
rz(-0.68651909) q[2];
rz(1.5416038) q[3];
sx q[3];
rz(-2.7053082) q[3];
sx q[3];
rz(-2.3655688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
