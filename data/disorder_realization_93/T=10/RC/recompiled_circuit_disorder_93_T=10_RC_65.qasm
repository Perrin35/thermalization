OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(-1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(-1.7927875) q[1];
sx q[1];
rz(-0.92372149) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1169352) q[0];
sx q[0];
rz(-0.57514656) q[0];
sx q[0];
rz(2.2206109) q[0];
rz(0.14416868) q[2];
sx q[2];
rz(-1.8509794) q[2];
sx q[2];
rz(0.35137128) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87286283) q[1];
sx q[1];
rz(-1.9951207) q[1];
sx q[1];
rz(-0.55117589) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6472595) q[3];
sx q[3];
rz(-1.9088073) q[3];
sx q[3];
rz(0.61275834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2279921) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(-2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0682003) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(-1.4555567) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992422) q[0];
sx q[0];
rz(-0.99827168) q[0];
sx q[0];
rz(0.19897977) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7669719) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(-2.3945216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96869722) q[1];
sx q[1];
rz(-0.52297938) q[1];
sx q[1];
rz(1.5437267) q[1];
rz(-pi) q[2];
rz(0.060828408) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(0.5371679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7130647) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(-2.9591566) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(-0.3119719) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75333726) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(-1.172539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4370255) q[0];
sx q[0];
rz(-0.66172681) q[0];
sx q[0];
rz(-2.6252803) q[0];
x q[1];
rz(1.2433979) q[2];
sx q[2];
rz(-2.2990169) q[2];
sx q[2];
rz(-0.74795216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6781569) q[1];
sx q[1];
rz(-2.0883745) q[1];
sx q[1];
rz(-0.28097681) q[1];
x q[2];
rz(0.058733744) q[3];
sx q[3];
rz(-1.3204832) q[3];
sx q[3];
rz(0.55610031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(0.91119901) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38055414) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(0.52350837) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9273705) q[0];
sx q[0];
rz(-2.4282051) q[0];
sx q[0];
rz(0.58332304) q[0];
x q[1];
rz(1.4857616) q[2];
sx q[2];
rz(-1.8997314) q[2];
sx q[2];
rz(2.0563682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6701811) q[1];
sx q[1];
rz(-1.3172611) q[1];
sx q[1];
rz(2.2815435) q[1];
x q[2];
rz(2.162096) q[3];
sx q[3];
rz(-1.4569836) q[3];
sx q[3];
rz(2.1554961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(2.5775487) q[2];
rz(-2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(-2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.4915285) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-0.99194828) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6742453) q[0];
sx q[0];
rz(-1.6773946) q[0];
sx q[0];
rz(0.17515134) q[0];
rz(3.132658) q[2];
sx q[2];
rz(-1.2679456) q[2];
sx q[2];
rz(-0.47700275) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3494898) q[1];
sx q[1];
rz(-1.2621242) q[1];
sx q[1];
rz(-1.5285138) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12985142) q[3];
sx q[3];
rz(-2.3240945) q[3];
sx q[3];
rz(1.070147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-2.8707855) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(0.38189608) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(-1.7165002) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3544918) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(-1.8255193) q[0];
rz(1.233333) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(0.5459107) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9565935) q[1];
sx q[1];
rz(-1.6374267) q[1];
sx q[1];
rz(2.0208298) q[1];
rz(-pi) q[2];
rz(0.53374966) q[3];
sx q[3];
rz(-1.0373877) q[3];
sx q[3];
rz(2.651754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(-2.9610736) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5086223) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(-1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-0.82180506) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9200631) q[0];
sx q[0];
rz(-2.1863345) q[0];
sx q[0];
rz(-2.2277742) q[0];
rz(-pi) q[1];
rz(0.84341151) q[2];
sx q[2];
rz(-1.9551829) q[2];
sx q[2];
rz(-2.9316528) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.070548363) q[1];
sx q[1];
rz(-2.3729635) q[1];
sx q[1];
rz(-0.32960906) q[1];
rz(-0.52845593) q[3];
sx q[3];
rz(-0.65675694) q[3];
sx q[3];
rz(2.9627851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3372779) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(2.896893) q[2];
rz(3.0120567) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(2.3186671) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(1.8364505) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37874052) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(0.61600323) q[0];
x q[1];
rz(-0.70456409) q[2];
sx q[2];
rz(-1.2837871) q[2];
sx q[2];
rz(1.2003843) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23755632) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(-1.1035641) q[1];
x q[2];
rz(3.1180624) q[3];
sx q[3];
rz(-1.2315893) q[3];
sx q[3];
rz(-0.48285218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(1.6513599) q[2];
rz(-2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(2.4386491) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212595) q[0];
sx q[0];
rz(-0.54176211) q[0];
sx q[0];
rz(2.3938177) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1019194) q[2];
sx q[2];
rz(-1.0618321) q[2];
sx q[2];
rz(2.2850125) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96303899) q[1];
sx q[1];
rz(-0.56561618) q[1];
sx q[1];
rz(1.8815243) q[1];
rz(-2.8076257) q[3];
sx q[3];
rz(-2.1595862) q[3];
sx q[3];
rz(0.86563084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(-1.3396324) q[2];
rz(-2.83589) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(-1.3302749) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16383485) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(-1.9158069) q[0];
rz(0.90351358) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9194473) q[0];
sx q[0];
rz(-1.7788017) q[0];
sx q[0];
rz(-2.749445) q[0];
x q[1];
rz(1.1852766) q[2];
sx q[2];
rz(-1.2670994) q[2];
sx q[2];
rz(0.38244837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2965282) q[1];
sx q[1];
rz(-1.5025286) q[1];
sx q[1];
rz(-1.8370085) q[1];
rz(-pi) q[2];
rz(-1.1115132) q[3];
sx q[3];
rz(-1.1772403) q[3];
sx q[3];
rz(0.72898385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6580711) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(-1.1432077) q[2];
rz(-3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(1.6121929) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(-0.62190965) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-0.8016349) q[2];
sx q[2];
rz(-2.2710706) q[2];
sx q[2];
rz(2.0588277) q[2];
rz(2.2223496) q[3];
sx q[3];
rz(-1.5083434) q[3];
sx q[3];
rz(-1.2283243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
