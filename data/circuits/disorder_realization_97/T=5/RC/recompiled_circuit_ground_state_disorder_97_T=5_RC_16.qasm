OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.360541) q[0];
sx q[0];
rz(-0.6147576) q[0];
sx q[0];
rz(2.0701261) q[0];
rz(-0.40845025) q[1];
sx q[1];
rz(4.0790494) q[1];
sx q[1];
rz(11.117878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0385742) q[0];
sx q[0];
rz(-1.9644613) q[0];
sx q[0];
rz(-1.0746687) q[0];
x q[1];
rz(-1.9117457) q[2];
sx q[2];
rz(-1.023479) q[2];
sx q[2];
rz(-2.7173619) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71943563) q[1];
sx q[1];
rz(-1.7195175) q[1];
sx q[1];
rz(1.4805111) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0759645) q[3];
sx q[3];
rz(-0.77407661) q[3];
sx q[3];
rz(0.89895338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94413269) q[2];
sx q[2];
rz(-2.2977273) q[2];
sx q[2];
rz(2.3485363) q[2];
rz(2.872725) q[3];
sx q[3];
rz(-1.3258508) q[3];
sx q[3];
rz(2.1813006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.31084138) q[0];
sx q[0];
rz(-0.38250592) q[0];
sx q[0];
rz(-2.261396) q[0];
rz(2.510732) q[1];
sx q[1];
rz(-2.3289101) q[1];
sx q[1];
rz(2.0426483) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55148516) q[0];
sx q[0];
rz(-1.5657358) q[0];
sx q[0];
rz(0.17781114) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3395202) q[2];
sx q[2];
rz(-2.510609) q[2];
sx q[2];
rz(-0.6501261) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4447282) q[1];
sx q[1];
rz(-0.46485365) q[1];
sx q[1];
rz(-1.9397199) q[1];
rz(-0.46862015) q[3];
sx q[3];
rz(-2.1543401) q[3];
sx q[3];
rz(-2.0035721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7546996) q[2];
sx q[2];
rz(-1.4852445) q[2];
sx q[2];
rz(2.5999542) q[2];
rz(-2.4382639) q[3];
sx q[3];
rz(-2.7024305) q[3];
sx q[3];
rz(-2.7642803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2331053) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(0.096906699) q[0];
rz(-0.58755177) q[1];
sx q[1];
rz(-2.2511626) q[1];
sx q[1];
rz(2.5403835) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158949) q[0];
sx q[0];
rz(-0.89381274) q[0];
sx q[0];
rz(1.8034794) q[0];
rz(0.16232441) q[2];
sx q[2];
rz(-1.5078203) q[2];
sx q[2];
rz(-1.8520825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2612348) q[1];
sx q[1];
rz(-2.6812892) q[1];
sx q[1];
rz(-0.017437915) q[1];
x q[2];
rz(-2.2279068) q[3];
sx q[3];
rz(-1.8538215) q[3];
sx q[3];
rz(1.0237657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2931557) q[2];
sx q[2];
rz(-1.4508672) q[2];
sx q[2];
rz(-2.4270774) q[2];
rz(-1.4826639) q[3];
sx q[3];
rz(-0.27479333) q[3];
sx q[3];
rz(0.37884918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191384) q[0];
sx q[0];
rz(-1.9062573) q[0];
sx q[0];
rz(1.6374913) q[0];
rz(-2.0488886) q[1];
sx q[1];
rz(-1.8605109) q[1];
sx q[1];
rz(2.5561996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4248237) q[0];
sx q[0];
rz(-0.94858525) q[0];
sx q[0];
rz(-1.9068002) q[0];
rz(-pi) q[1];
rz(-0.010527654) q[2];
sx q[2];
rz(-1.6750355) q[2];
sx q[2];
rz(-2.5679905) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5820739) q[1];
sx q[1];
rz(-2.2041956) q[1];
sx q[1];
rz(0.60057171) q[1];
x q[2];
rz(2.500072) q[3];
sx q[3];
rz(-1.975394) q[3];
sx q[3];
rz(1.4610964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66712159) q[2];
sx q[2];
rz(-0.67255628) q[2];
sx q[2];
rz(-2.7453864) q[2];
rz(0.19043663) q[3];
sx q[3];
rz(-2.2021144) q[3];
sx q[3];
rz(2.6207391) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1162972) q[0];
sx q[0];
rz(-2.681356) q[0];
sx q[0];
rz(-2.7468371) q[0];
rz(1.0074298) q[1];
sx q[1];
rz(-1.9837244) q[1];
sx q[1];
rz(1.5347068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3956063) q[0];
sx q[0];
rz(-1.2885546) q[0];
sx q[0];
rz(2.372118) q[0];
rz(-pi) q[1];
rz(-1.9362279) q[2];
sx q[2];
rz(-2.2534568) q[2];
sx q[2];
rz(-2.3548369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.04567178) q[1];
sx q[1];
rz(-1.878698) q[1];
sx q[1];
rz(-2.9363757) q[1];
rz(-2.0232361) q[3];
sx q[3];
rz(-1.1471738) q[3];
sx q[3];
rz(0.37003368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69636238) q[2];
sx q[2];
rz(-2.4424489) q[2];
sx q[2];
rz(2.9635079) q[2];
rz(0.94545025) q[3];
sx q[3];
rz(-0.33792308) q[3];
sx q[3];
rz(0.43615714) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6984542) q[0];
sx q[0];
rz(-0.66265023) q[0];
sx q[0];
rz(0.15289256) q[0];
rz(-1.0180417) q[1];
sx q[1];
rz(-0.94296229) q[1];
sx q[1];
rz(0.61000383) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40326443) q[0];
sx q[0];
rz(-0.77304196) q[0];
sx q[0];
rz(-0.18456642) q[0];
x q[1];
rz(2.1446682) q[2];
sx q[2];
rz(-1.5639493) q[2];
sx q[2];
rz(-2.8459542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.568012) q[1];
sx q[1];
rz(-1.3369155) q[1];
sx q[1];
rz(-2.3877445) q[1];
rz(-pi) q[2];
rz(-1.6527324) q[3];
sx q[3];
rz(-1.473236) q[3];
sx q[3];
rz(-0.20906258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6792004) q[2];
sx q[2];
rz(-1.8588763) q[2];
sx q[2];
rz(2.882615) q[2];
rz(0.14608832) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(2.5025388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4884969) q[0];
sx q[0];
rz(-2.275832) q[0];
sx q[0];
rz(-0.75188941) q[0];
rz(-2.409626) q[1];
sx q[1];
rz(-1.7611793) q[1];
sx q[1];
rz(-2.6123349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7649496) q[0];
sx q[0];
rz(-1.8762454) q[0];
sx q[0];
rz(-0.16965387) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13649444) q[2];
sx q[2];
rz(-1.1195099) q[2];
sx q[2];
rz(0.1642483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.098245278) q[1];
sx q[1];
rz(-1.6520513) q[1];
sx q[1];
rz(3.1286865) q[1];
rz(-pi) q[2];
rz(-2.1629754) q[3];
sx q[3];
rz(-0.53272034) q[3];
sx q[3];
rz(2.5902093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25972128) q[2];
sx q[2];
rz(-1.0292116) q[2];
sx q[2];
rz(-2.3465346) q[2];
rz(1.9184387) q[3];
sx q[3];
rz(-1.7984248) q[3];
sx q[3];
rz(-0.77643001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935598) q[0];
sx q[0];
rz(-1.5482276) q[0];
sx q[0];
rz(3.026631) q[0];
rz(-0.089275442) q[1];
sx q[1];
rz(-2.142579) q[1];
sx q[1];
rz(2.9866536) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45231014) q[0];
sx q[0];
rz(-1.4703106) q[0];
sx q[0];
rz(-1.4402251) q[0];
rz(-2.5903715) q[2];
sx q[2];
rz(-2.1586426) q[2];
sx q[2];
rz(1.8648704) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7748002) q[1];
sx q[1];
rz(-1.3143365) q[1];
sx q[1];
rz(-0.40089295) q[1];
x q[2];
rz(2.2959034) q[3];
sx q[3];
rz(-1.5165303) q[3];
sx q[3];
rz(2.3783663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0577008) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(2.4532301) q[2];
rz(-2.5734731) q[3];
sx q[3];
rz(-1.0024242) q[3];
sx q[3];
rz(2.4585371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532918) q[0];
sx q[0];
rz(-0.81128565) q[0];
sx q[0];
rz(1.2055093) q[0];
rz(1.0909117) q[1];
sx q[1];
rz(-0.58735192) q[1];
sx q[1];
rz(1.5376512) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0552222) q[0];
sx q[0];
rz(-1.1366664) q[0];
sx q[0];
rz(-0.40121292) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1411016) q[2];
sx q[2];
rz(-2.1386991) q[2];
sx q[2];
rz(2.9141324) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9133105) q[1];
sx q[1];
rz(-0.78332114) q[1];
sx q[1];
rz(-1.4128774) q[1];
x q[2];
rz(2.5439369) q[3];
sx q[3];
rz(-1.171805) q[3];
sx q[3];
rz(1.9982709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9566112) q[2];
sx q[2];
rz(-1.1143755) q[2];
sx q[2];
rz(1.1448917) q[2];
rz(1.4780809) q[3];
sx q[3];
rz(-0.76668113) q[3];
sx q[3];
rz(-1.112282) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0803273) q[0];
sx q[0];
rz(-0.63617951) q[0];
sx q[0];
rz(1.435745) q[0];
rz(-2.0044633) q[1];
sx q[1];
rz(-2.4500193) q[1];
sx q[1];
rz(0.93708509) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6244753) q[0];
sx q[0];
rz(-1.269425) q[0];
sx q[0];
rz(-1.6624438) q[0];
rz(-1.996973) q[2];
sx q[2];
rz(-2.0405031) q[2];
sx q[2];
rz(-1.5091648) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91793375) q[1];
sx q[1];
rz(-1.5964701) q[1];
sx q[1];
rz(0.77083807) q[1];
rz(-pi) q[2];
rz(0.61599515) q[3];
sx q[3];
rz(-1.2912036) q[3];
sx q[3];
rz(-0.071690138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3252141) q[2];
sx q[2];
rz(-0.64322317) q[2];
sx q[2];
rz(-0.31036672) q[2];
rz(0.54272932) q[3];
sx q[3];
rz(-2.2118745) q[3];
sx q[3];
rz(-0.31699666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.19001374) q[0];
sx q[0];
rz(-2.0484476) q[0];
sx q[0];
rz(0.53566326) q[0];
rz(-0.71470064) q[1];
sx q[1];
rz(-1.8938046) q[1];
sx q[1];
rz(-1.6751777) q[1];
rz(1.1874448) q[2];
sx q[2];
rz(-1.9786096) q[2];
sx q[2];
rz(-0.65828029) q[2];
rz(-1.6313995) q[3];
sx q[3];
rz(-0.84280673) q[3];
sx q[3];
rz(-1.6731586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
