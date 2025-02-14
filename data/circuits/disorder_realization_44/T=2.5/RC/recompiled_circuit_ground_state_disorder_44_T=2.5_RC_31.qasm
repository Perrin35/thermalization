OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7437416) q[0];
sx q[0];
rz(-0.97111094) q[0];
sx q[0];
rz(2.4308393) q[0];
rz(0.96762586) q[1];
sx q[1];
rz(-3.1275446) q[1];
sx q[1];
rz(0.79298055) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6502888) q[0];
sx q[0];
rz(-0.86435917) q[0];
sx q[0];
rz(-0.44235787) q[0];
x q[1];
rz(-1.6184801) q[2];
sx q[2];
rz(-1.5565201) q[2];
sx q[2];
rz(-0.88231495) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0147018) q[1];
sx q[1];
rz(-1.1387991) q[1];
sx q[1];
rz(-0.72662455) q[1];
rz(-2.1737384) q[3];
sx q[3];
rz(-1.0558898) q[3];
sx q[3];
rz(3.0647354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8636318) q[2];
sx q[2];
rz(-0.42676723) q[2];
sx q[2];
rz(-0.58941907) q[2];
rz(-2.3943118) q[3];
sx q[3];
rz(-0.092844754) q[3];
sx q[3];
rz(2.5417627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079161949) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(-2.2404501) q[0];
rz(-2.1597553) q[1];
sx q[1];
rz(-2.5603309) q[1];
sx q[1];
rz(0.99220401) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3311148) q[0];
sx q[0];
rz(-1.7077291) q[0];
sx q[0];
rz(2.1934319) q[0];
rz(-pi) q[1];
rz(-1.7478701) q[2];
sx q[2];
rz(-2.3858527) q[2];
sx q[2];
rz(2.9301639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76490957) q[1];
sx q[1];
rz(-2.1071395) q[1];
sx q[1];
rz(0.77002828) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76073356) q[3];
sx q[3];
rz(-2.2983716) q[3];
sx q[3];
rz(-2.6618877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.09484) q[2];
sx q[2];
rz(-0.53043008) q[2];
sx q[2];
rz(3.1165822) q[2];
rz(-2.181634) q[3];
sx q[3];
rz(-0.046006087) q[3];
sx q[3];
rz(-0.068232603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18441021) q[0];
sx q[0];
rz(-0.010951696) q[0];
sx q[0];
rz(0.72407323) q[0];
rz(3.030576) q[1];
sx q[1];
rz(-2.5357775) q[1];
sx q[1];
rz(-0.012880005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6136821) q[0];
sx q[0];
rz(-2.9146195) q[0];
sx q[0];
rz(1.2140973) q[0];
x q[1];
rz(2.3518042) q[2];
sx q[2];
rz(-0.26726535) q[2];
sx q[2];
rz(2.546026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4531969) q[1];
sx q[1];
rz(-1.7368142) q[1];
sx q[1];
rz(-2.9012381) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88408459) q[3];
sx q[3];
rz(-1.2753162) q[3];
sx q[3];
rz(-1.0425764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24590242) q[2];
sx q[2];
rz(-0.7220214) q[2];
sx q[2];
rz(2.8172909) q[2];
rz(-2.6445828) q[3];
sx q[3];
rz(-0.23247601) q[3];
sx q[3];
rz(2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67320353) q[0];
sx q[0];
rz(-0.16981801) q[0];
sx q[0];
rz(0.47624269) q[0];
rz(-0.4854804) q[1];
sx q[1];
rz(-2.6464033) q[1];
sx q[1];
rz(-0.21864299) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5056498) q[0];
sx q[0];
rz(-2.6723862) q[0];
sx q[0];
rz(-2.0903265) q[0];
x q[1];
rz(-0.013979023) q[2];
sx q[2];
rz(-2.4179248) q[2];
sx q[2];
rz(0.47823634) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8805595) q[1];
sx q[1];
rz(-1.1487242) q[1];
sx q[1];
rz(-2.5690394) q[1];
rz(-pi) q[2];
rz(-2.5095482) q[3];
sx q[3];
rz(-0.38292745) q[3];
sx q[3];
rz(0.51198375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.789088) q[2];
sx q[2];
rz(-2.3829491) q[2];
sx q[2];
rz(0.57141203) q[2];
rz(0.93829751) q[3];
sx q[3];
rz(-0.58656991) q[3];
sx q[3];
rz(2.9686484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57662302) q[0];
sx q[0];
rz(-2.7374856) q[0];
sx q[0];
rz(0.47082666) q[0];
rz(-0.91104031) q[1];
sx q[1];
rz(-0.43993479) q[1];
sx q[1];
rz(-2.7104673) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3456988) q[0];
sx q[0];
rz(-0.78471334) q[0];
sx q[0];
rz(0.037952947) q[0];
x q[1];
rz(1.9091102) q[2];
sx q[2];
rz(-2.4801284) q[2];
sx q[2];
rz(-2.4544883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87157618) q[1];
sx q[1];
rz(-1.5983014) q[1];
sx q[1];
rz(-1.639495) q[1];
x q[2];
rz(-0.2016368) q[3];
sx q[3];
rz(-1.7072225) q[3];
sx q[3];
rz(-2.0755943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33991459) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(-2.8002296) q[2];
rz(-2.7692128) q[3];
sx q[3];
rz(-2.5058993) q[3];
sx q[3];
rz(0.46975964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.85207283) q[0];
sx q[0];
rz(-2.2884123) q[0];
sx q[0];
rz(0.12292718) q[0];
rz(0.3854824) q[1];
sx q[1];
rz(-2.3434134) q[1];
sx q[1];
rz(0.72365671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9117377) q[0];
sx q[0];
rz(-2.9617972) q[0];
sx q[0];
rz(1.8263962) q[0];
rz(-pi) q[1];
rz(-0.18823024) q[2];
sx q[2];
rz(-1.9005132) q[2];
sx q[2];
rz(-2.4356869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3299218) q[1];
sx q[1];
rz(-0.99024665) q[1];
sx q[1];
rz(0.73402053) q[1];
rz(-pi) q[2];
rz(-2.9272363) q[3];
sx q[3];
rz(-1.949614) q[3];
sx q[3];
rz(-0.54730584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3896997) q[2];
sx q[2];
rz(-0.60201001) q[2];
sx q[2];
rz(0.67328084) q[2];
rz(-0.26345396) q[3];
sx q[3];
rz(-0.47502381) q[3];
sx q[3];
rz(-2.8365005) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70169705) q[0];
sx q[0];
rz(-0.79653996) q[0];
sx q[0];
rz(2.6252966) q[0];
rz(-0.75597489) q[1];
sx q[1];
rz(-0.95840234) q[1];
sx q[1];
rz(0.36852536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7411503) q[0];
sx q[0];
rz(-1.0545122) q[0];
sx q[0];
rz(-0.80908605) q[0];
rz(-pi) q[1];
rz(-1.6096949) q[2];
sx q[2];
rz(-2.2752011) q[2];
sx q[2];
rz(-0.6608085) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.75223756) q[1];
sx q[1];
rz(-2.9417243) q[1];
sx q[1];
rz(1.0081069) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51492274) q[3];
sx q[3];
rz(-2.2201) q[3];
sx q[3];
rz(0.89735555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64510173) q[2];
sx q[2];
rz(-3.0848905) q[2];
sx q[2];
rz(0.48218316) q[2];
rz(0.21969806) q[3];
sx q[3];
rz(-2.4441661) q[3];
sx q[3];
rz(0.68035948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7242303) q[0];
sx q[0];
rz(-2.9927232) q[0];
sx q[0];
rz(2.505488) q[0];
rz(-0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(0.87463921) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3134523) q[0];
sx q[0];
rz(-2.9813672) q[0];
sx q[0];
rz(-1.942722) q[0];
rz(-pi) q[1];
rz(0.12107559) q[2];
sx q[2];
rz(-1.4619451) q[2];
sx q[2];
rz(-1.7055275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88747178) q[1];
sx q[1];
rz(-1.1855108) q[1];
sx q[1];
rz(-2.0068568) q[1];
x q[2];
rz(-0.88040027) q[3];
sx q[3];
rz(-2.1244308) q[3];
sx q[3];
rz(2.2601489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8735698) q[2];
sx q[2];
rz(-2.7251254) q[2];
sx q[2];
rz(-0.91656172) q[2];
rz(0.77740866) q[3];
sx q[3];
rz(-0.55792266) q[3];
sx q[3];
rz(2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1189608) q[0];
sx q[0];
rz(-2.4038778) q[0];
sx q[0];
rz(0.23271261) q[0];
rz(0.64838707) q[1];
sx q[1];
rz(-0.62260038) q[1];
sx q[1];
rz(-0.56786215) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85405871) q[0];
sx q[0];
rz(-1.2191448) q[0];
sx q[0];
rz(-1.132347) q[0];
rz(-pi) q[1];
rz(0.21505996) q[2];
sx q[2];
rz(-1.8369007) q[2];
sx q[2];
rz(1.4799445) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1901924) q[1];
sx q[1];
rz(-1.7170329) q[1];
sx q[1];
rz(1.7758572) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81553163) q[3];
sx q[3];
rz(-0.99247265) q[3];
sx q[3];
rz(-0.92507833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3719172) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(-0.5522716) q[2];
rz(-0.28198379) q[3];
sx q[3];
rz(-2.7115287) q[3];
sx q[3];
rz(-0.10274398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.929739) q[0];
sx q[0];
rz(-3.0775098) q[0];
sx q[0];
rz(-3.0098359) q[0];
rz(2.6051104) q[1];
sx q[1];
rz(-0.15833144) q[1];
sx q[1];
rz(-2.2333455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.344438) q[0];
sx q[0];
rz(-0.69509172) q[0];
sx q[0];
rz(-1.6174843) q[0];
rz(-pi) q[1];
rz(-0.26226173) q[2];
sx q[2];
rz(-1.1168257) q[2];
sx q[2];
rz(-1.5020467) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40501198) q[1];
sx q[1];
rz(-0.93376505) q[1];
sx q[1];
rz(-0.65959658) q[1];
rz(-pi) q[2];
rz(0.0018168505) q[3];
sx q[3];
rz(-1.7967136) q[3];
sx q[3];
rz(-0.29669138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.50916719) q[2];
sx q[2];
rz(-0.91370344) q[2];
sx q[2];
rz(-0.2779648) q[2];
rz(2.5748504) q[3];
sx q[3];
rz(-2.9394579) q[3];
sx q[3];
rz(-0.91293269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.81000281) q[0];
sx q[0];
rz(-1.7362052) q[0];
sx q[0];
rz(2.195634) q[0];
rz(-2.6592061) q[1];
sx q[1];
rz(-1.2336029) q[1];
sx q[1];
rz(-0.69419669) q[1];
rz(1.214477) q[2];
sx q[2];
rz(-2.5209485) q[2];
sx q[2];
rz(-2.4232558) q[2];
rz(-2.0770155) q[3];
sx q[3];
rz(-1.5983221) q[3];
sx q[3];
rz(1.7479001) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
