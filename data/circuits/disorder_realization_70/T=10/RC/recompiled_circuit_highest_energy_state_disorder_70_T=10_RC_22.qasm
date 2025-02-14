OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3938067) q[0];
sx q[0];
rz(-0.90016794) q[0];
sx q[0];
rz(1.7322487) q[0];
rz(-2.4318168) q[1];
sx q[1];
rz(-0.33325279) q[1];
sx q[1];
rz(2.7594653) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20724498) q[0];
sx q[0];
rz(-0.74688321) q[0];
sx q[0];
rz(2.9435959) q[0];
rz(-pi) q[1];
rz(-1.0898148) q[2];
sx q[2];
rz(-0.28533563) q[2];
sx q[2];
rz(1.7794045) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5643049) q[1];
sx q[1];
rz(-1.2049985) q[1];
sx q[1];
rz(-1.1759961) q[1];
rz(-0.11973937) q[3];
sx q[3];
rz(-0.88212195) q[3];
sx q[3];
rz(-2.7514091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0673151) q[2];
sx q[2];
rz(-1.5019608) q[2];
sx q[2];
rz(-0.074946694) q[2];
rz(1.2429169) q[3];
sx q[3];
rz(-0.26151812) q[3];
sx q[3];
rz(-0.3956795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13260929) q[0];
sx q[0];
rz(-0.41985303) q[0];
sx q[0];
rz(-2.5123151) q[0];
rz(-2.3043326) q[1];
sx q[1];
rz(-2.3910797) q[1];
sx q[1];
rz(0.57964051) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64855498) q[0];
sx q[0];
rz(-0.84024197) q[0];
sx q[0];
rz(1.296098) q[0];
x q[1];
rz(0.57168269) q[2];
sx q[2];
rz(-0.46077706) q[2];
sx q[2];
rz(-2.1850815) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3677153) q[1];
sx q[1];
rz(-1.9019777) q[1];
sx q[1];
rz(-0.41172387) q[1];
rz(-1.5509255) q[3];
sx q[3];
rz(-1.7345034) q[3];
sx q[3];
rz(0.57396561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.89977449) q[2];
sx q[2];
rz(-2.9782229) q[2];
sx q[2];
rz(-2.3972798) q[2];
rz(0.36738473) q[3];
sx q[3];
rz(-1.2920047) q[3];
sx q[3];
rz(1.7261837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71017569) q[0];
sx q[0];
rz(-0.95040584) q[0];
sx q[0];
rz(-1.0071734) q[0];
rz(-0.011064359) q[1];
sx q[1];
rz(-2.8304351) q[1];
sx q[1];
rz(-0.99753582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57118124) q[0];
sx q[0];
rz(-1.4332674) q[0];
sx q[0];
rz(-0.026022051) q[0];
rz(-pi) q[1];
rz(-0.80834016) q[2];
sx q[2];
rz(-1.488121) q[2];
sx q[2];
rz(0.17553628) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.026873206) q[1];
sx q[1];
rz(-2.623154) q[1];
sx q[1];
rz(-3.1264831) q[1];
rz(-pi) q[2];
rz(-2.6297731) q[3];
sx q[3];
rz(-1.5580172) q[3];
sx q[3];
rz(0.46841533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1204388) q[2];
sx q[2];
rz(-2.8595371) q[2];
sx q[2];
rz(0.8417449) q[2];
rz(2.4833931) q[3];
sx q[3];
rz(-0.87116146) q[3];
sx q[3];
rz(-0.021520821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4206674) q[0];
sx q[0];
rz(-0.1249211) q[0];
sx q[0];
rz(-0.7290054) q[0];
rz(2.3782102) q[1];
sx q[1];
rz(-0.63012505) q[1];
sx q[1];
rz(-0.36300945) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1141407) q[0];
sx q[0];
rz(-1.318256) q[0];
sx q[0];
rz(0.30014287) q[0];
x q[1];
rz(0.82421397) q[2];
sx q[2];
rz(-1.9938278) q[2];
sx q[2];
rz(-1.3049478) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2263086) q[1];
sx q[1];
rz(-1.7123509) q[1];
sx q[1];
rz(2.841921) q[1];
rz(-pi) q[2];
rz(-0.12673817) q[3];
sx q[3];
rz(-1.5180382) q[3];
sx q[3];
rz(-2.3873752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18569854) q[2];
sx q[2];
rz(-2.318435) q[2];
sx q[2];
rz(1.6419179) q[2];
rz(-0.58756346) q[3];
sx q[3];
rz(-2.1524119) q[3];
sx q[3];
rz(2.6232918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8374306) q[0];
sx q[0];
rz(-0.70697933) q[0];
sx q[0];
rz(0.28513232) q[0];
rz(-0.25310165) q[1];
sx q[1];
rz(-1.0292091) q[1];
sx q[1];
rz(2.0957799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8469197) q[0];
sx q[0];
rz(-0.6631279) q[0];
sx q[0];
rz(0.66594932) q[0];
rz(-pi) q[1];
rz(1.990647) q[2];
sx q[2];
rz(-1.1641181) q[2];
sx q[2];
rz(0.57957725) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95438984) q[1];
sx q[1];
rz(-2.7437401) q[1];
sx q[1];
rz(0.14601645) q[1];
rz(1.771403) q[3];
sx q[3];
rz(-0.90292519) q[3];
sx q[3];
rz(0.29800668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4046341) q[2];
sx q[2];
rz(-2.0817231) q[2];
sx q[2];
rz(0.57445478) q[2];
rz(-2.5308841) q[3];
sx q[3];
rz(-2.6210531) q[3];
sx q[3];
rz(-1.0569388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33559281) q[0];
sx q[0];
rz(-1.3845504) q[0];
sx q[0];
rz(-2.7912676) q[0];
rz(2.8490745) q[1];
sx q[1];
rz(-3.0151093) q[1];
sx q[1];
rz(0.77936053) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3199917) q[0];
sx q[0];
rz(-1.6896833) q[0];
sx q[0];
rz(3.0827808) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33213385) q[2];
sx q[2];
rz(-2.6382425) q[2];
sx q[2];
rz(2.533874) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36679572) q[1];
sx q[1];
rz(-1.2506984) q[1];
sx q[1];
rz(1.3572378) q[1];
x q[2];
rz(2.8661714) q[3];
sx q[3];
rz(-2.6044327) q[3];
sx q[3];
rz(2.8583683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.902035) q[2];
sx q[2];
rz(-1.6263447) q[2];
sx q[2];
rz(2.1774192) q[2];
rz(-0.17298175) q[3];
sx q[3];
rz(-1.0123342) q[3];
sx q[3];
rz(-0.13121901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59703374) q[0];
sx q[0];
rz(-1.9439161) q[0];
sx q[0];
rz(-0.069393754) q[0];
rz(-1.767905) q[1];
sx q[1];
rz(-1.8927788) q[1];
sx q[1];
rz(-2.6232041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98250472) q[0];
sx q[0];
rz(-0.90686528) q[0];
sx q[0];
rz(-0.21592617) q[0];
rz(-pi) q[1];
rz(1.0807627) q[2];
sx q[2];
rz(-1.9980717) q[2];
sx q[2];
rz(3.0025122) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9897108) q[1];
sx q[1];
rz(-0.1406142) q[1];
sx q[1];
rz(-0.18528823) q[1];
x q[2];
rz(-1.5126918) q[3];
sx q[3];
rz(-1.4062506) q[3];
sx q[3];
rz(-1.8376415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18846506) q[2];
sx q[2];
rz(-1.1981107) q[2];
sx q[2];
rz(2.8971064) q[2];
rz(1.2178347) q[3];
sx q[3];
rz(-3.0471424) q[3];
sx q[3];
rz(0.34734669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.090488) q[0];
sx q[0];
rz(-0.29260391) q[0];
sx q[0];
rz(-0.69945139) q[0];
rz(-0.72169101) q[1];
sx q[1];
rz(-1.7803918) q[1];
sx q[1];
rz(1.1915421) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22868294) q[0];
sx q[0];
rz(-1.7828724) q[0];
sx q[0];
rz(-1.6985083) q[0];
x q[1];
rz(0.6416816) q[2];
sx q[2];
rz(-0.26452082) q[2];
sx q[2];
rz(1.3542287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9758975) q[1];
sx q[1];
rz(-1.6553406) q[1];
sx q[1];
rz(-2.1102474) q[1];
x q[2];
rz(0.64719836) q[3];
sx q[3];
rz(-2.5997926) q[3];
sx q[3];
rz(1.7447646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9222074) q[2];
sx q[2];
rz(-2.3228513) q[2];
sx q[2];
rz(-1.4156263) q[2];
rz(2.7042232) q[3];
sx q[3];
rz(-0.24961095) q[3];
sx q[3];
rz(-2.6387446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.9472083) q[0];
sx q[0];
rz(-1.2057065) q[0];
sx q[0];
rz(-2.6956287) q[0];
rz(-1.3392316) q[1];
sx q[1];
rz(-0.65473348) q[1];
sx q[1];
rz(1.1599249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8561607) q[0];
sx q[0];
rz(-2.2158379) q[0];
sx q[0];
rz(1.1354394) q[0];
x q[1];
rz(1.2959506) q[2];
sx q[2];
rz(-0.27121085) q[2];
sx q[2];
rz(1.7618568) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7801108) q[1];
sx q[1];
rz(-1.0358397) q[1];
sx q[1];
rz(1.8423716) q[1];
rz(-pi) q[2];
rz(2.2723012) q[3];
sx q[3];
rz(-0.79330595) q[3];
sx q[3];
rz(-1.3877804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40598536) q[2];
sx q[2];
rz(-1.0286101) q[2];
sx q[2];
rz(-0.73060161) q[2];
rz(2.2405911) q[3];
sx q[3];
rz(-0.42947072) q[3];
sx q[3];
rz(-3.1072531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.63448298) q[0];
sx q[0];
rz(-2.6431838) q[0];
sx q[0];
rz(0.46257567) q[0];
rz(-2.0787461) q[1];
sx q[1];
rz(-1.7741508) q[1];
sx q[1];
rz(-3.0752693) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1076814) q[0];
sx q[0];
rz(-0.66921355) q[0];
sx q[0];
rz(0.50743703) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0424329) q[2];
sx q[2];
rz(-2.7688469) q[2];
sx q[2];
rz(-0.78102222) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5328153) q[1];
sx q[1];
rz(-2.3022396) q[1];
sx q[1];
rz(2.3904496) q[1];
rz(-pi) q[2];
rz(1.9948694) q[3];
sx q[3];
rz(-1.3641285) q[3];
sx q[3];
rz(0.28462946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.037420951) q[2];
sx q[2];
rz(-0.53356844) q[2];
sx q[2];
rz(-1.9083692) q[2];
rz(2.6836266) q[3];
sx q[3];
rz(-0.27126867) q[3];
sx q[3];
rz(-2.7428194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11252277) q[0];
sx q[0];
rz(-1.4334913) q[0];
sx q[0];
rz(1.7423472) q[0];
rz(0.15432547) q[1];
sx q[1];
rz(-1.7185153) q[1];
sx q[1];
rz(-1.1850866) q[1];
rz(2.828601) q[2];
sx q[2];
rz(-1.9421158) q[2];
sx q[2];
rz(-2.4696642) q[2];
rz(-2.6575487) q[3];
sx q[3];
rz(-2.0917907) q[3];
sx q[3];
rz(1.2863822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
