OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39785102) q[0];
sx q[0];
rz(4.1127036) q[0];
sx q[0];
rz(10.135531) q[0];
rz(0.96762586) q[1];
sx q[1];
rz(-3.1275446) q[1];
sx q[1];
rz(0.79298055) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7638057) q[0];
sx q[0];
rz(-1.9024693) q[0];
sx q[0];
rz(2.3275231) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5231126) q[2];
sx q[2];
rz(-1.5850726) q[2];
sx q[2];
rz(2.2592777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88411623) q[1];
sx q[1];
rz(-0.82472825) q[1];
sx q[1];
rz(-0.60666879) q[1];
x q[2];
rz(0.6019117) q[3];
sx q[3];
rz(-1.0546393) q[3];
sx q[3];
rz(-1.3207787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8636318) q[2];
sx q[2];
rz(-0.42676723) q[2];
sx q[2];
rz(-2.5521736) q[2];
rz(-2.3943118) q[3];
sx q[3];
rz(-3.0487479) q[3];
sx q[3];
rz(-2.5417627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0624307) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(2.2404501) q[0];
rz(-2.1597553) q[1];
sx q[1];
rz(-2.5603309) q[1];
sx q[1];
rz(0.99220401) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81047786) q[0];
sx q[0];
rz(-1.4338636) q[0];
sx q[0];
rz(-0.94816072) q[0];
x q[1];
rz(-1.7478701) q[2];
sx q[2];
rz(-0.75573993) q[2];
sx q[2];
rz(-2.9301639) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7957669) q[1];
sx q[1];
rz(-2.2123033) q[1];
sx q[1];
rz(0.87916763) q[1];
rz(-pi) q[2];
rz(0.68282376) q[3];
sx q[3];
rz(-1.0299333) q[3];
sx q[3];
rz(1.4861432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.09484) q[2];
sx q[2];
rz(-0.53043008) q[2];
sx q[2];
rz(3.1165822) q[2];
rz(-2.181634) q[3];
sx q[3];
rz(-3.0955866) q[3];
sx q[3];
rz(-3.0733601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18441021) q[0];
sx q[0];
rz(-0.010951696) q[0];
sx q[0];
rz(-2.4175194) q[0];
rz(-3.030576) q[1];
sx q[1];
rz(-2.5357775) q[1];
sx q[1];
rz(0.012880005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52791053) q[0];
sx q[0];
rz(-0.22697313) q[0];
sx q[0];
rz(1.2140973) q[0];
rz(-1.7628646) q[2];
sx q[2];
rz(-1.3837866) q[2];
sx q[2];
rz(2.929304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47564214) q[1];
sx q[1];
rz(-2.8503909) q[1];
sx q[1];
rz(2.5282574) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7666088) q[3];
sx q[3];
rz(-0.91911784) q[3];
sx q[3];
rz(-2.3790272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8956902) q[2];
sx q[2];
rz(-0.7220214) q[2];
sx q[2];
rz(0.32430172) q[2];
rz(0.4970099) q[3];
sx q[3];
rz(-0.23247601) q[3];
sx q[3];
rz(2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683891) q[0];
sx q[0];
rz(-2.9717746) q[0];
sx q[0];
rz(-0.47624269) q[0];
rz(0.4854804) q[1];
sx q[1];
rz(-0.49518934) q[1];
sx q[1];
rz(2.9229497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46313686) q[0];
sx q[0];
rz(-1.3443724) q[0];
sx q[0];
rz(1.1562267) q[0];
rz(-pi) q[1];
rz(3.1276136) q[2];
sx q[2];
rz(-0.7236679) q[2];
sx q[2];
rz(-0.47823634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.051603457) q[1];
sx q[1];
rz(-2.0877503) q[1];
sx q[1];
rz(-1.0801143) q[1];
x q[2];
rz(0.31422887) q[3];
sx q[3];
rz(-1.7933729) q[3];
sx q[3];
rz(-1.4861388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.789088) q[2];
sx q[2];
rz(-0.75864351) q[2];
sx q[2];
rz(0.57141203) q[2];
rz(-0.93829751) q[3];
sx q[3];
rz(-0.58656991) q[3];
sx q[3];
rz(-2.9686484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5649696) q[0];
sx q[0];
rz(-2.7374856) q[0];
sx q[0];
rz(0.47082666) q[0];
rz(-2.2305523) q[1];
sx q[1];
rz(-2.7016579) q[1];
sx q[1];
rz(-2.7104673) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79589383) q[0];
sx q[0];
rz(-2.3568793) q[0];
sx q[0];
rz(3.1036397) q[0];
rz(-pi) q[1];
rz(1.9091102) q[2];
sx q[2];
rz(-0.66146427) q[2];
sx q[2];
rz(2.4544883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.70111245) q[1];
sx q[1];
rz(-1.639469) q[1];
sx q[1];
rz(-3.1140226) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7100077) q[3];
sx q[3];
rz(-1.7705342) q[3];
sx q[3];
rz(-2.6090007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8016781) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(0.3413631) q[2];
rz(-0.3723799) q[3];
sx q[3];
rz(-2.5058993) q[3];
sx q[3];
rz(-0.46975964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2895198) q[0];
sx q[0];
rz(-0.85318035) q[0];
sx q[0];
rz(3.0186655) q[0];
rz(-0.3854824) q[1];
sx q[1];
rz(-2.3434134) q[1];
sx q[1];
rz(2.4179359) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.052304) q[0];
sx q[0];
rz(-1.6160242) q[0];
sx q[0];
rz(-1.3967229) q[0];
rz(-pi) q[1];
rz(1.9060109) q[2];
sx q[2];
rz(-1.7487757) q[2];
sx q[2];
rz(2.2151057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.21863114) q[1];
sx q[1];
rz(-2.1653163) q[1];
sx q[1];
rz(-2.2943952) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9272363) q[3];
sx q[3];
rz(-1.1919787) q[3];
sx q[3];
rz(0.54730584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75189292) q[2];
sx q[2];
rz(-0.60201001) q[2];
sx q[2];
rz(-2.4683118) q[2];
rz(0.26345396) q[3];
sx q[3];
rz(-0.47502381) q[3];
sx q[3];
rz(2.8365005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4398956) q[0];
sx q[0];
rz(-0.79653996) q[0];
sx q[0];
rz(-0.51629603) q[0];
rz(2.3856178) q[1];
sx q[1];
rz(-2.1831903) q[1];
sx q[1];
rz(2.7730673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4004423) q[0];
sx q[0];
rz(-2.0870805) q[0];
sx q[0];
rz(2.3325066) q[0];
rz(1.6096949) q[2];
sx q[2];
rz(-2.2752011) q[2];
sx q[2];
rz(-2.4807841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8766899) q[1];
sx q[1];
rz(-1.4646834) q[1];
sx q[1];
rz(-1.4011032) q[1];
rz(-pi) q[2];
rz(-2.1462899) q[3];
sx q[3];
rz(-0.80484521) q[3];
sx q[3];
rz(0.14508776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4964909) q[2];
sx q[2];
rz(-3.0848905) q[2];
sx q[2];
rz(-2.6594095) q[2];
rz(-2.9218946) q[3];
sx q[3];
rz(-2.4441661) q[3];
sx q[3];
rz(-2.4612332) q[3];
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
rz(-2.7242303) q[0];
sx q[0];
rz(-2.9927232) q[0];
sx q[0];
rz(-2.505488) q[0];
rz(0.96027374) q[1];
sx q[1];
rz(-2.2033913) q[1];
sx q[1];
rz(0.87463921) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8281404) q[0];
sx q[0];
rz(-0.16022542) q[0];
sx q[0];
rz(-1.942722) q[0];
rz(-pi) q[1];
rz(1.6804439) q[2];
sx q[2];
rz(-1.4504408) q[2];
sx q[2];
rz(-2.9936444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1369774) q[1];
sx q[1];
rz(-2.5680313) q[1];
sx q[1];
rz(0.80570813) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88040027) q[3];
sx q[3];
rz(-2.1244308) q[3];
sx q[3];
rz(-0.88144377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26802289) q[2];
sx q[2];
rz(-2.7251254) q[2];
sx q[2];
rz(-2.2250309) q[2];
rz(-2.364184) q[3];
sx q[3];
rz(-0.55792266) q[3];
sx q[3];
rz(2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022631835) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(-2.90888) q[0];
rz(2.4932056) q[1];
sx q[1];
rz(-0.62260038) q[1];
sx q[1];
rz(0.56786215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85405871) q[0];
sx q[0];
rz(-1.9224478) q[0];
sx q[0];
rz(-1.132347) q[0];
x q[1];
rz(-2.9265327) q[2];
sx q[2];
rz(-1.304692) q[2];
sx q[2];
rz(-1.4799445) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.1901924) q[1];
sx q[1];
rz(-1.4245598) q[1];
sx q[1];
rz(-1.3657354) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81553163) q[3];
sx q[3];
rz(-2.14912) q[3];
sx q[3];
rz(2.2165143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3719172) q[2];
sx q[2];
rz(-0.19352517) q[2];
sx q[2];
rz(-0.5522716) q[2];
rz(-0.28198379) q[3];
sx q[3];
rz(-2.7115287) q[3];
sx q[3];
rz(3.0388487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2118537) q[0];
sx q[0];
rz(-0.064082853) q[0];
sx q[0];
rz(0.1317568) q[0];
rz(0.53648221) q[1];
sx q[1];
rz(-2.9832612) q[1];
sx q[1];
rz(-2.2333455) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73777481) q[0];
sx q[0];
rz(-1.5409011) q[0];
sx q[0];
rz(0.87624082) q[0];
rz(-2.0591878) q[2];
sx q[2];
rz(-0.51967144) q[2];
sx q[2];
rz(0.95275098) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51146347) q[1];
sx q[1];
rz(-0.88246934) q[1];
sx q[1];
rz(0.87911112) q[1];
rz(-1.3448787) q[3];
sx q[3];
rz(-1.5690256) q[3];
sx q[3];
rz(1.2745119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.50916719) q[2];
sx q[2];
rz(-0.91370344) q[2];
sx q[2];
rz(0.2779648) q[2];
rz(2.5748504) q[3];
sx q[3];
rz(-2.9394579) q[3];
sx q[3];
rz(2.22866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3315898) q[0];
sx q[0];
rz(-1.7362052) q[0];
sx q[0];
rz(2.195634) q[0];
rz(0.48238659) q[1];
sx q[1];
rz(-1.2336029) q[1];
sx q[1];
rz(-0.69419669) q[1];
rz(-2.8972068) q[2];
sx q[2];
rz(-2.1472211) q[2];
sx q[2];
rz(0.28923464) q[2];
rz(0.031470555) q[3];
sx q[3];
rz(-2.0768055) q[3];
sx q[3];
rz(0.16184645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
