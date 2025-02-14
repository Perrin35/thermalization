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
rz(-0.32831353) q[0];
sx q[0];
rz(-1.1348731) q[0];
sx q[0];
rz(-2.5757134) q[0];
rz(-2.8973051) q[1];
sx q[1];
rz(-1.7841508) q[1];
sx q[1];
rz(-0.53425962) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8847647) q[0];
sx q[0];
rz(-2.3951911) q[0];
sx q[0];
rz(-1.166626) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5466629) q[2];
sx q[2];
rz(-1.317136) q[2];
sx q[2];
rz(-0.67753917) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.104887) q[1];
sx q[1];
rz(-1.5944696) q[1];
sx q[1];
rz(2.8410683) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47655063) q[3];
sx q[3];
rz(-0.50560564) q[3];
sx q[3];
rz(1.4534392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1285706) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(2.8924083) q[2];
rz(-1.7976044) q[3];
sx q[3];
rz(-1.598571) q[3];
sx q[3];
rz(-0.070076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6473815) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(0.98980728) q[0];
rz(0.79065943) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(-0.31165037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30325952) q[0];
sx q[0];
rz(-1.0619535) q[0];
sx q[0];
rz(-0.80755393) q[0];
rz(-3.0333854) q[2];
sx q[2];
rz(-1.3851067) q[2];
sx q[2];
rz(2.3600887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6271882) q[1];
sx q[1];
rz(-1.7429105) q[1];
sx q[1];
rz(-1.3283511) q[1];
x q[2];
rz(2.8007224) q[3];
sx q[3];
rz(-2.1580527) q[3];
sx q[3];
rz(2.1275008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0011757294) q[2];
sx q[2];
rz(-1.9813462) q[2];
sx q[2];
rz(2.1902093) q[2];
rz(-2.9133255) q[3];
sx q[3];
rz(-2.164866) q[3];
sx q[3];
rz(1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839639) q[0];
sx q[0];
rz(-1.5887337) q[0];
sx q[0];
rz(-3.0271295) q[0];
rz(-1.0022256) q[1];
sx q[1];
rz(-2.7148235) q[1];
sx q[1];
rz(-0.1618298) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9660258) q[0];
sx q[0];
rz(-1.9593666) q[0];
sx q[0];
rz(-2.2184664) q[0];
x q[1];
rz(0.21185565) q[2];
sx q[2];
rz(-0.39688928) q[2];
sx q[2];
rz(-0.8241764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.075167478) q[1];
sx q[1];
rz(-0.6614162) q[1];
sx q[1];
rz(0.39638806) q[1];
rz(-pi) q[2];
rz(1.9818241) q[3];
sx q[3];
rz(-0.81361249) q[3];
sx q[3];
rz(0.66853722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3942922) q[2];
sx q[2];
rz(-2.6960399) q[2];
sx q[2];
rz(-0.14075819) q[2];
rz(0.55177871) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(-0.75132918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20525876) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(-2.4667013) q[0];
rz(0.15448013) q[1];
sx q[1];
rz(-1.4090425) q[1];
sx q[1];
rz(2.6836269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58139153) q[0];
sx q[0];
rz(-2.0580252) q[0];
sx q[0];
rz(-2.013252) q[0];
x q[1];
rz(2.1764917) q[2];
sx q[2];
rz(-2.6332246) q[2];
sx q[2];
rz(-1.8338211) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.612466) q[1];
sx q[1];
rz(-1.7548029) q[1];
sx q[1];
rz(0.37671582) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6776039) q[3];
sx q[3];
rz(-1.065101) q[3];
sx q[3];
rz(3.0042574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13350479) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(-2.8221455) q[2];
rz(1.8114629) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(0.5602079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.9820246) q[0];
sx q[0];
rz(-1.1363131) q[0];
sx q[0];
rz(-1.2215479) q[0];
rz(2.9087861) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(-0.98308841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.074551) q[0];
sx q[0];
rz(-0.75387757) q[0];
sx q[0];
rz(0.75410944) q[0];
rz(-1.7295425) q[2];
sx q[2];
rz(-2.5149143) q[2];
sx q[2];
rz(1.9959297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5447588) q[1];
sx q[1];
rz(-1.4248409) q[1];
sx q[1];
rz(-0.18827855) q[1];
rz(-pi) q[2];
rz(-2.9748989) q[3];
sx q[3];
rz(-1.1536666) q[3];
sx q[3];
rz(-2.6670109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7588707) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(2.7663084) q[2];
rz(2.0659857) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(-0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7406834) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(1.4373454) q[0];
rz(-3.0628693) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(0.54814235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1931216) q[0];
sx q[0];
rz(-1.7821494) q[0];
sx q[0];
rz(-2.4708492) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70912045) q[2];
sx q[2];
rz(-0.68163423) q[2];
sx q[2];
rz(-0.97709419) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.50997616) q[1];
sx q[1];
rz(-2.2215034) q[1];
sx q[1];
rz(-1.0751616) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1918212) q[3];
sx q[3];
rz(-0.91663893) q[3];
sx q[3];
rz(-1.1266687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4569725) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(-0.6753298) q[2];
rz(1.5572549) q[3];
sx q[3];
rz(-1.3074338) q[3];
sx q[3];
rz(2.5244782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-2.6455998) q[0];
rz(-1.4542106) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(-2.0248263) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6750609) q[0];
sx q[0];
rz(-1.9218311) q[0];
sx q[0];
rz(2.768633) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1363813) q[2];
sx q[2];
rz(-2.5602753) q[2];
sx q[2];
rz(2.8535247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7796554) q[1];
sx q[1];
rz(-1.2314267) q[1];
sx q[1];
rz(-1.2324059) q[1];
rz(-0.12567606) q[3];
sx q[3];
rz(-0.45914859) q[3];
sx q[3];
rz(-2.7530376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3108959) q[2];
sx q[2];
rz(-0.98429698) q[2];
sx q[2];
rz(2.7867219) q[2];
rz(2.3762459) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(3.0520458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68882051) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(0.76559693) q[0];
rz(0.87144026) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(-1.3227468) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8091178) q[0];
sx q[0];
rz(-1.5189384) q[0];
sx q[0];
rz(-0.26082592) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9051612) q[2];
sx q[2];
rz(-0.90518236) q[2];
sx q[2];
rz(2.1588391) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31963667) q[1];
sx q[1];
rz(-1.0624891) q[1];
sx q[1];
rz(1.5896569) q[1];
rz(-1.345383) q[3];
sx q[3];
rz(-1.7724049) q[3];
sx q[3];
rz(-0.20568337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65779296) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(2.9198809) q[2];
rz(2.0486369) q[3];
sx q[3];
rz(-1.7287798) q[3];
sx q[3];
rz(-2.4634585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0677277) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(2.8644323) q[0];
rz(2.3932338) q[1];
sx q[1];
rz(-2.6942418) q[1];
sx q[1];
rz(-1.1433196) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049452) q[0];
sx q[0];
rz(-1.7713393) q[0];
sx q[0];
rz(-1.8822549) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3242993) q[2];
sx q[2];
rz(-0.72028226) q[2];
sx q[2];
rz(1.5718854) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9710247) q[1];
sx q[1];
rz(-1.2538365) q[1];
sx q[1];
rz(-0.16522878) q[1];
x q[2];
rz(-2.7128814) q[3];
sx q[3];
rz(-2.747251) q[3];
sx q[3];
rz(1.1129781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9144168) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(0.048967036) q[2];
rz(1.0660727) q[3];
sx q[3];
rz(-2.5132096) q[3];
sx q[3];
rz(-0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-0.12298909) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(-3.1220806) q[0];
rz(-0.77876577) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(-1.5628409) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922958) q[0];
sx q[0];
rz(-1.6598633) q[0];
sx q[0];
rz(-1.7013447) q[0];
rz(-pi) q[1];
rz(2.1807272) q[2];
sx q[2];
rz(-1.7210519) q[2];
sx q[2];
rz(-0.021266887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.409118) q[1];
sx q[1];
rz(-0.98608398) q[1];
sx q[1];
rz(0.52701108) q[1];
rz(-pi) q[2];
rz(0.83694038) q[3];
sx q[3];
rz(-1.6815974) q[3];
sx q[3];
rz(1.5821804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19266263) q[2];
sx q[2];
rz(-2.2781585) q[2];
sx q[2];
rz(-2.9760402) q[2];
rz(2.5340269) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(-2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3792569) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(-1.634585) q[1];
sx q[1];
rz(-1.4611117) q[1];
sx q[1];
rz(-1.5932105) q[1];
rz(-0.77539879) q[2];
sx q[2];
rz(-1.4682894) q[2];
sx q[2];
rz(-0.52649047) q[2];
rz(2.8589917) q[3];
sx q[3];
rz(-0.62487124) q[3];
sx q[3];
rz(1.665653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
