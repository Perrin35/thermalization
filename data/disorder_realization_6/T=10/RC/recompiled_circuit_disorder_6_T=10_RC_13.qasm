OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(4.5594112) q[0];
sx q[0];
rz(8.863908) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(-1.2150432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96072223) q[0];
sx q[0];
rz(-1.7579494) q[0];
sx q[0];
rz(-0.10589177) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4807793) q[2];
sx q[2];
rz(-1.3408957) q[2];
sx q[2];
rz(1.7125318) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6601662) q[1];
sx q[1];
rz(-0.62454849) q[1];
sx q[1];
rz(0.98510965) q[1];
x q[2];
rz(0.14903544) q[3];
sx q[3];
rz(-1.8334853) q[3];
sx q[3];
rz(-1.5460154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(-2.7637774) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(3.0644754) q[0];
rz(0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.6024626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9248283) q[0];
sx q[0];
rz(-0.98470682) q[0];
sx q[0];
rz(-3.0490962) q[0];
rz(-pi) q[1];
rz(1.3537172) q[2];
sx q[2];
rz(-2.2976544) q[2];
sx q[2];
rz(-0.94490563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96670818) q[1];
sx q[1];
rz(-1.0918573) q[1];
sx q[1];
rz(0.19660463) q[1];
rz(-pi) q[2];
rz(-0.42373557) q[3];
sx q[3];
rz(-1.7859922) q[3];
sx q[3];
rz(-2.7620897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-2.4831333) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-2.4466799) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(-2.5862397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99291891) q[0];
sx q[0];
rz(-2.2745471) q[0];
sx q[0];
rz(2.2253195) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0981512) q[2];
sx q[2];
rz(-2.522509) q[2];
sx q[2];
rz(1.667779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66545031) q[1];
sx q[1];
rz(-0.94201554) q[1];
sx q[1];
rz(2.9412342) q[1];
x q[2];
rz(1.8954574) q[3];
sx q[3];
rz(-2.5723296) q[3];
sx q[3];
rz(-3.0666921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.4499433) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(-1.3793777) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(0.25517685) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53084757) q[0];
sx q[0];
rz(-1.2753715) q[0];
sx q[0];
rz(0.3770963) q[0];
rz(-1.2404664) q[2];
sx q[2];
rz(-0.9579881) q[2];
sx q[2];
rz(-1.4532879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55528044) q[1];
sx q[1];
rz(-1.3603856) q[1];
sx q[1];
rz(-1.8754688) q[1];
rz(-pi) q[2];
x q[2];
rz(0.051024036) q[3];
sx q[3];
rz(-2.0912366) q[3];
sx q[3];
rz(0.40363064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(0.17318428) q[2];
rz(-2.611768) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(-1.3758855) q[0];
rz(-1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(-0.056093562) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.010904) q[0];
sx q[0];
rz(-2.4405257) q[0];
sx q[0];
rz(-2.5581193) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34727879) q[2];
sx q[2];
rz(-1.1355073) q[2];
sx q[2];
rz(1.6821282) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5406815) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(1.5555698) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0523473) q[3];
sx q[3];
rz(-1.0122932) q[3];
sx q[3];
rz(-0.23111471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0284001) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(0.64754852) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19676767) q[0];
sx q[0];
rz(-1.7260572) q[0];
sx q[0];
rz(-1.0374271) q[0];
rz(-3.1068222) q[2];
sx q[2];
rz(-0.94049847) q[2];
sx q[2];
rz(0.45331732) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2973605) q[1];
sx q[1];
rz(-0.73207049) q[1];
sx q[1];
rz(0.75433235) q[1];
rz(2.9783863) q[3];
sx q[3];
rz(-1.8582134) q[3];
sx q[3];
rz(0.67374574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(1.6339533) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9305206) q[0];
sx q[0];
rz(-2.0606344) q[0];
sx q[0];
rz(1.8563489) q[0];
rz(-3.0257752) q[2];
sx q[2];
rz(-0.91997416) q[2];
sx q[2];
rz(-1.7388369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66079084) q[1];
sx q[1];
rz(-1.5457488) q[1];
sx q[1];
rz(0.90344306) q[1];
x q[2];
rz(-2.0537297) q[3];
sx q[3];
rz(-0.42290877) q[3];
sx q[3];
rz(-1.203323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(-0.91840333) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36088762) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.5135182) q[0];
rz(2.6121415) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(2.4050074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442473) q[0];
sx q[0];
rz(-1.559942) q[0];
sx q[0];
rz(1.6008475) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.893115) q[2];
sx q[2];
rz(-1.6741447) q[2];
sx q[2];
rz(1.5600187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.089162) q[1];
sx q[1];
rz(-0.4757291) q[1];
sx q[1];
rz(0.22389852) q[1];
rz(-0.96111091) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(-1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(2.4576808) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(0.99705237) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(2.396778) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1574402) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(1.1293344) q[0];
rz(-pi) q[1];
rz(1.1698193) q[2];
sx q[2];
rz(-0.81911659) q[2];
sx q[2];
rz(0.70105201) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77764952) q[1];
sx q[1];
rz(-1.793982) q[1];
sx q[1];
rz(1.2121483) q[1];
rz(-pi) q[2];
rz(-2.3585988) q[3];
sx q[3];
rz(-1.8564312) q[3];
sx q[3];
rz(2.1180958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.194681) q[2];
rz(-2.1448994) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(0.031127302) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(-1.9706479) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084328018) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(-1.5945934) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3002214) q[2];
sx q[2];
rz(-0.8136533) q[2];
sx q[2];
rz(-0.68683456) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78466258) q[1];
sx q[1];
rz(-1.5481879) q[1];
sx q[1];
rz(-0.008634062) q[1];
rz(-pi) q[2];
rz(2.714614) q[3];
sx q[3];
rz(-1.3037762) q[3];
sx q[3];
rz(-3.0468575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(-0.5919624) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82407172) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(0.099427632) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(0.912491) q[2];
sx q[2];
rz(-0.9463263) q[2];
sx q[2];
rz(1.8778388) q[2];
rz(3.1254461) q[3];
sx q[3];
rz(-1.2373677) q[3];
sx q[3];
rz(2.2131372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
