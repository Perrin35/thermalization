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
rz(3.1556407) q[1];
sx q[1];
rz(10.217759) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86092615) q[0];
sx q[0];
rz(-2.3287822) q[0];
sx q[0];
rz(2.0357512) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5231126) q[2];
sx q[2];
rz(-1.5850726) q[2];
sx q[2];
rz(2.2592777) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0539581) q[1];
sx q[1];
rz(-2.2184508) q[1];
sx q[1];
rz(-2.1235076) q[1];
x q[2];
rz(-2.1737384) q[3];
sx q[3];
rz(-2.0857028) q[3];
sx q[3];
rz(0.076857278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8636318) q[2];
sx q[2];
rz(-0.42676723) q[2];
sx q[2];
rz(2.5521736) q[2];
rz(2.3943118) q[3];
sx q[3];
rz(-0.092844754) q[3];
sx q[3];
rz(-2.5417627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079161949) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(-2.2404501) q[0];
rz(-0.98183739) q[1];
sx q[1];
rz(-2.5603309) q[1];
sx q[1];
rz(-0.99220401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81047786) q[0];
sx q[0];
rz(-1.7077291) q[0];
sx q[0];
rz(-2.1934319) q[0];
rz(0.82291863) q[2];
sx q[2];
rz(-1.4496921) q[2];
sx q[2];
rz(1.9117282) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7957669) q[1];
sx q[1];
rz(-2.2123033) q[1];
sx q[1];
rz(0.87916763) q[1];
rz(-0.76073356) q[3];
sx q[3];
rz(-2.2983716) q[3];
sx q[3];
rz(0.47970495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.09484) q[2];
sx q[2];
rz(-0.53043008) q[2];
sx q[2];
rz(-0.025010427) q[2];
rz(0.95995861) q[3];
sx q[3];
rz(-0.046006087) q[3];
sx q[3];
rz(-0.068232603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18441021) q[0];
sx q[0];
rz(-0.010951696) q[0];
sx q[0];
rz(-0.72407323) q[0];
rz(-3.030576) q[1];
sx q[1];
rz(-0.60581517) q[1];
sx q[1];
rz(-0.012880005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89319481) q[0];
sx q[0];
rz(-1.7832558) q[0];
sx q[0];
rz(-0.080470632) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.951163) q[2];
sx q[2];
rz(-1.3821175) q[2];
sx q[2];
rz(-1.7469454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6659505) q[1];
sx q[1];
rz(-2.8503909) q[1];
sx q[1];
rz(2.5282574) q[1];
rz(-pi) q[2];
rz(-2.0184085) q[3];
sx q[3];
rz(-0.73799282) q[3];
sx q[3];
rz(0.18692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24590242) q[2];
sx q[2];
rz(-2.4195713) q[2];
sx q[2];
rz(0.32430172) q[2];
rz(0.4970099) q[3];
sx q[3];
rz(-2.9091166) q[3];
sx q[3];
rz(-2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683891) q[0];
sx q[0];
rz(-2.9717746) q[0];
sx q[0];
rz(-0.47624269) q[0];
rz(-0.4854804) q[1];
sx q[1];
rz(-0.49518934) q[1];
sx q[1];
rz(0.21864299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63594288) q[0];
sx q[0];
rz(-0.46920645) q[0];
sx q[0];
rz(2.0903265) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4179732) q[2];
sx q[2];
rz(-1.5800522) q[2];
sx q[2];
rz(-1.1030359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8805595) q[1];
sx q[1];
rz(-1.1487242) q[1];
sx q[1];
rz(2.5690394) q[1];
x q[2];
rz(-1.8044295) q[3];
sx q[3];
rz(-1.8770185) q[3];
sx q[3];
rz(-2.985317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.789088) q[2];
sx q[2];
rz(-2.3829491) q[2];
sx q[2];
rz(0.57141203) q[2];
rz(2.2032951) q[3];
sx q[3];
rz(-0.58656991) q[3];
sx q[3];
rz(0.17294426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5649696) q[0];
sx q[0];
rz(-2.7374856) q[0];
sx q[0];
rz(0.47082666) q[0];
rz(2.2305523) q[1];
sx q[1];
rz(-0.43993479) q[1];
sx q[1];
rz(-2.7104673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.39931) q[0];
sx q[0];
rz(-2.3547908) q[0];
sx q[0];
rz(1.6086701) q[0];
rz(2.8887553) q[2];
sx q[2];
rz(-0.95273861) q[2];
sx q[2];
rz(-1.1065266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87157618) q[1];
sx q[1];
rz(-1.5983014) q[1];
sx q[1];
rz(-1.639495) q[1];
rz(0.2016368) q[3];
sx q[3];
rz(-1.7072225) q[3];
sx q[3];
rz(-1.0659983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8016781) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(2.8002296) q[2];
rz(-0.3723799) q[3];
sx q[3];
rz(-2.5058993) q[3];
sx q[3];
rz(-0.46975964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.85207283) q[0];
sx q[0];
rz(-2.2884123) q[0];
sx q[0];
rz(-0.12292718) q[0];
rz(0.3854824) q[1];
sx q[1];
rz(-0.79817927) q[1];
sx q[1];
rz(-0.72365671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22985499) q[0];
sx q[0];
rz(-2.9617972) q[0];
sx q[0];
rz(-1.3151965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18823024) q[2];
sx q[2];
rz(-1.2410795) q[2];
sx q[2];
rz(-2.4356869) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21863114) q[1];
sx q[1];
rz(-0.9762763) q[1];
sx q[1];
rz(0.84719744) q[1];
rz(-2.0615929) q[3];
sx q[3];
rz(-0.43269581) q[3];
sx q[3];
rz(-3.1263292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75189292) q[2];
sx q[2];
rz(-2.5395826) q[2];
sx q[2];
rz(0.67328084) q[2];
rz(-2.8781387) q[3];
sx q[3];
rz(-2.6665688) q[3];
sx q[3];
rz(-2.8365005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4398956) q[0];
sx q[0];
rz(-2.3450527) q[0];
sx q[0];
rz(-0.51629603) q[0];
rz(0.75597489) q[1];
sx q[1];
rz(-2.1831903) q[1];
sx q[1];
rz(-2.7730673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4004423) q[0];
sx q[0];
rz(-2.0870805) q[0];
sx q[0];
rz(-2.3325066) q[0];
rz(-pi) q[1];
rz(-0.045727878) q[2];
sx q[2];
rz(-0.70529443) q[2];
sx q[2];
rz(-0.7208342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3893551) q[1];
sx q[1];
rz(-0.19986831) q[1];
sx q[1];
rz(-1.0081069) q[1];
rz(-pi) q[2];
rz(-0.51492274) q[3];
sx q[3];
rz(-0.92149261) q[3];
sx q[3];
rz(2.2442371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4964909) q[2];
sx q[2];
rz(-0.056702159) q[2];
sx q[2];
rz(-0.48218316) q[2];
rz(-0.21969806) q[3];
sx q[3];
rz(-2.4441661) q[3];
sx q[3];
rz(-0.68035948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7242303) q[0];
sx q[0];
rz(-0.14886947) q[0];
sx q[0];
rz(2.505488) q[0];
rz(-0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(-2.2669534) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2518305) q[0];
sx q[0];
rz(-1.5127851) q[0];
sx q[0];
rz(-1.4213575) q[0];
rz(-pi) q[1];
rz(0.12107559) q[2];
sx q[2];
rz(-1.4619451) q[2];
sx q[2];
rz(-1.7055275) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1369774) q[1];
sx q[1];
rz(-0.57356131) q[1];
sx q[1];
rz(-2.3358845) q[1];
rz(-pi) q[2];
rz(2.2611924) q[3];
sx q[3];
rz(-2.1244308) q[3];
sx q[3];
rz(2.2601489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26802289) q[2];
sx q[2];
rz(-0.41646725) q[2];
sx q[2];
rz(2.2250309) q[2];
rz(0.77740866) q[3];
sx q[3];
rz(-2.58367) q[3];
sx q[3];
rz(0.53434813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022631835) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(0.23271261) q[0];
rz(-2.4932056) q[1];
sx q[1];
rz(-0.62260038) q[1];
sx q[1];
rz(-0.56786215) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647301) q[0];
sx q[0];
rz(-1.1608539) q[0];
sx q[0];
rz(-0.38500824) q[0];
rz(-pi) q[1];
rz(-1.2987192) q[2];
sx q[2];
rz(-1.7781742) q[2];
sx q[2];
rz(-3.1081215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4109013) q[1];
sx q[1];
rz(-1.7736378) q[1];
sx q[1];
rz(-2.9922725) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80983272) q[3];
sx q[3];
rz(-2.2264635) q[3];
sx q[3];
rz(1.9698433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3719172) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(-0.5522716) q[2];
rz(-2.8596089) q[3];
sx q[3];
rz(-0.43006399) q[3];
sx q[3];
rz(3.0388487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2118537) q[0];
sx q[0];
rz(-0.064082853) q[0];
sx q[0];
rz(0.1317568) q[0];
rz(-0.53648221) q[1];
sx q[1];
rz(-0.15833144) q[1];
sx q[1];
rz(-2.2333455) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2836766) q[0];
sx q[0];
rz(-2.2649797) q[0];
sx q[0];
rz(3.1026918) q[0];
rz(-pi) q[1];
rz(-0.26226173) q[2];
sx q[2];
rz(-1.1168257) q[2];
sx q[2];
rz(-1.5020467) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5436369) q[1];
sx q[1];
rz(-2.0858793) q[1];
sx q[1];
rz(-2.3233344) q[1];
rz(-pi) q[2];
rz(-3.1397758) q[3];
sx q[3];
rz(-1.7967136) q[3];
sx q[3];
rz(-0.29669138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50916719) q[2];
sx q[2];
rz(-0.91370344) q[2];
sx q[2];
rz(-0.2779648) q[2];
rz(-0.5667423) q[3];
sx q[3];
rz(-0.20213474) q[3];
sx q[3];
rz(-2.22866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3315898) q[0];
sx q[0];
rz(-1.7362052) q[0];
sx q[0];
rz(2.195634) q[0];
rz(-2.6592061) q[1];
sx q[1];
rz(-1.2336029) q[1];
sx q[1];
rz(-0.69419669) q[1];
rz(-2.161088) q[2];
sx q[2];
rz(-1.7750778) q[2];
sx q[2];
rz(1.9951174) q[2];
rz(1.514074) q[3];
sx q[3];
rz(-2.6346907) q[3];
sx q[3];
rz(-2.9148867) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
