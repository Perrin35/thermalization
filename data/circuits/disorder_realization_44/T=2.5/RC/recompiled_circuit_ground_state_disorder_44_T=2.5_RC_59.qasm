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
rz(-2.1704817) q[0];
sx q[0];
rz(0.71075332) q[0];
rz(-2.1739668) q[1];
sx q[1];
rz(-0.014048014) q[1];
sx q[1];
rz(-0.79298055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86092615) q[0];
sx q[0];
rz(-2.3287822) q[0];
sx q[0];
rz(1.1058415) q[0];
rz(-0.014292467) q[2];
sx q[2];
rz(-1.6184752) q[2];
sx q[2];
rz(-2.45243) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2574764) q[1];
sx q[1];
rz(-2.3168644) q[1];
sx q[1];
rz(-2.5349239) q[1];
rz(2.3550911) q[3];
sx q[3];
rz(-2.3701027) q[3];
sx q[3];
rz(0.87302333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8636318) q[2];
sx q[2];
rz(-2.7148254) q[2];
sx q[2];
rz(0.58941907) q[2];
rz(-2.3943118) q[3];
sx q[3];
rz(-0.092844754) q[3];
sx q[3];
rz(-0.59982991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079161949) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(-2.2404501) q[0];
rz(2.1597553) q[1];
sx q[1];
rz(-2.5603309) q[1];
sx q[1];
rz(-0.99220401) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3311148) q[0];
sx q[0];
rz(-1.4338636) q[0];
sx q[0];
rz(-0.94816072) q[0];
x q[1];
rz(2.9770933) q[2];
sx q[2];
rz(-2.3119011) q[2];
sx q[2];
rz(-2.6890597) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34582576) q[1];
sx q[1];
rz(-2.2123033) q[1];
sx q[1];
rz(-0.87916763) q[1];
rz(-pi) q[2];
rz(-2.2295932) q[3];
sx q[3];
rz(-0.9992632) q[3];
sx q[3];
rz(2.6603572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.046752669) q[2];
sx q[2];
rz(-2.6111626) q[2];
sx q[2];
rz(0.025010427) q[2];
rz(0.95995861) q[3];
sx q[3];
rz(-3.0955866) q[3];
sx q[3];
rz(-3.0733601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18441021) q[0];
sx q[0];
rz(-0.010951696) q[0];
sx q[0];
rz(2.4175194) q[0];
rz(0.11101668) q[1];
sx q[1];
rz(-2.5357775) q[1];
sx q[1];
rz(-3.1287126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89319481) q[0];
sx q[0];
rz(-1.7832558) q[0];
sx q[0];
rz(0.080470632) q[0];
rz(-pi) q[1];
rz(0.19042966) q[2];
sx q[2];
rz(-1.7594751) q[2];
sx q[2];
rz(-1.3946472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6659505) q[1];
sx q[1];
rz(-0.29120177) q[1];
sx q[1];
rz(-0.61333527) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1231842) q[3];
sx q[3];
rz(-2.4035998) q[3];
sx q[3];
rz(-2.9546705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8956902) q[2];
sx q[2];
rz(-2.4195713) q[2];
sx q[2];
rz(-2.8172909) q[2];
rz(-2.6445828) q[3];
sx q[3];
rz(-2.9091166) q[3];
sx q[3];
rz(-2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4683891) q[0];
sx q[0];
rz(-2.9717746) q[0];
sx q[0];
rz(-2.66535) q[0];
rz(-2.6561123) q[1];
sx q[1];
rz(-0.49518934) q[1];
sx q[1];
rz(-0.21864299) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9354585) q[0];
sx q[0];
rz(-1.1674321) q[0];
sx q[0];
rz(2.8950188) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72361941) q[2];
sx q[2];
rz(-1.5800522) q[2];
sx q[2];
rz(2.0385567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87590295) q[1];
sx q[1];
rz(-2.4445718) q[1];
sx q[1];
rz(2.4494996) q[1];
rz(-1.3371631) q[3];
sx q[3];
rz(-1.2645742) q[3];
sx q[3];
rz(-2.985317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.789088) q[2];
sx q[2];
rz(-0.75864351) q[2];
sx q[2];
rz(-2.5701806) q[2];
rz(2.2032951) q[3];
sx q[3];
rz(-2.5550227) q[3];
sx q[3];
rz(2.9686484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(2.5649696) q[0];
sx q[0];
rz(-2.7374856) q[0];
sx q[0];
rz(-0.47082666) q[0];
rz(0.91104031) q[1];
sx q[1];
rz(-0.43993479) q[1];
sx q[1];
rz(2.7104673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8017641) q[0];
sx q[0];
rz(-1.5439811) q[0];
sx q[0];
rz(-2.3572395) q[0];
rz(0.25283739) q[2];
sx q[2];
rz(-0.95273861) q[2];
sx q[2];
rz(-2.035066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.70111245) q[1];
sx q[1];
rz(-1.639469) q[1];
sx q[1];
rz(0.027570034) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9399559) q[3];
sx q[3];
rz(-1.7072225) q[3];
sx q[3];
rz(1.0659983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8016781) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(-2.8002296) q[2];
rz(0.3723799) q[3];
sx q[3];
rz(-0.63569331) q[3];
sx q[3];
rz(-0.46975964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85207283) q[0];
sx q[0];
rz(-2.2884123) q[0];
sx q[0];
rz(-3.0186655) q[0];
rz(0.3854824) q[1];
sx q[1];
rz(-2.3434134) q[1];
sx q[1];
rz(0.72365671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.052304) q[0];
sx q[0];
rz(-1.6160242) q[0];
sx q[0];
rz(-1.3967229) q[0];
rz(1.2355818) q[2];
sx q[2];
rz(-1.3928169) q[2];
sx q[2];
rz(2.2151057) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3544729) q[1];
sx q[1];
rz(-2.2406089) q[1];
sx q[1];
rz(-0.77490999) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0615929) q[3];
sx q[3];
rz(-2.7088968) q[3];
sx q[3];
rz(3.1263292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3896997) q[2];
sx q[2];
rz(-0.60201001) q[2];
sx q[2];
rz(2.4683118) q[2];
rz(2.8781387) q[3];
sx q[3];
rz(-0.47502381) q[3];
sx q[3];
rz(-2.8365005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70169705) q[0];
sx q[0];
rz(-0.79653996) q[0];
sx q[0];
rz(-2.6252966) q[0];
rz(-2.3856178) q[1];
sx q[1];
rz(-2.1831903) q[1];
sx q[1];
rz(0.36852536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26966345) q[0];
sx q[0];
rz(-0.9270398) q[0];
sx q[0];
rz(-2.4764349) q[0];
x q[1];
rz(1.6096949) q[2];
sx q[2];
rz(-0.86639154) q[2];
sx q[2];
rz(2.4807841) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26490274) q[1];
sx q[1];
rz(-1.4646834) q[1];
sx q[1];
rz(-1.7404895) q[1];
rz(-2.1462899) q[3];
sx q[3];
rz(-0.80484521) q[3];
sx q[3];
rz(0.14508776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64510173) q[2];
sx q[2];
rz(-0.056702159) q[2];
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
x q[1];
sx q[2];
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
rz(-2.2669534) q[1];
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
x q[1];
rz(-3.0205171) q[2];
sx q[2];
rz(-1.6796475) q[2];
sx q[2];
rz(-1.4360652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.85669163) q[1];
sx q[1];
rz(-1.1686348) q[1];
sx q[1];
rz(0.42070893) q[1];
x q[2];
rz(2.465807) q[3];
sx q[3];
rz(-0.99832557) q[3];
sx q[3];
rz(-1.0990717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26802289) q[2];
sx q[2];
rz(-2.7251254) q[2];
sx q[2];
rz(2.2250309) q[2];
rz(2.364184) q[3];
sx q[3];
rz(-0.55792266) q[3];
sx q[3];
rz(0.53434813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1189608) q[0];
sx q[0];
rz(-2.4038778) q[0];
sx q[0];
rz(-2.90888) q[0];
rz(0.64838707) q[1];
sx q[1];
rz(-0.62260038) q[1];
sx q[1];
rz(-0.56786215) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083160087) q[0];
sx q[0];
rz(-0.55480236) q[0];
sx q[0];
rz(-0.85809274) q[0];
x q[1];
rz(-0.21505996) q[2];
sx q[2];
rz(-1.304692) q[2];
sx q[2];
rz(-1.6616481) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.372124) q[1];
sx q[1];
rz(-0.25126496) q[1];
sx q[1];
rz(-0.94458802) q[1];
rz(-pi) q[2];
rz(2.4106823) q[3];
sx q[3];
rz(-0.95940351) q[3];
sx q[3];
rz(-2.9711593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76967543) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(-2.589321) q[2];
rz(2.8596089) q[3];
sx q[3];
rz(-0.43006399) q[3];
sx q[3];
rz(0.10274398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.2118537) q[0];
sx q[0];
rz(-3.0775098) q[0];
sx q[0];
rz(3.0098359) q[0];
rz(0.53648221) q[1];
sx q[1];
rz(-2.9832612) q[1];
sx q[1];
rz(0.90824711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73777481) q[0];
sx q[0];
rz(-1.5409011) q[0];
sx q[0];
rz(-0.87624082) q[0];
x q[1];
rz(-2.0591878) q[2];
sx q[2];
rz(-2.6219212) q[2];
sx q[2];
rz(2.1888417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7365807) q[1];
sx q[1];
rz(-0.93376505) q[1];
sx q[1];
rz(0.65959658) q[1];
x q[2];
rz(-1.3448787) q[3];
sx q[3];
rz(-1.572567) q[3];
sx q[3];
rz(1.8670807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6324255) q[2];
sx q[2];
rz(-2.2278892) q[2];
sx q[2];
rz(-2.8636279) q[2];
rz(0.5667423) q[3];
sx q[3];
rz(-2.9394579) q[3];
sx q[3];
rz(0.91293269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81000281) q[0];
sx q[0];
rz(-1.4053874) q[0];
sx q[0];
rz(-0.94595861) q[0];
rz(2.6592061) q[1];
sx q[1];
rz(-1.9079897) q[1];
sx q[1];
rz(2.447396) q[1];
rz(-2.8972068) q[2];
sx q[2];
rz(-2.1472211) q[2];
sx q[2];
rz(0.28923464) q[2];
rz(2.0770155) q[3];
sx q[3];
rz(-1.5432705) q[3];
sx q[3];
rz(-1.3936925) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
