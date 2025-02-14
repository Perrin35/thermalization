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
rz(0.1306611) q[0];
sx q[0];
rz(4.9424439) q[0];
sx q[0];
rz(10.688936) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(-1.0882508) q[1];
sx q[1];
rz(-3.1286214) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1931827) q[0];
sx q[0];
rz(-1.7275066) q[0];
sx q[0];
rz(2.66314) q[0];
rz(-pi) q[1];
rz(2.5532236) q[2];
sx q[2];
rz(-0.19489637) q[2];
sx q[2];
rz(0.68258679) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88718677) q[1];
sx q[1];
rz(-1.9209083) q[1];
sx q[1];
rz(1.9113879) q[1];
x q[2];
rz(0.31855984) q[3];
sx q[3];
rz(-2.5346642) q[3];
sx q[3];
rz(-1.3598809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.26244792) q[2];
sx q[2];
rz(-1.4538572) q[2];
sx q[2];
rz(3.0461779) q[2];
rz(2.6573507) q[3];
sx q[3];
rz(-0.80317322) q[3];
sx q[3];
rz(1.6835083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5937623) q[0];
sx q[0];
rz(-1.7050803) q[0];
sx q[0];
rz(2.4875212) q[0];
rz(-1.7895128) q[1];
sx q[1];
rz(-2.4857931) q[1];
sx q[1];
rz(0.10993122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1760565) q[0];
sx q[0];
rz(-1.5586462) q[0];
sx q[0];
rz(-1.5506844) q[0];
rz(-pi) q[1];
x q[1];
rz(2.201033) q[2];
sx q[2];
rz(-1.1917243) q[2];
sx q[2];
rz(0.68507588) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9603574) q[1];
sx q[1];
rz(-1.7434038) q[1];
sx q[1];
rz(2.6828275) q[1];
rz(-pi) q[2];
rz(0.55628784) q[3];
sx q[3];
rz(-2.6835052) q[3];
sx q[3];
rz(2.213495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.316651) q[2];
sx q[2];
rz(-1.2085088) q[2];
sx q[2];
rz(1.6744772) q[2];
rz(1.0351099) q[3];
sx q[3];
rz(-2.3605774) q[3];
sx q[3];
rz(1.3181814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255945) q[0];
sx q[0];
rz(-2.9075629) q[0];
sx q[0];
rz(0.79063928) q[0];
rz(0.74360338) q[1];
sx q[1];
rz(-1.5433106) q[1];
sx q[1];
rz(-2.5620983) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7121237) q[0];
sx q[0];
rz(-2.3050024) q[0];
sx q[0];
rz(0.81487687) q[0];
rz(-2.2966301) q[2];
sx q[2];
rz(-2.7562592) q[2];
sx q[2];
rz(-2.8266738) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0708376) q[1];
sx q[1];
rz(-2.8808314) q[1];
sx q[1];
rz(-1.7160796) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0813339) q[3];
sx q[3];
rz(-2.0474985) q[3];
sx q[3];
rz(-2.7923194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6774595) q[2];
sx q[2];
rz(-0.83659283) q[2];
sx q[2];
rz(-0.38132384) q[2];
rz(0.85121202) q[3];
sx q[3];
rz(-2.2150453) q[3];
sx q[3];
rz(0.73494953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5320936) q[0];
sx q[0];
rz(-2.5272326) q[0];
sx q[0];
rz(-0.36439782) q[0];
rz(2.9413307) q[1];
sx q[1];
rz(-1.4118782) q[1];
sx q[1];
rz(-3.0416378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1276949) q[0];
sx q[0];
rz(-1.5911926) q[0];
sx q[0];
rz(1.6457895) q[0];
rz(-1.3447138) q[2];
sx q[2];
rz(-0.43856171) q[2];
sx q[2];
rz(1.3815224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8980838) q[1];
sx q[1];
rz(-1.5484515) q[1];
sx q[1];
rz(0.0925272) q[1];
rz(-0.95619802) q[3];
sx q[3];
rz(-1.0805939) q[3];
sx q[3];
rz(-0.7008926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.072307) q[2];
sx q[2];
rz(-1.0682718) q[2];
sx q[2];
rz(0.11492534) q[2];
rz(-2.0011486) q[3];
sx q[3];
rz(-1.2830257) q[3];
sx q[3];
rz(-2.1663402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045778) q[0];
sx q[0];
rz(-0.95252043) q[0];
sx q[0];
rz(0.17247795) q[0];
rz(1.0109488) q[1];
sx q[1];
rz(-0.5609678) q[1];
sx q[1];
rz(0.15484658) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2456937) q[0];
sx q[0];
rz(-1.8666728) q[0];
sx q[0];
rz(-0.22801836) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92645433) q[2];
sx q[2];
rz(-1.3727034) q[2];
sx q[2];
rz(-1.4938172) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6754969) q[1];
sx q[1];
rz(-1.4267529) q[1];
sx q[1];
rz(0.41464582) q[1];
rz(-pi) q[2];
rz(1.6825292) q[3];
sx q[3];
rz(-1.3936685) q[3];
sx q[3];
rz(-0.88645173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.054472063) q[2];
sx q[2];
rz(-0.27721578) q[2];
sx q[2];
rz(-2.7692914) q[2];
rz(2.7436658) q[3];
sx q[3];
rz(-1.9979265) q[3];
sx q[3];
rz(-3.1411689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.953124) q[0];
sx q[0];
rz(-2.5690881) q[0];
sx q[0];
rz(1.5484126) q[0];
rz(0.36987034) q[1];
sx q[1];
rz(-1.3984171) q[1];
sx q[1];
rz(-0.63873783) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8014378) q[0];
sx q[0];
rz(-1.7246913) q[0];
sx q[0];
rz(1.811054) q[0];
x q[1];
rz(2.0962786) q[2];
sx q[2];
rz(-1.3945701) q[2];
sx q[2];
rz(0.8187364) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46203278) q[1];
sx q[1];
rz(-2.5278494) q[1];
sx q[1];
rz(-0.8590974) q[1];
rz(-2.3494598) q[3];
sx q[3];
rz(-1.9895305) q[3];
sx q[3];
rz(-1.2611539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0729735) q[2];
sx q[2];
rz(-1.0516473) q[2];
sx q[2];
rz(0.2529141) q[2];
rz(-2.2933293) q[3];
sx q[3];
rz(-1.6631283) q[3];
sx q[3];
rz(-0.084913582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0406168) q[0];
sx q[0];
rz(-0.210013) q[0];
sx q[0];
rz(2.9926391) q[0];
rz(1.8303998) q[1];
sx q[1];
rz(-1.4102178) q[1];
sx q[1];
rz(0.054100903) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.427092) q[0];
sx q[0];
rz(-0.50760554) q[0];
sx q[0];
rz(-0.18113329) q[0];
x q[1];
rz(-0.04345036) q[2];
sx q[2];
rz(-0.80294007) q[2];
sx q[2];
rz(2.174365) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91840345) q[1];
sx q[1];
rz(-1.4610054) q[1];
sx q[1];
rz(0.70007433) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6160433) q[3];
sx q[3];
rz(-2.6658969) q[3];
sx q[3];
rz(3.0632116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80186239) q[2];
sx q[2];
rz(-1.9714377) q[2];
sx q[2];
rz(-2.1745963) q[2];
rz(-2.4407834) q[3];
sx q[3];
rz(-2.1613224) q[3];
sx q[3];
rz(2.7907659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0130149) q[0];
sx q[0];
rz(-1.1739434) q[0];
sx q[0];
rz(1.5933734) q[0];
rz(2.7633527) q[1];
sx q[1];
rz(-0.75235569) q[1];
sx q[1];
rz(-0.38984782) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10791099) q[0];
sx q[0];
rz(-2.0706688) q[0];
sx q[0];
rz(1.3477911) q[0];
rz(-pi) q[1];
x q[1];
rz(0.016640113) q[2];
sx q[2];
rz(-1.2667873) q[2];
sx q[2];
rz(-1.9250411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2521542) q[1];
sx q[1];
rz(-0.61823003) q[1];
sx q[1];
rz(-0.63207027) q[1];
rz(2.8865783) q[3];
sx q[3];
rz(-1.3552291) q[3];
sx q[3];
rz(2.4795585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0491911) q[2];
sx q[2];
rz(-2.0202049) q[2];
sx q[2];
rz(1.258705) q[2];
rz(0.24426584) q[3];
sx q[3];
rz(-1.786307) q[3];
sx q[3];
rz(0.99610966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8213537) q[0];
sx q[0];
rz(-1.2293674) q[0];
sx q[0];
rz(1.356333) q[0];
rz(-2.9755196) q[1];
sx q[1];
rz(-0.95172721) q[1];
sx q[1];
rz(-0.43620268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22101519) q[0];
sx q[0];
rz(-1.2599143) q[0];
sx q[0];
rz(0.47500821) q[0];
rz(-1.9848083) q[2];
sx q[2];
rz(-1.8824667) q[2];
sx q[2];
rz(-3.0729635) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2959152) q[1];
sx q[1];
rz(-1.3406173) q[1];
sx q[1];
rz(-0.796411) q[1];
rz(-pi) q[2];
rz(-0.43667359) q[3];
sx q[3];
rz(-1.5092351) q[3];
sx q[3];
rz(1.7583282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0210375) q[2];
sx q[2];
rz(-1.1641116) q[2];
sx q[2];
rz(-2.5980921) q[2];
rz(-0.049526878) q[3];
sx q[3];
rz(-1.6560358) q[3];
sx q[3];
rz(1.4878976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5307584) q[0];
sx q[0];
rz(-0.13228358) q[0];
sx q[0];
rz(-0.079205967) q[0];
rz(-0.40009701) q[1];
sx q[1];
rz(-0.85405093) q[1];
sx q[1];
rz(-1.0265464) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.067018) q[0];
sx q[0];
rz(-1.3368784) q[0];
sx q[0];
rz(2.3090655) q[0];
rz(-1.6658989) q[2];
sx q[2];
rz(-1.9355023) q[2];
sx q[2];
rz(2.022418) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3803097) q[1];
sx q[1];
rz(-1.6649462) q[1];
sx q[1];
rz(0.32431079) q[1];
rz(-pi) q[2];
rz(-0.35086326) q[3];
sx q[3];
rz(-1.4961582) q[3];
sx q[3];
rz(-0.78001744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93572271) q[2];
sx q[2];
rz(-1.8368072) q[2];
sx q[2];
rz(1.4523466) q[2];
rz(2.3116889) q[3];
sx q[3];
rz(-0.87703505) q[3];
sx q[3];
rz(-0.95480603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6185388) q[0];
sx q[0];
rz(-1.9970311) q[0];
sx q[0];
rz(1.507623) q[0];
rz(1.2670831) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(-0.53955033) q[2];
sx q[2];
rz(-2.087941) q[2];
sx q[2];
rz(2.6826442) q[2];
rz(1.2290365) q[3];
sx q[3];
rz(-1.0377586) q[3];
sx q[3];
rz(0.23317045) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
