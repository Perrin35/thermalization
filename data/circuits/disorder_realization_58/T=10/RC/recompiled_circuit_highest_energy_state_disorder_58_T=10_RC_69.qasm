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
rz(-0.11197055) q[0];
sx q[0];
rz(-2.114871) q[0];
sx q[0];
rz(1.8486899) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(-0.81967241) q[1];
sx q[1];
rz(1.1991062) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0831868) q[0];
sx q[0];
rz(-1.633344) q[0];
sx q[0];
rz(0.10857554) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67352301) q[2];
sx q[2];
rz(-1.4269514) q[2];
sx q[2];
rz(2.0604482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0173229) q[1];
sx q[1];
rz(-2.5376179) q[1];
sx q[1];
rz(-1.225597) q[1];
rz(-pi) q[2];
rz(0.85924826) q[3];
sx q[3];
rz(-1.0785127) q[3];
sx q[3];
rz(1.8627721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5045515) q[2];
sx q[2];
rz(-1.7402288) q[2];
sx q[2];
rz(-1.6999647) q[2];
rz(-2.1291034) q[3];
sx q[3];
rz(-1.9952521) q[3];
sx q[3];
rz(2.6700524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3016894) q[0];
sx q[0];
rz(-2.58044) q[0];
sx q[0];
rz(2.1121693) q[0];
rz(-1.3599716) q[1];
sx q[1];
rz(-1.9554892) q[1];
sx q[1];
rz(-0.50672466) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6678807) q[0];
sx q[0];
rz(-0.10598826) q[0];
sx q[0];
rz(1.7539658) q[0];
rz(-0.87575715) q[2];
sx q[2];
rz(-1.8570261) q[2];
sx q[2];
rz(1.8935668) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4480865) q[1];
sx q[1];
rz(-0.58219456) q[1];
sx q[1];
rz(-1.4513998) q[1];
rz(0.84552879) q[3];
sx q[3];
rz(-1.3432028) q[3];
sx q[3];
rz(-1.1448154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8257251) q[2];
sx q[2];
rz(-2.5434912) q[2];
sx q[2];
rz(0.21710795) q[2];
rz(1.0202967) q[3];
sx q[3];
rz(-1.3291357) q[3];
sx q[3];
rz(0.98062688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8751136) q[0];
sx q[0];
rz(-1.7784235) q[0];
sx q[0];
rz(-2.3734221) q[0];
rz(2.8504596) q[1];
sx q[1];
rz(-1.365064) q[1];
sx q[1];
rz(1.1870144) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465837) q[0];
sx q[0];
rz(-1.6542985) q[0];
sx q[0];
rz(-2.9299987) q[0];
x q[1];
rz(1.6441543) q[2];
sx q[2];
rz(-1.5268451) q[2];
sx q[2];
rz(-2.8188561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1817081) q[1];
sx q[1];
rz(-1.925029) q[1];
sx q[1];
rz(-1.6003057) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1891425) q[3];
sx q[3];
rz(-0.72126167) q[3];
sx q[3];
rz(-0.14474584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27176303) q[2];
sx q[2];
rz(-2.2972079) q[2];
sx q[2];
rz(-2.1738906) q[2];
rz(-2.9076231) q[3];
sx q[3];
rz(-2.440019) q[3];
sx q[3];
rz(-0.35681891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135603) q[0];
sx q[0];
rz(-1.5910633) q[0];
sx q[0];
rz(-2.0618942) q[0];
rz(0.30403852) q[1];
sx q[1];
rz(-1.1540776) q[1];
sx q[1];
rz(1.8255723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42342351) q[0];
sx q[0];
rz(-1.4470513) q[0];
sx q[0];
rz(-2.607671) q[0];
rz(-pi) q[1];
x q[1];
rz(1.448668) q[2];
sx q[2];
rz(-2.5425362) q[2];
sx q[2];
rz(2.0640822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58319891) q[1];
sx q[1];
rz(-2.2474562) q[1];
sx q[1];
rz(2.3469902) q[1];
x q[2];
rz(3.0791829) q[3];
sx q[3];
rz(-1.1452598) q[3];
sx q[3];
rz(-1.1264914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6198081) q[2];
sx q[2];
rz(-1.9079756) q[2];
sx q[2];
rz(3.0002777) q[2];
rz(-2.5610793) q[3];
sx q[3];
rz(-1.6655191) q[3];
sx q[3];
rz(0.22835246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2572131) q[0];
sx q[0];
rz(-2.3452106) q[0];
sx q[0];
rz(1.4759395) q[0];
rz(-2.2085704) q[1];
sx q[1];
rz(-2.4304183) q[1];
sx q[1];
rz(-1.3915541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2403657) q[0];
sx q[0];
rz(-2.3963442) q[0];
sx q[0];
rz(1.3979493) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7180842) q[2];
sx q[2];
rz(-1.4299586) q[2];
sx q[2];
rz(0.25885669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7793624) q[1];
sx q[1];
rz(-2.0915151) q[1];
sx q[1];
rz(0.034265072) q[1];
x q[2];
rz(0.33008195) q[3];
sx q[3];
rz(-2.0795855) q[3];
sx q[3];
rz(-1.1710081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2453904) q[2];
sx q[2];
rz(-2.7707477) q[2];
sx q[2];
rz(-2.89213) q[2];
rz(-1.4977411) q[3];
sx q[3];
rz(-1.4062873) q[3];
sx q[3];
rz(-0.41516414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6849218) q[0];
sx q[0];
rz(-2.6455854) q[0];
sx q[0];
rz(1.7556835) q[0];
rz(1.8798401) q[1];
sx q[1];
rz(-1.2077786) q[1];
sx q[1];
rz(-0.52880803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7700634) q[0];
sx q[0];
rz(-2.0160107) q[0];
sx q[0];
rz(3.1337409) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4086921) q[2];
sx q[2];
rz(-0.7846047) q[2];
sx q[2];
rz(-0.32151383) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4836135) q[1];
sx q[1];
rz(-0.93320642) q[1];
sx q[1];
rz(0.083483551) q[1];
rz(-pi) q[2];
rz(2.4566023) q[3];
sx q[3];
rz(-1.6719975) q[3];
sx q[3];
rz(-1.2736959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1888107) q[2];
sx q[2];
rz(-2.6691801) q[2];
sx q[2];
rz(-1.7702276) q[2];
rz(-2.93907) q[3];
sx q[3];
rz(-1.6731508) q[3];
sx q[3];
rz(-1.1892345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33223575) q[0];
sx q[0];
rz(-1.3487331) q[0];
sx q[0];
rz(0.15383823) q[0];
rz(0.73506749) q[1];
sx q[1];
rz(-2.0497597) q[1];
sx q[1];
rz(-0.14911266) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64407883) q[0];
sx q[0];
rz(-0.91877684) q[0];
sx q[0];
rz(0.10494167) q[0];
x q[1];
rz(2.9682069) q[2];
sx q[2];
rz(-1.9299704) q[2];
sx q[2];
rz(-2.2196291) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11483773) q[1];
sx q[1];
rz(-0.85875612) q[1];
sx q[1];
rz(0.65178443) q[1];
rz(-1.0668236) q[3];
sx q[3];
rz(-1.6503346) q[3];
sx q[3];
rz(-0.37802896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4601124) q[2];
sx q[2];
rz(-1.6296547) q[2];
sx q[2];
rz(-2.6849449) q[2];
rz(1.418142) q[3];
sx q[3];
rz(-0.8816312) q[3];
sx q[3];
rz(2.452623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.5099455) q[0];
sx q[0];
rz(-2.8896285) q[0];
sx q[0];
rz(3.0889567) q[0];
rz(-1.3098199) q[1];
sx q[1];
rz(-2.8927264) q[1];
sx q[1];
rz(1.2059258) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9249155) q[0];
sx q[0];
rz(-1.9572568) q[0];
sx q[0];
rz(-1.0189702) q[0];
rz(-2.7810026) q[2];
sx q[2];
rz(-2.5506488) q[2];
sx q[2];
rz(0.49345582) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0585021) q[1];
sx q[1];
rz(-1.7278226) q[1];
sx q[1];
rz(-2.462802) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8926335) q[3];
sx q[3];
rz(-0.6475237) q[3];
sx q[3];
rz(2.3237438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.84411088) q[2];
sx q[2];
rz(-1.7606807) q[2];
sx q[2];
rz(2.5816176) q[2];
rz(2.0848134) q[3];
sx q[3];
rz(-2.4323075) q[3];
sx q[3];
rz(2.3287676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63447222) q[0];
sx q[0];
rz(-1.3993323) q[0];
sx q[0];
rz(-1.5647474) q[0];
rz(-1.5722081) q[1];
sx q[1];
rz(-1.1610616) q[1];
sx q[1];
rz(-0.80894583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3472373) q[0];
sx q[0];
rz(-2.6449727) q[0];
sx q[0];
rz(-0.079023208) q[0];
rz(0.81124146) q[2];
sx q[2];
rz(-1.277207) q[2];
sx q[2];
rz(0.21101235) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2613649) q[1];
sx q[1];
rz(-2.5305809) q[1];
sx q[1];
rz(-1.5935807) q[1];
rz(-1.7478757) q[3];
sx q[3];
rz(-1.5054323) q[3];
sx q[3];
rz(-0.040249947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.481015) q[2];
sx q[2];
rz(-2.1889841) q[2];
sx q[2];
rz(3.0926404) q[2];
rz(1.281721) q[3];
sx q[3];
rz(-2.6661524) q[3];
sx q[3];
rz(1.6767282) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.319741) q[0];
sx q[0];
rz(-0.26092437) q[0];
sx q[0];
rz(-1.2224181) q[0];
rz(-1.6659196) q[1];
sx q[1];
rz(-1.2011352) q[1];
sx q[1];
rz(-0.45752057) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0973061) q[0];
sx q[0];
rz(-0.73839085) q[0];
sx q[0];
rz(-3.0074658) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7399212) q[2];
sx q[2];
rz(-1.0980415) q[2];
sx q[2];
rz(0.061701802) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4508707) q[1];
sx q[1];
rz(-1.9198981) q[1];
sx q[1];
rz(-1.3551177) q[1];
x q[2];
rz(1.1531257) q[3];
sx q[3];
rz(-0.92781298) q[3];
sx q[3];
rz(-2.2394799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7484625) q[2];
sx q[2];
rz(-2.8905383) q[2];
sx q[2];
rz(2.1012696) q[2];
rz(-2.6535502) q[3];
sx q[3];
rz(-1.4405684) q[3];
sx q[3];
rz(0.53781167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77448612) q[0];
sx q[0];
rz(-2.0699061) q[0];
sx q[0];
rz(0.62181428) q[0];
rz(1.294301) q[1];
sx q[1];
rz(-2.558567) q[1];
sx q[1];
rz(-0.8898215) q[1];
rz(2.1637259) q[2];
sx q[2];
rz(-1.3203096) q[2];
sx q[2];
rz(1.6537501) q[2];
rz(-1.9593976) q[3];
sx q[3];
rz(-0.91872707) q[3];
sx q[3];
rz(1.6891458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
