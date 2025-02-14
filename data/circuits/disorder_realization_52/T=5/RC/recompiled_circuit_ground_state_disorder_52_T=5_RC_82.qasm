OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0527394) q[0];
sx q[0];
rz(-2.7272447) q[0];
sx q[0];
rz(-2.156884) q[0];
rz(2.2766278) q[1];
sx q[1];
rz(-2.2948761) q[1];
sx q[1];
rz(2.4196978) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5295204) q[0];
sx q[0];
rz(-1.9871166) q[0];
sx q[0];
rz(0.4952391) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4059695) q[2];
sx q[2];
rz(-1.3072701) q[2];
sx q[2];
rz(-1.6014388) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5492603) q[1];
sx q[1];
rz(-1.6417786) q[1];
sx q[1];
rz(-3.1123726) q[1];
x q[2];
rz(-1.1837051) q[3];
sx q[3];
rz(-1.4287881) q[3];
sx q[3];
rz(-1.3660688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3856753) q[2];
sx q[2];
rz(-2.6188681) q[2];
sx q[2];
rz(-0.44504607) q[2];
rz(-1.2612777) q[3];
sx q[3];
rz(-1.4538366) q[3];
sx q[3];
rz(1.5104347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5234914) q[0];
sx q[0];
rz(-1.0502676) q[0];
sx q[0];
rz(-0.65446788) q[0];
rz(2.0596313) q[1];
sx q[1];
rz(-1.315821) q[1];
sx q[1];
rz(1.4644324) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0102756) q[0];
sx q[0];
rz(-0.47476381) q[0];
sx q[0];
rz(-1.4608338) q[0];
rz(-pi) q[1];
rz(-3.0686321) q[2];
sx q[2];
rz(-2.1734997) q[2];
sx q[2];
rz(-0.51034865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46772114) q[1];
sx q[1];
rz(-2.4003211) q[1];
sx q[1];
rz(0.060972496) q[1];
rz(2.6707021) q[3];
sx q[3];
rz(-1.8857368) q[3];
sx q[3];
rz(-2.5641914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64708853) q[2];
sx q[2];
rz(-0.096179811) q[2];
sx q[2];
rz(0.75867009) q[2];
rz(1.8132973) q[3];
sx q[3];
rz(-2.0601065) q[3];
sx q[3];
rz(1.7648511) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7741622) q[0];
sx q[0];
rz(-1.6663015) q[0];
sx q[0];
rz(-0.037516315) q[0];
rz(1.0388177) q[1];
sx q[1];
rz(-0.33375868) q[1];
sx q[1];
rz(-2.2822101) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0276703) q[0];
sx q[0];
rz(-1.3827208) q[0];
sx q[0];
rz(-2.5187224) q[0];
x q[1];
rz(-0.1350624) q[2];
sx q[2];
rz(-2.3762868) q[2];
sx q[2];
rz(-1.0926334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9199279) q[1];
sx q[1];
rz(-1.2499735) q[1];
sx q[1];
rz(-2.0179835) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5435197) q[3];
sx q[3];
rz(-1.155164) q[3];
sx q[3];
rz(1.9126836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3694156) q[2];
sx q[2];
rz(-0.93537664) q[2];
sx q[2];
rz(0.76425648) q[2];
rz(0.94591004) q[3];
sx q[3];
rz(-1.2063682) q[3];
sx q[3];
rz(0.8980155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94364828) q[0];
sx q[0];
rz(-2.2479489) q[0];
sx q[0];
rz(-1.8290895) q[0];
rz(-2.8370044) q[1];
sx q[1];
rz(-2.1423788) q[1];
sx q[1];
rz(-2.8056858) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6524624) q[0];
sx q[0];
rz(-1.5565598) q[0];
sx q[0];
rz(1.4984309) q[0];
rz(-0.96520378) q[2];
sx q[2];
rz(-1.9450257) q[2];
sx q[2];
rz(2.7363914) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87631851) q[1];
sx q[1];
rz(-0.66702038) q[1];
sx q[1];
rz(-0.69885175) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8943039) q[3];
sx q[3];
rz(-0.22824057) q[3];
sx q[3];
rz(0.3397371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37561068) q[2];
sx q[2];
rz(-0.42372647) q[2];
sx q[2];
rz(1.3099526) q[2];
rz(1.8464108) q[3];
sx q[3];
rz(-2.3056307) q[3];
sx q[3];
rz(-1.0999934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.497371) q[0];
sx q[0];
rz(-0.27232429) q[0];
sx q[0];
rz(-0.8154794) q[0];
rz(0.71715912) q[1];
sx q[1];
rz(-1.4445137) q[1];
sx q[1];
rz(1.1517634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64989728) q[0];
sx q[0];
rz(-2.5743352) q[0];
sx q[0];
rz(-0.41383608) q[0];
x q[1];
rz(-0.30381153) q[2];
sx q[2];
rz(-1.3388435) q[2];
sx q[2];
rz(-1.0396921) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6155225) q[1];
sx q[1];
rz(-0.87837362) q[1];
sx q[1];
rz(0.39151616) q[1];
rz(-0.30188876) q[3];
sx q[3];
rz(-1.535901) q[3];
sx q[3];
rz(-0.85012415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0723476) q[2];
sx q[2];
rz(-2.3657511) q[2];
sx q[2];
rz(-1.5838712) q[2];
rz(0.97258687) q[3];
sx q[3];
rz(-1.2310622) q[3];
sx q[3];
rz(-1.7044273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.4438181) q[0];
sx q[0];
rz(-0.2727209) q[0];
sx q[0];
rz(0.66993237) q[0];
rz(2.2386235) q[1];
sx q[1];
rz(-2.4882856) q[1];
sx q[1];
rz(-0.088937581) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020821559) q[0];
sx q[0];
rz(-2.1820314) q[0];
sx q[0];
rz(-0.24586785) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3207664) q[2];
sx q[2];
rz(-2.2019115) q[2];
sx q[2];
rz(-0.64455143) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.34076128) q[1];
sx q[1];
rz(-0.25909875) q[1];
sx q[1];
rz(1.8482261) q[1];
x q[2];
rz(0.1286147) q[3];
sx q[3];
rz(-2.0484784) q[3];
sx q[3];
rz(2.0292119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4778121) q[2];
sx q[2];
rz(-2.0974443) q[2];
sx q[2];
rz(-3.1204379) q[2];
rz(-2.2145005) q[3];
sx q[3];
rz(-1.6095716) q[3];
sx q[3];
rz(-2.482614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053719036) q[0];
sx q[0];
rz(-1.3608195) q[0];
sx q[0];
rz(-1.1370283) q[0];
rz(-2.673705) q[1];
sx q[1];
rz(-1.7158022) q[1];
sx q[1];
rz(-0.55380026) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1374986) q[0];
sx q[0];
rz(-1.8589673) q[0];
sx q[0];
rz(-2.1779446) q[0];
rz(2.0050029) q[2];
sx q[2];
rz(-1.2469262) q[2];
sx q[2];
rz(-2.6153274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1688331) q[1];
sx q[1];
rz(-1.0565896) q[1];
sx q[1];
rz(2.7943816) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6662237) q[3];
sx q[3];
rz(-0.57668873) q[3];
sx q[3];
rz(-0.63729034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0566473) q[2];
sx q[2];
rz(-0.84359303) q[2];
sx q[2];
rz(2.4455369) q[2];
rz(2.7737235) q[3];
sx q[3];
rz(-1.1080247) q[3];
sx q[3];
rz(-2.8554816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87357658) q[0];
sx q[0];
rz(-1.3883256) q[0];
sx q[0];
rz(-2.3263113) q[0];
rz(-2.6878327) q[1];
sx q[1];
rz(-0.90180698) q[1];
sx q[1];
rz(0.59660965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2011947) q[0];
sx q[0];
rz(-2.0773281) q[0];
sx q[0];
rz(-0.80123637) q[0];
x q[1];
rz(-1.6758133) q[2];
sx q[2];
rz(-2.3987282) q[2];
sx q[2];
rz(-0.32203963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4441159) q[1];
sx q[1];
rz(-1.2713065) q[1];
sx q[1];
rz(-1.0018574) q[1];
rz(0.96489789) q[3];
sx q[3];
rz(-2.4526074) q[3];
sx q[3];
rz(2.8953972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0348908) q[2];
sx q[2];
rz(-1.7009578) q[2];
sx q[2];
rz(-1.9218669) q[2];
rz(1.8969511) q[3];
sx q[3];
rz(-1.605426) q[3];
sx q[3];
rz(-2.9419148) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.225746) q[0];
sx q[0];
rz(-2.2016278) q[0];
sx q[0];
rz(-1.4858656) q[0];
rz(2.0303717) q[1];
sx q[1];
rz(-1.8467555) q[1];
sx q[1];
rz(-1.3908386) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6538453) q[0];
sx q[0];
rz(-1.5731205) q[0];
sx q[0];
rz(2.8043158) q[0];
rz(-pi) q[1];
rz(-2.7971091) q[2];
sx q[2];
rz(-1.6793669) q[2];
sx q[2];
rz(2.6950633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35090205) q[1];
sx q[1];
rz(-1.6054253) q[1];
sx q[1];
rz(1.778855) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2038571) q[3];
sx q[3];
rz(-1.2436359) q[3];
sx q[3];
rz(-3.0373472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86604649) q[2];
sx q[2];
rz(-1.486843) q[2];
sx q[2];
rz(-0.1869959) q[2];
rz(-2.9343119) q[3];
sx q[3];
rz(-0.63396251) q[3];
sx q[3];
rz(-1.0435587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.2279376) q[0];
sx q[0];
rz(-0.74988237) q[0];
sx q[0];
rz(0.0023181152) q[0];
rz(-1.2222611) q[1];
sx q[1];
rz(-1.5536676) q[1];
sx q[1];
rz(1.7244171) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.173439) q[0];
sx q[0];
rz(-1.7002956) q[0];
sx q[0];
rz(0.20682516) q[0];
x q[1];
rz(1.8882636) q[2];
sx q[2];
rz(-2.0892882) q[2];
sx q[2];
rz(-0.41050342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.05310381) q[1];
sx q[1];
rz(-1.3657059) q[1];
sx q[1];
rz(-3.0916832) q[1];
rz(-pi) q[2];
rz(0.98448344) q[3];
sx q[3];
rz(-1.7147736) q[3];
sx q[3];
rz(1.3534304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8861387) q[2];
sx q[2];
rz(-1.6243434) q[2];
sx q[2];
rz(0.21190602) q[2];
rz(-0.40031561) q[3];
sx q[3];
rz(-2.2369657) q[3];
sx q[3];
rz(-1.8109842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0935852) q[0];
sx q[0];
rz(-1.2496017) q[0];
sx q[0];
rz(-1.6954419) q[0];
rz(2.6493337) q[1];
sx q[1];
rz(-1.4795563) q[1];
sx q[1];
rz(1.5835887) q[1];
rz(2.3401101) q[2];
sx q[2];
rz(-2.8469606) q[2];
sx q[2];
rz(-1.4599232) q[2];
rz(2.0782804) q[3];
sx q[3];
rz(-1.1125269) q[3];
sx q[3];
rz(-0.0081882523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
