OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(1.1428042) q[0];
sx q[0];
rz(11.354843) q[0];
rz(-0.22663528) q[1];
sx q[1];
rz(-1.5770788) q[1];
sx q[1];
rz(0.29830631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5601215) q[0];
sx q[0];
rz(-1.6257964) q[0];
sx q[0];
rz(-2.4762857) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0478893) q[2];
sx q[2];
rz(-1.9022577) q[2];
sx q[2];
rz(-1.1790438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6023941) q[1];
sx q[1];
rz(-2.1089923) q[1];
sx q[1];
rz(2.5799275) q[1];
x q[2];
rz(1.467642) q[3];
sx q[3];
rz(-1.4872331) q[3];
sx q[3];
rz(2.8148357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(0.99386627) q[2];
rz(2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(-0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64269972) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(2.6699064) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90549201) q[0];
sx q[0];
rz(-0.24370757) q[0];
sx q[0];
rz(-2.7729176) q[0];
rz(-pi) q[1];
rz(-2.6366028) q[2];
sx q[2];
rz(-1.8171176) q[2];
sx q[2];
rz(1.9976975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7029611) q[1];
sx q[1];
rz(-1.303181) q[1];
sx q[1];
rz(-0.57499927) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7945812) q[3];
sx q[3];
rz(-0.61962485) q[3];
sx q[3];
rz(2.1686045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5919684) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(-0.96898752) q[2];
rz(2.5668868) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(-2.1000752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4330924) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(0.18181268) q[0];
rz(-2.0388942) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(-1.4556494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239379) q[0];
sx q[0];
rz(-1.0929937) q[0];
sx q[0];
rz(2.2295203) q[0];
rz(-pi) q[1];
rz(-0.74080148) q[2];
sx q[2];
rz(-1.2503137) q[2];
sx q[2];
rz(-1.6522811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94565832) q[1];
sx q[1];
rz(-0.12609005) q[1];
sx q[1];
rz(0.083421589) q[1];
rz(-0.62526838) q[3];
sx q[3];
rz(-1.4843974) q[3];
sx q[3];
rz(-2.3141935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2640947) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(2.1739615) q[2];
rz(2.4140221) q[3];
sx q[3];
rz(-1.2604159) q[3];
sx q[3];
rz(2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(-2.0137285) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.9086054) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9905332) q[0];
sx q[0];
rz(-2.6418243) q[0];
sx q[0];
rz(0.14453669) q[0];
x q[1];
rz(-2.1190676) q[2];
sx q[2];
rz(-2.0892482) q[2];
sx q[2];
rz(-0.099629121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9150881) q[1];
sx q[1];
rz(-2.9595032) q[1];
sx q[1];
rz(-1.7712797) q[1];
x q[2];
rz(-0.87807699) q[3];
sx q[3];
rz(-1.746776) q[3];
sx q[3];
rz(0.39719492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91810742) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(-2.0984086) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(-1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84213132) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(2.2498851) q[0];
rz(-1.8978329) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(-2.9290501) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1067057) q[0];
sx q[0];
rz(-1.5910774) q[0];
sx q[0];
rz(-0.01625343) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5151305) q[2];
sx q[2];
rz(-2.5362483) q[2];
sx q[2];
rz(0.84673131) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7674539) q[1];
sx q[1];
rz(-0.98857388) q[1];
sx q[1];
rz(0.9064845) q[1];
rz(2.6796474) q[3];
sx q[3];
rz(-2.2432703) q[3];
sx q[3];
rz(-0.40601054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1084958) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(-0.5212211) q[2];
rz(1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(1.2683755) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8591156) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(2.916472) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(0.37757847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9281884) q[0];
sx q[0];
rz(-1.3047991) q[0];
sx q[0];
rz(-0.022797419) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16212459) q[2];
sx q[2];
rz(-2.8505278) q[2];
sx q[2];
rz(1.383701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4695417) q[1];
sx q[1];
rz(-2.6968323) q[1];
sx q[1];
rz(1.3252392) q[1];
rz(-0.029929786) q[3];
sx q[3];
rz(-0.79513351) q[3];
sx q[3];
rz(-3.0628169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(-2.6605576) q[2];
rz(0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(-0.19454923) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(2.887168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1628111) q[0];
sx q[0];
rz(-0.33919507) q[0];
sx q[0];
rz(1.9041054) q[0];
x q[1];
rz(2.9872586) q[2];
sx q[2];
rz(-0.96417226) q[2];
sx q[2];
rz(2.2301205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3212657) q[1];
sx q[1];
rz(-1.1202381) q[1];
sx q[1];
rz(-1.8727559) q[1];
x q[2];
rz(0.45215996) q[3];
sx q[3];
rz(-1.9981355) q[3];
sx q[3];
rz(2.6219581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.97757942) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(-1.8257726) q[2];
rz(2.2655462) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(-2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(2.3176106) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(-1.5751858) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67714171) q[0];
sx q[0];
rz(-1.8832708) q[0];
sx q[0];
rz(0.72014767) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1112061) q[2];
sx q[2];
rz(-1.724913) q[2];
sx q[2];
rz(2.266778) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52289256) q[1];
sx q[1];
rz(-0.52378264) q[1];
sx q[1];
rz(2.3431542) q[1];
rz(-pi) q[2];
rz(2.4616562) q[3];
sx q[3];
rz(-0.61938647) q[3];
sx q[3];
rz(0.030127545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(0.27754647) q[2];
rz(1.4510441) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9797416) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.45654) q[0];
rz(0.62943554) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(-2.004752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0946227) q[0];
sx q[0];
rz(-1.2828865) q[0];
sx q[0];
rz(2.1043491) q[0];
x q[1];
rz(-0.37461899) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(-1.7916726) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7965664) q[1];
sx q[1];
rz(-1.6654286) q[1];
sx q[1];
rz(1.8887397) q[1];
rz(-pi) q[2];
rz(1.9433446) q[3];
sx q[3];
rz(-2.2714943) q[3];
sx q[3];
rz(-0.31853279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(-2.1006404) q[2];
rz(-3.1395636) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(2.8295529) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390389) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(2.1955406) q[0];
rz(-2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(-2.5295703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8004868) q[0];
sx q[0];
rz(-2.2414811) q[0];
sx q[0];
rz(2.1675046) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34531784) q[2];
sx q[2];
rz(-1.7974263) q[2];
sx q[2];
rz(2.741284) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7358688) q[1];
sx q[1];
rz(-1.6940261) q[1];
sx q[1];
rz(-0.65995364) q[1];
rz(0.53491433) q[3];
sx q[3];
rz(-1.2423008) q[3];
sx q[3];
rz(-1.9096149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(-0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(2.4408834) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54031298) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(2.3868949) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-0.66773141) q[2];
sx q[2];
rz(-2.3420391) q[2];
sx q[2];
rz(2.5426368) q[2];
rz(1.3238974) q[3];
sx q[3];
rz(-0.3331475) q[3];
sx q[3];
rz(-3.0039136) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
