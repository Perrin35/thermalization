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
rz(-2.8432863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053781833) q[0];
sx q[0];
rz(-2.2349173) q[0];
sx q[0];
rz(-1.640663) q[0];
x q[1];
rz(-1.9036129) q[2];
sx q[2];
rz(-1.4822072) q[2];
sx q[2];
rz(0.36117902) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71516192) q[1];
sx q[1];
rz(-0.75725812) q[1];
sx q[1];
rz(-2.2992579) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0575849) q[3];
sx q[3];
rz(-1.4680032) q[3];
sx q[3];
rz(-1.2526797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(2.1477264) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(-2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(0.64269972) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(-0.018218536) q[0];
rz(-2.3253564) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(2.6699064) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6149711) q[0];
sx q[0];
rz(-1.3437628) q[0];
sx q[0];
rz(1.4814266) q[0];
x q[1];
rz(0.50498982) q[2];
sx q[2];
rz(-1.3244751) q[2];
sx q[2];
rz(-1.9976975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3968518) q[1];
sx q[1];
rz(-2.5137915) q[1];
sx q[1];
rz(0.46698924) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3470115) q[3];
sx q[3];
rz(-0.61962485) q[3];
sx q[3];
rz(-0.97298813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5919684) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(2.1726051) q[2];
rz(-0.5747059) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(2.95978) q[0];
rz(1.1026985) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(1.6859432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2244814) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(-2.2729421) q[0];
x q[1];
rz(0.74080148) q[2];
sx q[2];
rz(-1.2503137) q[2];
sx q[2];
rz(-1.4893116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0297444) q[1];
sx q[1];
rz(-1.6964456) q[1];
sx q[1];
rz(-1.5602342) q[1];
rz(-1.4643747) q[3];
sx q[3];
rz(-0.94821804) q[3];
sx q[3];
rz(-0.6811844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2640947) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(-0.96763119) q[2];
rz(2.4140221) q[3];
sx q[3];
rz(-1.2604159) q[3];
sx q[3];
rz(-0.23770604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(-0.63823429) q[0];
rz(2.0137285) q[1];
sx q[1];
rz(-2.3141839) q[1];
sx q[1];
rz(1.9086054) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1283135) q[0];
sx q[0];
rz(-1.0767125) q[0];
sx q[0];
rz(-1.6492776) q[0];
x q[1];
rz(1.0225251) q[2];
sx q[2];
rz(-1.0523445) q[2];
sx q[2];
rz(-3.0419635) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4302664) q[1];
sx q[1];
rz(-1.3923936) q[1];
sx q[1];
rz(-3.1049411) q[1];
rz(-pi) q[2];
rz(-2.9144985) q[3];
sx q[3];
rz(-0.89082754) q[3];
sx q[3];
rz(2.1122776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91810742) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(0.81400648) q[2];
rz(1.043184) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(-1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84213132) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(-2.2498851) q[0];
rz(1.2437598) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(-2.9290501) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5355797) q[0];
sx q[0];
rz(-1.5870464) q[0];
sx q[0];
rz(1.5505126) q[0];
rz(-1.5151305) q[2];
sx q[2];
rz(-2.5362483) q[2];
sx q[2];
rz(-2.2948613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37413874) q[1];
sx q[1];
rz(-0.98857388) q[1];
sx q[1];
rz(-0.9064845) q[1];
rz(-pi) q[2];
rz(1.0605293) q[3];
sx q[3];
rz(-2.3465996) q[3];
sx q[3];
rz(-2.0612962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1084958) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(2.6203716) q[2];
rz(1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(-1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824771) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(2.916472) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(2.7640142) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8416653) q[0];
sx q[0];
rz(-0.26694926) q[0];
sx q[0];
rz(1.4873234) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6191143) q[2];
sx q[2];
rz(-1.2836576) q[2];
sx q[2];
rz(-1.552812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0203591) q[1];
sx q[1];
rz(-1.4660144) q[1];
sx q[1];
rz(-1.1377513) q[1];
x q[2];
rz(1.6013006) q[3];
sx q[3];
rz(-0.77611938) q[3];
sx q[3];
rz(-3.0200849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(0.48103508) q[2];
rz(-2.7379819) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0525381) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(0.19454923) q[0];
rz(2.9220707) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(0.25442466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90760224) q[0];
sx q[0];
rz(-1.4617209) q[0];
sx q[0];
rz(1.8926189) q[0];
rz(-2.1830325) q[2];
sx q[2];
rz(-1.6974291) q[2];
sx q[2];
rz(0.74778344) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82032694) q[1];
sx q[1];
rz(-2.0213545) q[1];
sx q[1];
rz(-1.2688368) q[1];
rz(-0.80612225) q[3];
sx q[3];
rz(-0.61180173) q[3];
sx q[3];
rz(-2.7968189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.97757942) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(1.3158201) q[2];
rz(2.2655462) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(1.0036489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106237) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(2.3176106) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(-1.5664068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9109089) q[0];
sx q[0];
rz(-2.3678603) q[0];
sx q[0];
rz(-0.45549972) q[0];
rz(-pi) q[1];
rz(2.9623716) q[2];
sx q[2];
rz(-2.1041098) q[2];
sx q[2];
rz(-2.3537677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6187001) q[1];
sx q[1];
rz(-0.52378264) q[1];
sx q[1];
rz(-0.79843847) q[1];
rz(-2.4616562) q[3];
sx q[3];
rz(-0.61938647) q[3];
sx q[3];
rz(-0.030127545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87551293) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(-1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16185109) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(1.6850527) q[0];
rz(-0.62943554) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(2.004752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0946227) q[0];
sx q[0];
rz(-1.8587062) q[0];
sx q[0];
rz(-1.0372435) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37461899) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(1.7916726) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5053619) q[1];
sx q[1];
rz(-0.33126918) q[1];
sx q[1];
rz(-1.865571) q[1];
x q[2];
rz(-0.40739079) q[3];
sx q[3];
rz(-2.3630777) q[3];
sx q[3];
rz(0.22637573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(1.0409522) q[2];
rz(0.0020290931) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(2.229915) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(0.61202234) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6529022) q[0];
sx q[0];
rz(-2.2757747) q[0];
sx q[0];
rz(-0.61625723) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5435796) q[2];
sx q[2];
rz(-0.41053718) q[2];
sx q[2];
rz(1.41278) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0079841) q[1];
sx q[1];
rz(-2.4719279) q[1];
sx q[1];
rz(0.19934166) q[1];
rz(-2.5522334) q[3];
sx q[3];
rz(-0.61925626) q[3];
sx q[3];
rz(-2.9818231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(-2.4882312) q[2];
rz(2.7907794) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(0.70070926) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6012797) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(-0.75469771) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(0.67946658) q[2];
sx q[2];
rz(-1.1107399) q[2];
sx q[2];
rz(-2.6723292) q[2];
rz(3.0572206) q[3];
sx q[3];
rz(-1.8934688) q[3];
sx q[3];
rz(-2.7432751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
