OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(-0.20110826) q[0];
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(2.5147658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7456822) q[0];
sx q[0];
rz(-0.22326176) q[0];
sx q[0];
rz(-2.1366871) q[0];
rz(-pi) q[1];
rz(-1.8148412) q[2];
sx q[2];
rz(-1.913117) q[2];
sx q[2];
rz(-2.2133881) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99170291) q[1];
sx q[1];
rz(-2.457379) q[1];
sx q[1];
rz(1.0690772) q[1];
rz(1.2960035) q[3];
sx q[3];
rz(-1.7405675) q[3];
sx q[3];
rz(-2.3613289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(-2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(0.75479341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24793808) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(0.98518103) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(-3.0342297) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3283008) q[0];
sx q[0];
rz(-2.11587) q[0];
sx q[0];
rz(-2.6585447) q[0];
x q[1];
rz(2.7254843) q[2];
sx q[2];
rz(-0.62972087) q[2];
sx q[2];
rz(-2.1925418) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3423691) q[1];
sx q[1];
rz(-0.93042004) q[1];
sx q[1];
rz(0.46701042) q[1];
rz(-1.4025877) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5169516) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(-0.66169935) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702883) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(3.1012428) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(-1.0823762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45535116) q[0];
sx q[0];
rz(-2.4903114) q[0];
sx q[0];
rz(2.5097646) q[0];
rz(-1.1788549) q[2];
sx q[2];
rz(-2.760663) q[2];
sx q[2];
rz(0.71289635) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7554453) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(-0.2281245) q[1];
rz(0.97614395) q[3];
sx q[3];
rz(-2.757132) q[3];
sx q[3];
rz(2.773075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(2.5877) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937113) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(-0.035382263) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(-0.48809537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29823829) q[0];
sx q[0];
rz(-2.6697864) q[0];
sx q[0];
rz(1.7422373) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1610291) q[2];
sx q[2];
rz(-1.0733986) q[2];
sx q[2];
rz(2.892708) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6285703) q[1];
sx q[1];
rz(-1.7254618) q[1];
sx q[1];
rz(-0.13048529) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.475708) q[3];
sx q[3];
rz(-2.4443691) q[3];
sx q[3];
rz(2.9149885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(1.1821702) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(-1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(-0.3702634) q[0];
rz(0.26369357) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(0.19651861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.927658) q[0];
sx q[0];
rz(-1.2688583) q[0];
sx q[0];
rz(-1.2769804) q[0];
rz(-pi) q[1];
rz(-1.1097124) q[2];
sx q[2];
rz(-1.5421252) q[2];
sx q[2];
rz(-1.2210786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4938441) q[1];
sx q[1];
rz(-1.7462247) q[1];
sx q[1];
rz(-2.8791048) q[1];
rz(-pi) q[2];
rz(1.7408095) q[3];
sx q[3];
rz(-1.8432872) q[3];
sx q[3];
rz(-1.221399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1281517) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(-0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(2.939558) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0840543) q[0];
sx q[0];
rz(-1.6469427) q[0];
sx q[0];
rz(2.8094859) q[0];
rz(-pi) q[1];
rz(-1.2932111) q[2];
sx q[2];
rz(-2.7536014) q[2];
sx q[2];
rz(-2.9953621) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.85642568) q[1];
sx q[1];
rz(-1.3060096) q[1];
sx q[1];
rz(1.9734288) q[1];
x q[2];
rz(2.6447547) q[3];
sx q[3];
rz(-2.3416069) q[3];
sx q[3];
rz(2.5132688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(1.0874282) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(-1.6585763) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6120646) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(-2.0102603) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-0.21025118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41482571) q[0];
sx q[0];
rz(-1.6370156) q[0];
sx q[0];
rz(1.4365175) q[0];
rz(-pi) q[1];
rz(-2.7619751) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(-1.1773674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6107632) q[1];
sx q[1];
rz(-1.3794583) q[1];
sx q[1];
rz(0.32236871) q[1];
rz(-pi) q[2];
rz(-0.57742124) q[3];
sx q[3];
rz(-0.93772674) q[3];
sx q[3];
rz(-3.0949743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.38114) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(-0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(2.6953221) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(2.1648724) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2864838) q[0];
sx q[0];
rz(-2.6628462) q[0];
sx q[0];
rz(0.32740645) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.080539695) q[2];
sx q[2];
rz(-1.409568) q[2];
sx q[2];
rz(-1.5921519) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50423065) q[1];
sx q[1];
rz(-2.0081875) q[1];
sx q[1];
rz(0.38423844) q[1];
rz(-2.1920491) q[3];
sx q[3];
rz(-0.81406677) q[3];
sx q[3];
rz(0.66063389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(-2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66822806) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(-0.92064944) q[0];
rz(-3.0905511) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(-2.2424973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5024922) q[0];
sx q[0];
rz(-0.99196767) q[0];
sx q[0];
rz(2.6274908) q[0];
x q[1];
rz(-2.8092438) q[2];
sx q[2];
rz(-2.1045448) q[2];
sx q[2];
rz(0.74408434) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2194849) q[1];
sx q[1];
rz(-1.2412294) q[1];
sx q[1];
rz(3.0199514) q[1];
rz(1.2006239) q[3];
sx q[3];
rz(-1.5705974) q[3];
sx q[3];
rz(1.9716594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1961394) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(-2.3032522) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33912441) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(0.71240187) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-2.7446279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1966232) q[0];
sx q[0];
rz(-2.0485749) q[0];
sx q[0];
rz(-2.1150132) q[0];
rz(1.9372349) q[2];
sx q[2];
rz(-1.8574742) q[2];
sx q[2];
rz(2.3057126) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9576479) q[1];
sx q[1];
rz(-2.7937104) q[1];
sx q[1];
rz(1.3089048) q[1];
rz(-pi) q[2];
rz(-2.0685365) q[3];
sx q[3];
rz(-1.2101678) q[3];
sx q[3];
rz(-1.3090759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1931856) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5554572) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(1.4861698) q[2];
sx q[2];
rz(-1.463269) q[2];
sx q[2];
rz(-0.29381897) q[2];
rz(1.312064) q[3];
sx q[3];
rz(-1.4553634) q[3];
sx q[3];
rz(2.8298557) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
