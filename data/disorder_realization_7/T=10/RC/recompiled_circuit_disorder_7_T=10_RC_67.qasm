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
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7456822) q[0];
sx q[0];
rz(-2.9183309) q[0];
sx q[0];
rz(-2.1366871) q[0];
rz(-1.8148412) q[2];
sx q[2];
rz(-1.913117) q[2];
sx q[2];
rz(-2.2133881) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.99170291) q[1];
sx q[1];
rz(-0.68421364) q[1];
sx q[1];
rz(-2.0725155) q[1];
rz(2.1342282) q[3];
sx q[3];
rz(-2.8197188) q[3];
sx q[3];
rz(-2.8909627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25519249) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.2935151) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24793808) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(2.1564116) q[0];
rz(2.7424116) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(-0.10736297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1194612) q[0];
sx q[0];
rz(-0.71170002) q[0];
sx q[0];
rz(2.2244781) q[0];
x q[1];
rz(-1.2843578) q[2];
sx q[2];
rz(-1.0019433) q[2];
sx q[2];
rz(1.4494277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0774035) q[1];
sx q[1];
rz(-1.2014376) q[1];
sx q[1];
rz(2.2662152) q[1];
rz(1.739005) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(-0.66169935) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645638) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-3.1012428) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(1.0823762) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4989935) q[0];
sx q[0];
rz(-1.2046308) q[0];
sx q[0];
rz(-2.5901592) q[0];
rz(-pi) q[1];
rz(2.9897887) q[2];
sx q[2];
rz(-1.2200583) q[2];
sx q[2];
rz(-1.1317859) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(-0.2281245) q[1];
rz(-1.8941746) q[3];
sx q[3];
rz(-1.3591027) q[3];
sx q[3];
rz(2.4993103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0064156) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(-2.5877) q[2];
rz(-1.4931549) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(2.5890787) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(3.1062104) q[0];
rz(-0.9961876) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-2.6534973) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0220538) q[0];
sx q[0];
rz(-1.6484123) q[0];
sx q[0];
rz(1.104943) q[0];
rz(-0.79779063) q[2];
sx q[2];
rz(-2.3893223) q[2];
sx q[2];
rz(0.70300245) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51302233) q[1];
sx q[1];
rz(-1.7254618) q[1];
sx q[1];
rz(3.0111074) q[1];
x q[2];
rz(-0.66588464) q[3];
sx q[3];
rz(-2.4443691) q[3];
sx q[3];
rz(0.22660412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35486832) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(-1.1821702) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(-1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.5421211) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(2.7713293) q[0];
rz(-0.26369357) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(2.945074) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949961) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(0.3145991) q[0];
x q[1];
rz(-1.5064266) q[2];
sx q[2];
rz(-2.6796821) q[2];
sx q[2];
rz(-0.40735746) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4938441) q[1];
sx q[1];
rz(-1.3953679) q[1];
sx q[1];
rz(-2.8791048) q[1];
x q[2];
rz(0.27627857) q[3];
sx q[3];
rz(-1.4071138) q[3];
sx q[3];
rz(0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1281517) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(0.91709843) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(0.20203461) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0575384) q[0];
sx q[0];
rz(-1.6469427) q[0];
sx q[0];
rz(0.3321068) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2932111) q[2];
sx q[2];
rz(-0.38799122) q[2];
sx q[2];
rz(-0.1462305) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.85642568) q[1];
sx q[1];
rz(-1.3060096) q[1];
sx q[1];
rz(1.9734288) q[1];
rz(2.0270258) q[3];
sx q[3];
rz(-2.2531415) q[3];
sx q[3];
rz(1.8519459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1918734) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(-2.2992772) q[2];
rz(-2.0541644) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(-1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52952805) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(2.1202309) q[0];
rz(2.0102603) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(2.9313415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1649095) q[0];
sx q[0];
rz(-1.4368136) q[0];
sx q[0];
rz(3.0747736) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7619751) q[2];
sx q[2];
rz(-0.175975) q[2];
sx q[2];
rz(1.9642252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0381895) q[1];
sx q[1];
rz(-1.2545171) q[1];
sx q[1];
rz(-1.3693387) q[1];
x q[2];
rz(-0.93123575) q[3];
sx q[3];
rz(-0.82914549) q[3];
sx q[3];
rz(0.88013807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82688275) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(1.38114) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(-0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(-0.44627055) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(-2.1648724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92102345) q[0];
sx q[0];
rz(-1.1194293) q[0];
sx q[0];
rz(-1.736182) q[0];
rz(-pi) q[1];
rz(1.7325399) q[2];
sx q[2];
rz(-1.6502893) q[2];
sx q[2];
rz(0.0083991945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8969438) q[1];
sx q[1];
rz(-1.2243425) q[1];
sx q[1];
rz(1.1036554) q[1];
rz(-pi) q[2];
rz(-0.94954357) q[3];
sx q[3];
rz(-0.81406677) q[3];
sx q[3];
rz(2.4809588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66822806) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(-2.2209432) q[0];
rz(-0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(-0.89909536) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2313401) q[0];
sx q[0];
rz(-1.1465342) q[0];
sx q[0];
rz(2.2146241) q[0];
rz(-pi) q[1];
rz(-1.012072) q[2];
sx q[2];
rz(-1.8555102) q[2];
sx q[2];
rz(0.65288359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5324459) q[1];
sx q[1];
rz(-1.4557314) q[1];
sx q[1];
rz(-1.902641) q[1];
x q[2];
rz(1.5702463) q[3];
sx q[3];
rz(-0.37017248) q[3];
sx q[3];
rz(0.40035029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(0.10458874) q[2];
rz(0.83834046) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33912441) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(2.0237645) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-0.39696473) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94496942) q[0];
sx q[0];
rz(-2.0485749) q[0];
sx q[0];
rz(-2.1150132) q[0];
x q[1];
rz(1.2043578) q[2];
sx q[2];
rz(-1.2841184) q[2];
sx q[2];
rz(-0.8358801) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46170235) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(-0.093613503) q[1];
rz(2.7361717) q[3];
sx q[3];
rz(-1.1076895) q[3];
sx q[3];
rz(-0.4511569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-3.0449384) q[2];
rz(1.4139253) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5554572) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(1.6554228) q[2];
sx q[2];
rz(-1.6783236) q[2];
sx q[2];
rz(2.8477737) q[2];
rz(-1.312064) q[3];
sx q[3];
rz(-1.6862292) q[3];
sx q[3];
rz(-0.31173691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];