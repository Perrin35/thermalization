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
rz(-2.5459557) q[0];
sx q[0];
rz(-0.41002265) q[0];
sx q[0];
rz(-2.3469143) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(-0.95212189) q[1];
sx q[1];
rz(-1.9917537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4021989) q[0];
sx q[0];
rz(-1.5350685) q[0];
sx q[0];
rz(1.8648007) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8914521) q[2];
sx q[2];
rz(-2.3806664) q[2];
sx q[2];
rz(0.93283949) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9147891) q[1];
sx q[1];
rz(-2.7046596) q[1];
sx q[1];
rz(2.7257257) q[1];
rz(1.8273699) q[3];
sx q[3];
rz(-1.1169516) q[3];
sx q[3];
rz(-0.75405771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4475693) q[2];
sx q[2];
rz(-2.1489216) q[2];
sx q[2];
rz(-0.64219323) q[2];
rz(-1.0688952) q[3];
sx q[3];
rz(-1.5725719) q[3];
sx q[3];
rz(1.4890891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25817961) q[0];
sx q[0];
rz(-3.1386107) q[0];
sx q[0];
rz(-0.51134837) q[0];
rz(1.5072352) q[1];
sx q[1];
rz(-1.8664482) q[1];
sx q[1];
rz(1.3828329) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94039932) q[0];
sx q[0];
rz(-0.44006881) q[0];
sx q[0];
rz(2.0725193) q[0];
rz(-pi) q[1];
rz(-1.9963919) q[2];
sx q[2];
rz(-0.71768242) q[2];
sx q[2];
rz(-1.8459783) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1630711) q[1];
sx q[1];
rz(-0.51632612) q[1];
sx q[1];
rz(2.2341245) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6182943) q[3];
sx q[3];
rz(-1.4804258) q[3];
sx q[3];
rz(-0.096631526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6766659) q[2];
sx q[2];
rz(-1.2049462) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(-2.0788705) q[3];
sx q[3];
rz(-1.1954185) q[3];
sx q[3];
rz(-3.0265871) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17387834) q[0];
sx q[0];
rz(-2.102484) q[0];
sx q[0];
rz(-0.1097196) q[0];
rz(0.32490718) q[1];
sx q[1];
rz(-1.9903245) q[1];
sx q[1];
rz(-0.5330162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1905381) q[0];
sx q[0];
rz(-1.2742654) q[0];
sx q[0];
rz(-1.426479) q[0];
x q[1];
rz(2.7092878) q[2];
sx q[2];
rz(-3.084331) q[2];
sx q[2];
rz(-0.5106411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.24025) q[1];
sx q[1];
rz(-0.91687119) q[1];
sx q[1];
rz(2.4010977) q[1];
x q[2];
rz(0.80954062) q[3];
sx q[3];
rz(-1.8969826) q[3];
sx q[3];
rz(2.1128138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.816232) q[2];
sx q[2];
rz(-1.8133138) q[2];
sx q[2];
rz(-1.4556966) q[2];
rz(0.085518941) q[3];
sx q[3];
rz(-2.072008) q[3];
sx q[3];
rz(-2.1948309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7607255) q[0];
sx q[0];
rz(-0.65422288) q[0];
sx q[0];
rz(2.4192659) q[0];
rz(1.5707312) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(-3.0991203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0642424) q[0];
sx q[0];
rz(-1.1288861) q[0];
sx q[0];
rz(1.3674221) q[0];
x q[1];
rz(-3.1257079) q[2];
sx q[2];
rz(-1.4804891) q[2];
sx q[2];
rz(1.3549733) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3610246) q[1];
sx q[1];
rz(-2.3454822) q[1];
sx q[1];
rz(-1.3033504) q[1];
x q[2];
rz(2.9080584) q[3];
sx q[3];
rz(-0.81538768) q[3];
sx q[3];
rz(-2.5366207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8954358) q[2];
sx q[2];
rz(-1.7652721) q[2];
sx q[2];
rz(0.24406544) q[2];
rz(0.79143381) q[3];
sx q[3];
rz(-1.3118728) q[3];
sx q[3];
rz(0.70146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1112082) q[0];
sx q[0];
rz(-3.0272439) q[0];
sx q[0];
rz(2.9682888) q[0];
rz(-0.93470848) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(-0.25397837) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0365216) q[0];
sx q[0];
rz(-1.4582086) q[0];
sx q[0];
rz(-2.3238411) q[0];
rz(-pi) q[1];
rz(2.6298275) q[2];
sx q[2];
rz(-2.4825271) q[2];
sx q[2];
rz(-2.2816254) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18561126) q[1];
sx q[1];
rz(-0.69359457) q[1];
sx q[1];
rz(3.1090956) q[1];
x q[2];
rz(1.41296) q[3];
sx q[3];
rz(-0.39208347) q[3];
sx q[3];
rz(0.81721701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0678593) q[2];
sx q[2];
rz(-0.24418712) q[2];
sx q[2];
rz(1.7013288) q[2];
rz(-2.4654147) q[3];
sx q[3];
rz(-1.9192326) q[3];
sx q[3];
rz(3.0618099) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9192231) q[0];
sx q[0];
rz(-1.8310522) q[0];
sx q[0];
rz(-2.4134912) q[0];
rz(-1.1657731) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(0.56070915) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3325719) q[0];
sx q[0];
rz(-0.14223465) q[0];
sx q[0];
rz(-0.47551544) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54975551) q[2];
sx q[2];
rz(-1.4696331) q[2];
sx q[2];
rz(1.2899952) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47885671) q[1];
sx q[1];
rz(-2.6881892) q[1];
sx q[1];
rz(-1.0591566) q[1];
x q[2];
rz(-0.097499121) q[3];
sx q[3];
rz(-1.6336771) q[3];
sx q[3];
rz(-1.9815552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59549436) q[2];
sx q[2];
rz(-1.6807154) q[2];
sx q[2];
rz(1.072849) q[2];
rz(2.8192375) q[3];
sx q[3];
rz(-0.41283804) q[3];
sx q[3];
rz(-2.7505007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76674616) q[0];
sx q[0];
rz(-1.640919) q[0];
sx q[0];
rz(-2.7333976) q[0];
rz(-0.32632581) q[1];
sx q[1];
rz(-2.7944481) q[1];
sx q[1];
rz(1.4761285) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6107619) q[0];
sx q[0];
rz(-0.179804) q[0];
sx q[0];
rz(-1.7925949) q[0];
x q[1];
rz(-1.2017952) q[2];
sx q[2];
rz(-2.546406) q[2];
sx q[2];
rz(-0.068142224) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0261111) q[1];
sx q[1];
rz(-0.12899765) q[1];
sx q[1];
rz(1.1327101) q[1];
rz(-pi) q[2];
rz(2.9093698) q[3];
sx q[3];
rz(-1.7360064) q[3];
sx q[3];
rz(2.8286407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.25492111) q[2];
sx q[2];
rz(-1.6559867) q[2];
sx q[2];
rz(-0.86611789) q[2];
rz(2.4050889) q[3];
sx q[3];
rz(-2.0564506) q[3];
sx q[3];
rz(2.2542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0199468) q[0];
sx q[0];
rz(-3.096088) q[0];
sx q[0];
rz(1.2787) q[0];
rz(2.1389029) q[1];
sx q[1];
rz(-1.8808782) q[1];
sx q[1];
rz(-0.444828) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747975) q[0];
sx q[0];
rz(-0.54682362) q[0];
sx q[0];
rz(-1.2840136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1110531) q[2];
sx q[2];
rz(-0.90219775) q[2];
sx q[2];
rz(-0.3102613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5096566) q[1];
sx q[1];
rz(-0.77101123) q[1];
sx q[1];
rz(1.7421175) q[1];
x q[2];
rz(-0.37737198) q[3];
sx q[3];
rz(-1.8131953) q[3];
sx q[3];
rz(-0.19211543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43850809) q[2];
sx q[2];
rz(-1.727641) q[2];
sx q[2];
rz(0.50602305) q[2];
rz(-0.94432008) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(0.73616141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29189062) q[0];
sx q[0];
rz(-0.99221748) q[0];
sx q[0];
rz(-3.131409) q[0];
rz(-2.5306375) q[1];
sx q[1];
rz(-2.1612031) q[1];
sx q[1];
rz(0.082286509) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60264325) q[0];
sx q[0];
rz(-1.9398488) q[0];
sx q[0];
rz(-0.51765387) q[0];
rz(-pi) q[1];
rz(-0.23046979) q[2];
sx q[2];
rz(-0.90417143) q[2];
sx q[2];
rz(0.33483349) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4520784) q[1];
sx q[1];
rz(-1.5094712) q[1];
sx q[1];
rz(-0.21558) q[1];
rz(-1.786834) q[3];
sx q[3];
rz(-0.91751614) q[3];
sx q[3];
rz(2.9263484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33132195) q[2];
sx q[2];
rz(-1.4484582) q[2];
sx q[2];
rz(-2.4324379) q[2];
rz(1.6592615) q[3];
sx q[3];
rz(-0.63153657) q[3];
sx q[3];
rz(-2.6813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.822478) q[0];
sx q[0];
rz(-2.3111486) q[0];
sx q[0];
rz(-2.5122232) q[0];
rz(1.4156226) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(-1.1782882) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30692682) q[0];
sx q[0];
rz(-2.3884058) q[0];
sx q[0];
rz(-0.96145328) q[0];
x q[1];
rz(-2.7897116) q[2];
sx q[2];
rz(-2.0469672) q[2];
sx q[2];
rz(1.9989524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0337035) q[1];
sx q[1];
rz(-0.7830355) q[1];
sx q[1];
rz(1.4153925) q[1];
rz(2.4201336) q[3];
sx q[3];
rz(-2.0155523) q[3];
sx q[3];
rz(0.12192425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2668931) q[2];
sx q[2];
rz(-1.5723672) q[2];
sx q[2];
rz(-0.81864041) q[2];
rz(-1.5107752) q[3];
sx q[3];
rz(-1.2997593) q[3];
sx q[3];
rz(2.7394133) q[3];
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
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64869399) q[0];
sx q[0];
rz(-2.0317827) q[0];
sx q[0];
rz(-1.6775525) q[0];
rz(1.2661487) q[1];
sx q[1];
rz(-1.1504953) q[1];
sx q[1];
rz(-1.6082416) q[1];
rz(-0.89546236) q[2];
sx q[2];
rz(-1.3668899) q[2];
sx q[2];
rz(-1.5820506) q[2];
rz(-0.18465445) q[3];
sx q[3];
rz(-1.481537) q[3];
sx q[3];
rz(2.8164345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
