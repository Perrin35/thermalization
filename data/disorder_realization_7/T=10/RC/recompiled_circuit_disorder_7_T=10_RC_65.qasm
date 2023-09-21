OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.1768567) q[0];
sx q[0];
rz(-0.86369792) q[0];
sx q[0];
rz(0.20110826) q[0];
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
rz(2.7456822) q[0];
sx q[0];
rz(-2.9183309) q[0];
sx q[0];
rz(-2.1366871) q[0];
rz(-1.8148412) q[2];
sx q[2];
rz(-1.913117) q[2];
sx q[2];
rz(0.92820456) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7658753) q[1];
sx q[1];
rz(-0.98343508) q[1];
sx q[1];
rz(2.767763) q[1];
rz(-pi) q[2];
rz(-0.17625531) q[3];
sx q[3];
rz(-1.3000543) q[3];
sx q[3];
rz(2.3034629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(-1.8480776) q[2];
rz(-2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023022368) q[0];
sx q[0];
rz(-1.1623628) q[0];
sx q[0];
rz(-0.9704216) q[0];
rz(-pi) q[1];
rz(1.8572349) q[2];
sx q[2];
rz(-1.0019433) q[2];
sx q[2];
rz(1.4494277) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7992236) q[1];
sx q[1];
rz(-2.2111726) q[1];
sx q[1];
rz(2.6745822) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3932088) q[3];
sx q[3];
rz(-0.41718306) q[3];
sx q[3];
rz(1.8821017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.3890283) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(3.1012428) q[0];
rz(-0.57998747) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(2.0592164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6425991) q[0];
sx q[0];
rz(-1.2046308) q[0];
sx q[0];
rz(0.55143349) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9897887) q[2];
sx q[2];
rz(-1.9215343) q[2];
sx q[2];
rz(2.0098067) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.97835097) q[1];
sx q[1];
rz(-0.32838531) q[1];
sx q[1];
rz(-2.3204625) q[1];
rz(-0.97614395) q[3];
sx q[3];
rz(-0.38446063) q[3];
sx q[3];
rz(2.773075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0064156) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(-0.55389261) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
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
rz(0.9961876) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-0.48809537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0220538) q[0];
sx q[0];
rz(-1.6484123) q[0];
sx q[0];
rz(-1.104943) q[0];
x q[1];
rz(-2.5627665) q[2];
sx q[2];
rz(-1.0597214) q[2];
sx q[2];
rz(1.6312815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0779873) q[1];
sx q[1];
rz(-1.4418774) q[1];
sx q[1];
rz(-1.7267666) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0933236) q[3];
sx q[3];
rz(-1.0415047) q[3];
sx q[3];
rz(2.570591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35486832) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(-1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(-0.3702634) q[0];
rz(-0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(-2.945074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58013433) q[0];
sx q[0];
rz(-0.41813865) q[0];
sx q[0];
rz(-0.74905507) q[0];
rz(2.0318803) q[2];
sx q[2];
rz(-1.5994674) q[2];
sx q[2];
rz(1.2210786) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6477485) q[1];
sx q[1];
rz(-1.3953679) q[1];
sx q[1];
rz(-2.8791048) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5444364) q[3];
sx q[3];
rz(-2.8215373) q[3];
sx q[3];
rz(-0.65359945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1281517) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(-2.2244942) q[2];
rz(2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-2.0687885) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(-1.3938168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48702792) q[0];
sx q[0];
rz(-1.9019039) q[0];
sx q[0];
rz(-1.4902671) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1962842) q[2];
sx q[2];
rz(-1.6746582) q[2];
sx q[2];
rz(1.9749157) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85642568) q[1];
sx q[1];
rz(-1.835583) q[1];
sx q[1];
rz(-1.1681639) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1145669) q[3];
sx q[3];
rz(-2.2531415) q[3];
sx q[3];
rz(1.2896468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(2.2992772) q[2];
rz(2.0541644) q[3];
sx q[3];
rz(-0.73173404) q[3];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6120646) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(1.0213617) q[0];
rz(-1.1313324) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(2.9313415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410941) q[0];
sx q[0];
rz(-0.14963089) q[0];
sx q[0];
rz(1.11087) q[0];
rz(-pi) q[1];
rz(1.636593) q[2];
sx q[2];
rz(-1.734126) q[2];
sx q[2];
rz(0.79236275) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6639858) q[1];
sx q[1];
rz(-0.37316445) q[1];
sx q[1];
rz(-0.54877703) q[1];
rz(-pi) q[2];
rz(0.93123575) q[3];
sx q[3];
rz(-2.3124472) q[3];
sx q[3];
rz(-2.2614546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(-2.6953221) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(-0.97672021) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57709549) q[0];
sx q[0];
rz(-1.4221039) q[0];
sx q[0];
rz(-0.45678267) q[0];
rz(-pi) q[1];
rz(0.080539695) q[2];
sx q[2];
rz(-1.7320247) q[2];
sx q[2];
rz(-1.5921519) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.50423065) q[1];
sx q[1];
rz(-2.0081875) q[1];
sx q[1];
rz(-0.38423844) q[1];
rz(0.85985698) q[3];
sx q[3];
rz(-1.1338187) q[3];
sx q[3];
rz(2.6881998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7405159) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(2.5891499) q[2];
rz(2.6756514) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(2.1140816) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(-2.2209432) q[0];
rz(-3.0905511) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(2.2424973) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5024922) q[0];
sx q[0];
rz(-2.149625) q[0];
sx q[0];
rz(-2.6274908) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1295206) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(-0.65288359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2194849) q[1];
sx q[1];
rz(-1.9003632) q[1];
sx q[1];
rz(3.0199514) q[1];
x q[2];
rz(1.5713463) q[3];
sx q[3];
rz(-0.37017248) q[3];
sx q[3];
rz(2.7412424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1961394) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(-2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33912441) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(-0.39696473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1966232) q[0];
sx q[0];
rz(-2.0485749) q[0];
sx q[0];
rz(-1.0265795) q[0];
rz(-pi) q[1];
rz(1.9372349) q[2];
sx q[2];
rz(-1.2841184) q[2];
sx q[2];
rz(-2.3057126) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1399999) q[1];
sx q[1];
rz(-1.4824176) q[1];
sx q[1];
rz(-1.2338851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0730562) q[3];
sx q[3];
rz(-1.2101678) q[3];
sx q[3];
rz(1.8325168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94840702) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5554572) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(3.0336822) q[2];
sx q[2];
rz(-1.6549329) q[2];
sx q[2];
rz(-1.8555117) q[2];
rz(0.11937033) q[3];
sx q[3];
rz(-1.3138249) q[3];
sx q[3];
rz(-1.8520595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
