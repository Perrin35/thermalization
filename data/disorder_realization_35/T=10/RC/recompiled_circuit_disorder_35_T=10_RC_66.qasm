OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(-2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(2.6095698) q[1];
sx q[1];
rz(11.397059) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5380733) q[0];
sx q[0];
rz(-2.3581714) q[0];
sx q[0];
rz(-0.51622434) q[0];
x q[1];
rz(-0.92584445) q[2];
sx q[2];
rz(-2.1241509) q[2];
sx q[2];
rz(2.4990621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6572666) q[1];
sx q[1];
rz(-2.6063759) q[1];
sx q[1];
rz(2.3372997) q[1];
rz(-pi) q[2];
rz(0.3124247) q[3];
sx q[3];
rz(-1.6426622) q[3];
sx q[3];
rz(1.3952599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(1.8189836) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.009636119) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(1.0377201) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1514725) q[0];
sx q[0];
rz(-1.7904141) q[0];
sx q[0];
rz(-2.5734076) q[0];
x q[1];
rz(-1.0449045) q[2];
sx q[2];
rz(-1.0962588) q[2];
sx q[2];
rz(-1.431312) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2270826) q[1];
sx q[1];
rz(-2.7343379) q[1];
sx q[1];
rz(-1.6362908) q[1];
x q[2];
rz(-0.59229895) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(0.084687106) q[2];
rz(2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5304853) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(0.95570046) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37704913) q[0];
sx q[0];
rz(-2.0007613) q[0];
sx q[0];
rz(0.069805108) q[0];
rz(1.8884044) q[2];
sx q[2];
rz(-2.3177958) q[2];
sx q[2];
rz(-0.93712805) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9925476) q[1];
sx q[1];
rz(-2.6547975) q[1];
sx q[1];
rz(2.4942314) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7335448) q[3];
sx q[3];
rz(-1.1432262) q[3];
sx q[3];
rz(1.5023295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0009784) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-2.9476681) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(-2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.8594584) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(-0.55066806) q[0];
rz(-2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(0.36270025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6837316) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(1.213221) q[0];
rz(-pi) q[1];
rz(-1.0417468) q[2];
sx q[2];
rz(-3.1051271) q[2];
sx q[2];
rz(-2.1349825) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6994233) q[1];
sx q[1];
rz(-1.6740834) q[1];
sx q[1];
rz(-2.5073754) q[1];
rz(-pi) q[2];
rz(1.4705212) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68391934) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(3.1385699) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(1.4594706) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97809726) q[0];
sx q[0];
rz(-0.28143829) q[0];
sx q[0];
rz(1.0174169) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95894496) q[2];
sx q[2];
rz(-0.78275567) q[2];
sx q[2];
rz(-2.3252955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0046878) q[1];
sx q[1];
rz(-1.3740174) q[1];
sx q[1];
rz(-0.23006567) q[1];
rz(-pi) q[2];
rz(-0.32096433) q[3];
sx q[3];
rz(-1.1043613) q[3];
sx q[3];
rz(-2.5209559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(2.2996976) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.5766597) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(-0.58419624) q[0];
rz(-1.9110514) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(3.0029283) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4045227) q[0];
sx q[0];
rz(-1.5504019) q[0];
sx q[0];
rz(1.5741332) q[0];
rz(-pi) q[1];
rz(-2.9406592) q[2];
sx q[2];
rz(-1.4954508) q[2];
sx q[2];
rz(-2.2546525) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8530635) q[1];
sx q[1];
rz(-1.001734) q[1];
sx q[1];
rz(-2.2437375) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6824529) q[3];
sx q[3];
rz(-0.79677478) q[3];
sx q[3];
rz(-1.9285942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(1.8072051) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(-0.62430635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1034575) q[0];
sx q[0];
rz(-0.75224829) q[0];
sx q[0];
rz(2.8318846) q[0];
x q[1];
rz(-0.46632669) q[2];
sx q[2];
rz(-1.000324) q[2];
sx q[2];
rz(1.8898659) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99164167) q[1];
sx q[1];
rz(-2.2551564) q[1];
sx q[1];
rz(-0.73189862) q[1];
rz(-0.94675605) q[3];
sx q[3];
rz(-1.9692) q[3];
sx q[3];
rz(2.1053673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6162993) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-2.7590511) q[2];
rz(3.110102) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265825) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(3.0126742) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77057225) q[0];
sx q[0];
rz(-2.7311374) q[0];
sx q[0];
rz(1.8380941) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89959896) q[2];
sx q[2];
rz(-2.0313782) q[2];
sx q[2];
rz(-0.34924289) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0205295) q[1];
sx q[1];
rz(-1.9951092) q[1];
sx q[1];
rz(1.2709649) q[1];
rz(1.586732) q[3];
sx q[3];
rz(-0.83105479) q[3];
sx q[3];
rz(-2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(2.9028153) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(0.016013913) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(1.790766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5577561) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(1.7519959) q[0];
rz(-pi) q[1];
rz(-2.5596041) q[2];
sx q[2];
rz(-0.66170035) q[2];
sx q[2];
rz(2.2616507) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13547922) q[1];
sx q[1];
rz(-1.0891501) q[1];
sx q[1];
rz(0.7559795) q[1];
rz(-pi) q[2];
rz(2.8899367) q[3];
sx q[3];
rz(-2.3725384) q[3];
sx q[3];
rz(-2.3264399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0525557) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.982622) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(2.9283438) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(2.5949123) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0944259) q[0];
sx q[0];
rz(-1.7921653) q[0];
sx q[0];
rz(2.3012327) q[0];
rz(-0.47537739) q[2];
sx q[2];
rz(-1.7065085) q[2];
sx q[2];
rz(-0.00098468653) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.039779546) q[1];
sx q[1];
rz(-1.8598078) q[1];
sx q[1];
rz(1.6092369) q[1];
rz(-pi) q[2];
rz(2.91594) q[3];
sx q[3];
rz(-1.7060346) q[3];
sx q[3];
rz(1.7547363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13835779) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(-0.004301087) q[2];
rz(0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(2.6314541) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1417086) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.6179686) q[2];
sx q[2];
rz(-1.5856367) q[2];
sx q[2];
rz(-0.71935364) q[2];
rz(-1.0989582) q[3];
sx q[3];
rz(-1.8660587) q[3];
sx q[3];
rz(2.3793424) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];