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
rz(4.9464524) q[1];
sx q[1];
rz(10.051605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39591046) q[0];
sx q[0];
rz(-0.22326176) q[0];
sx q[0];
rz(-1.0049055) q[0];
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
rz(pi/2) q[0];
rz(-2.1605692) q[1];
sx q[1];
rz(-1.8796646) q[1];
sx q[1];
rz(-0.94998756) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8455891) q[3];
sx q[3];
rz(-1.4010251) q[3];
sx q[3];
rz(-0.78026375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(-1.8480776) q[2];
rz(-1.1039929) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(2.1564116) q[0];
rz(-2.7424116) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(0.10736297) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0221314) q[0];
sx q[0];
rz(-2.4298926) q[0];
sx q[0];
rz(2.2244781) q[0];
rz(0.58789247) q[2];
sx q[2];
rz(-1.8111472) q[2];
sx q[2];
rz(-0.27871486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7992236) q[1];
sx q[1];
rz(-2.2111726) q[1];
sx q[1];
rz(0.46701042) q[1];
rz(-pi) q[2];
x q[2];
rz(1.739005) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5169516) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-3.1012428) q[0];
rz(0.57998747) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(-2.0592164) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6862415) q[0];
sx q[0];
rz(-0.65128122) q[0];
sx q[0];
rz(-2.5097646) q[0];
x q[1];
rz(1.9627377) q[2];
sx q[2];
rz(-2.760663) q[2];
sx q[2];
rz(-2.4286963) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(2.9134681) q[1];
rz(-2.1654487) q[3];
sx q[3];
rz(-2.757132) q[3];
sx q[3];
rz(-0.36851766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0064156) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(-2.5877) q[2];
rz(-1.6484377) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(0.035382263) q[0];
rz(-0.9961876) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(-0.48809537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8433544) q[0];
sx q[0];
rz(-0.47180628) q[0];
sx q[0];
rz(-1.3993553) q[0];
rz(-pi) q[1];
rz(2.1610291) q[2];
sx q[2];
rz(-2.0681941) q[2];
sx q[2];
rz(-0.24888466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6285703) q[1];
sx q[1];
rz(-1.7254618) q[1];
sx q[1];
rz(0.13048529) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.048269) q[3];
sx q[3];
rz(-1.0415047) q[3];
sx q[3];
rz(-0.57100163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7867243) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(-2.6341237) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5421211) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(-2.7713293) q[0];
rz(0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(2.945074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5614583) q[0];
sx q[0];
rz(-0.41813865) q[0];
sx q[0];
rz(-0.74905507) q[0];
rz(-0.032012149) q[2];
sx q[2];
rz(-2.0316761) q[2];
sx q[2];
rz(-0.33547685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6477485) q[1];
sx q[1];
rz(-1.7462247) q[1];
sx q[1];
rz(-2.8791048) q[1];
rz(1.4007832) q[3];
sx q[3];
rz(-1.2983054) q[3];
sx q[3];
rz(-1.221399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(-0.91709843) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(-2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31345263) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(2.939558) q[0];
rz(1.0728041) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(1.7477759) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0575384) q[0];
sx q[0];
rz(-1.6469427) q[0];
sx q[0];
rz(-0.3321068) q[0];
rz(-pi) q[1];
rz(1.9453085) q[2];
sx q[2];
rz(-1.4669344) q[2];
sx q[2];
rz(1.166677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60336941) q[1];
sx q[1];
rz(-1.1829611) q[1];
sx q[1];
rz(-2.8549854) q[1];
rz(-1.1145669) q[3];
sx q[3];
rz(-0.8884512) q[3];
sx q[3];
rz(-1.8519459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(2.2992772) q[2];
rz(-1.0874282) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(-1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6120646) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(-1.1313324) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(-0.21025118) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267669) q[0];
sx q[0];
rz(-1.504577) q[0];
sx q[0];
rz(-1.4365175) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5049997) q[2];
sx q[2];
rz(-1.4074667) q[2];
sx q[2];
rz(-0.79236275) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6639858) q[1];
sx q[1];
rz(-2.7684282) q[1];
sx q[1];
rz(0.54877703) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.290091) q[3];
sx q[3];
rz(-1.1151033) q[3];
sx q[3];
rz(-1.9853026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82688275) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(1.7604527) q[2];
rz(-2.3794877) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(2.7281318) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.64637) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(-0.58724171) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(2.1648724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2205692) q[0];
sx q[0];
rz(-2.0221634) q[0];
sx q[0];
rz(1.4054106) q[0];
x q[1];
rz(-1.4090528) q[2];
sx q[2];
rz(-1.4913034) q[2];
sx q[2];
rz(-0.0083991945) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8752012) q[1];
sx q[1];
rz(-2.5677263) q[1];
sx q[1];
rz(-2.2465474) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2817357) q[3];
sx q[3];
rz(-2.007774) q[3];
sx q[3];
rz(2.6881998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7405159) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-0.66822806) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(0.051041516) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(2.2424973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5024922) q[0];
sx q[0];
rz(-2.149625) q[0];
sx q[0];
rz(-0.51410189) q[0];
rz(-0.33234889) q[2];
sx q[2];
rz(-1.0370478) q[2];
sx q[2];
rz(0.74408434) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60914674) q[1];
sx q[1];
rz(-1.4557314) q[1];
sx q[1];
rz(-1.2389517) q[1];
x q[2];
rz(-0.00021342834) q[3];
sx q[3];
rz(-1.2006239) q[3];
sx q[3];
rz(2.7406524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(-0.10458874) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(-2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(2.7446279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1966232) q[0];
sx q[0];
rz(-1.0930177) q[0];
sx q[0];
rz(-2.1150132) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9372349) q[2];
sx q[2];
rz(-1.2841184) q[2];
sx q[2];
rz(0.8358801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18394477) q[1];
sx q[1];
rz(-0.34788222) q[1];
sx q[1];
rz(1.3089048) q[1];
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
rz(0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(-1.7276673) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(-2.0146501) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(-3.0336822) q[2];
sx q[2];
rz(-1.4866598) q[2];
sx q[2];
rz(1.286081) q[2];
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
