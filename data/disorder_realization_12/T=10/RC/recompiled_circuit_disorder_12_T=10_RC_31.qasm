OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.84500399) q[0];
sx q[0];
rz(-0.72405976) q[0];
sx q[0];
rz(1.4847423) q[0];
rz(1.2031263) q[1];
sx q[1];
rz(-0.523518) q[1];
sx q[1];
rz(2.2533921) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523525) q[0];
sx q[0];
rz(-1.7416735) q[0];
sx q[0];
rz(1.5255552) q[0];
rz(-pi) q[1];
rz(2.1076803) q[2];
sx q[2];
rz(-1.7134588) q[2];
sx q[2];
rz(-0.85096525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8623212) q[1];
sx q[1];
rz(-0.16695484) q[1];
sx q[1];
rz(0.0032940666) q[1];
rz(1.5376484) q[3];
sx q[3];
rz(-2.2691233) q[3];
sx q[3];
rz(0.87227277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28329864) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(-1.1179914) q[2];
rz(0.14532146) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(3.0956691) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730826) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(-0.65482393) q[0];
rz(1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-0.12589802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4100285) q[0];
sx q[0];
rz(-1.0985939) q[0];
sx q[0];
rz(-0.41759755) q[0];
x q[1];
rz(-0.84463859) q[2];
sx q[2];
rz(-2.3181049) q[2];
sx q[2];
rz(-0.57074947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8996846) q[1];
sx q[1];
rz(-0.13026127) q[1];
sx q[1];
rz(-1.3379407) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1308261) q[3];
sx q[3];
rz(-0.84732238) q[3];
sx q[3];
rz(2.7303498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(2.753567) q[2];
rz(-1.4240501) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.5858784) q[0];
sx q[0];
rz(-0.782574) q[0];
sx q[0];
rz(-0.079332381) q[0];
rz(-3.0575867) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(1.1598587) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0214329) q[0];
sx q[0];
rz(-0.5172356) q[0];
sx q[0];
rz(-1.4034127) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1986387) q[2];
sx q[2];
rz(-1.9938333) q[2];
sx q[2];
rz(1.1689651) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68418903) q[1];
sx q[1];
rz(-2.0261129) q[1];
sx q[1];
rz(1.1475569) q[1];
rz(0.056315259) q[3];
sx q[3];
rz(-1.0259797) q[3];
sx q[3];
rz(-0.69308263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40538654) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(2.4978499) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.165034) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.4820341) q[0];
rz(2.4404793) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(2.8569417) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4581504) q[0];
sx q[0];
rz(-1.0939286) q[0];
sx q[0];
rz(-0.083806888) q[0];
x q[1];
rz(0.20547159) q[2];
sx q[2];
rz(-1.2887508) q[2];
sx q[2];
rz(-3.0021283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.804867) q[1];
sx q[1];
rz(-1.6972099) q[1];
sx q[1];
rz(-3.0625507) q[1];
rz(-pi) q[2];
rz(0.90162006) q[3];
sx q[3];
rz(-1.1812783) q[3];
sx q[3];
rz(-0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9099137) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(2.5740734) q[2];
rz(0.41401687) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48859566) q[0];
sx q[0];
rz(-0.96019205) q[0];
sx q[0];
rz(-1.3265142) q[0];
rz(1.1524221) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(2.2036536) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4612761) q[0];
sx q[0];
rz(-0.094705908) q[0];
sx q[0];
rz(-1.8041496) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1300975) q[2];
sx q[2];
rz(-2.0213631) q[2];
sx q[2];
rz(0.089103854) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.282498) q[1];
sx q[1];
rz(-1.5364093) q[1];
sx q[1];
rz(3.1153468) q[1];
rz(0.35477562) q[3];
sx q[3];
rz(-1.4491023) q[3];
sx q[3];
rz(2.0718758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9662629) q[2];
sx q[2];
rz(-2.9523409) q[2];
sx q[2];
rz(1.0413292) q[2];
rz(-1.3025618) q[3];
sx q[3];
rz(-1.1338736) q[3];
sx q[3];
rz(1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(-2.5842216) q[0];
rz(0.56466651) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(-2.5851137) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59069955) q[0];
sx q[0];
rz(-1.3644344) q[0];
sx q[0];
rz(2.7569689) q[0];
rz(1.3715903) q[2];
sx q[2];
rz(-2.4161985) q[2];
sx q[2];
rz(0.2142011) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0007799) q[1];
sx q[1];
rz(-1.8218578) q[1];
sx q[1];
rz(-0.18745969) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3658386) q[3];
sx q[3];
rz(-1.6072825) q[3];
sx q[3];
rz(-3.0690103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53211987) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(1.9167985) q[2];
rz(-1.5103643) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0790134) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(-0.042908948) q[0];
rz(-0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(2.6409805) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6734877) q[0];
sx q[0];
rz(-1.3930495) q[0];
sx q[0];
rz(-0.034981473) q[0];
x q[1];
rz(2.1130354) q[2];
sx q[2];
rz(-1.8430029) q[2];
sx q[2];
rz(-2.3173463) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7350406) q[1];
sx q[1];
rz(-2.7222689) q[1];
sx q[1];
rz(-2.8419644) q[1];
rz(-pi) q[2];
rz(-0.41916267) q[3];
sx q[3];
rz(-0.8014285) q[3];
sx q[3];
rz(-2.6475346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6033972) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(-0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(-1.0674397) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(-1.6759466) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38795234) q[0];
sx q[0];
rz(-2.6590829) q[0];
sx q[0];
rz(1.2778736) q[0];
rz(-pi) q[1];
rz(-2.6162716) q[2];
sx q[2];
rz(-2.9267731) q[2];
sx q[2];
rz(0.36052442) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9626179) q[1];
sx q[1];
rz(-1.0777506) q[1];
sx q[1];
rz(-2.7589873) q[1];
x q[2];
rz(2.500324) q[3];
sx q[3];
rz(-1.0675758) q[3];
sx q[3];
rz(-0.36676952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8026768) q[2];
sx q[2];
rz(-2.2854476) q[2];
sx q[2];
rz(-1.476293) q[2];
rz(0.87351292) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840238) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(-2.4043758) q[0];
rz(3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(0.92528701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7172456) q[0];
sx q[0];
rz(-2.2204917) q[0];
sx q[0];
rz(2.2109277) q[0];
rz(-pi) q[1];
rz(-1.4383573) q[2];
sx q[2];
rz(-2.2297511) q[2];
sx q[2];
rz(3.1140285) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8268938) q[1];
sx q[1];
rz(-1.3876545) q[1];
sx q[1];
rz(-1.7032743) q[1];
rz(1.0989283) q[3];
sx q[3];
rz(-0.77583757) q[3];
sx q[3];
rz(-2.8692506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3725738) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(-2.1378689) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-3.1019679) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(-1.8819303) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(-2.3840747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.789061) q[0];
sx q[0];
rz(-1.7081982) q[0];
sx q[0];
rz(1.765541) q[0];
rz(-0.9805571) q[2];
sx q[2];
rz(-1.440181) q[2];
sx q[2];
rz(-0.57934258) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1265035) q[1];
sx q[1];
rz(-1.1474097) q[1];
sx q[1];
rz(0.24366118) q[1];
rz(0.032978756) q[3];
sx q[3];
rz(-0.71027256) q[3];
sx q[3];
rz(-0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99047986) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(-2.5027067) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(1.0206153) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6486075) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(-1.6408625) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(2.5095148) q[2];
sx q[2];
rz(-2.6161604) q[2];
sx q[2];
rz(-2.5642774) q[2];
rz(-2.1502675) q[3];
sx q[3];
rz(-1.2990868) q[3];
sx q[3];
rz(2.2494863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];