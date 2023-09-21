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
rz(5.4194874) q[0];
sx q[0];
rz(9.6258862) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(4.9464524) q[1];
sx q[1];
rz(10.051605) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5212844) q[0];
sx q[0];
rz(-1.4518019) q[0];
sx q[0];
rz(1.3814397) q[0];
x q[1];
rz(0.5958545) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(-2.8516304) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7658753) q[1];
sx q[1];
rz(-2.1581576) q[1];
sx q[1];
rz(0.37382965) q[1];
rz(1.8455891) q[3];
sx q[3];
rz(-1.4010251) q[3];
sx q[3];
rz(0.78026375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(-1.8480776) q[2];
rz(1.1039929) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-2.3867992) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-1.2626516) q[1];
sx q[1];
rz(0.10736297) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1194612) q[0];
sx q[0];
rz(-2.4298926) q[0];
sx q[0];
rz(-0.91711451) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41610833) q[2];
sx q[2];
rz(-0.62972087) q[2];
sx q[2];
rz(0.94905084) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.098298479) q[1];
sx q[1];
rz(-0.77273332) q[1];
sx q[1];
rz(1.0272825) q[1];
rz(1.739005) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5169516) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(0.055756904) q[3];
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
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(-1.0823762) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530281) q[0];
sx q[0];
rz(-2.0819426) q[0];
sx q[0];
rz(-1.1477863) q[0];
rz(-2.9897887) q[2];
sx q[2];
rz(-1.2200583) q[2];
sx q[2];
rz(-2.0098067) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.12990002) q[1];
sx q[1];
rz(-1.7923647) q[1];
sx q[1];
rz(1.8151912) q[1];
x q[2];
rz(1.8941746) q[3];
sx q[3];
rz(-1.3591027) q[3];
sx q[3];
rz(-2.4993103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(0.55389261) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(-0.55251399) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4937113) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(0.035382263) q[0];
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
rz(1.1195388) q[0];
sx q[0];
rz(-1.4931803) q[0];
sx q[0];
rz(2.0366497) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5627665) q[2];
sx q[2];
rz(-1.0597214) q[2];
sx q[2];
rz(-1.5103112) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9490337) q[1];
sx q[1];
rz(-2.9395736) q[1];
sx q[1];
rz(-0.87534027) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.048269) q[3];
sx q[3];
rz(-1.0415047) q[3];
sx q[3];
rz(2.570591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35486832) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(-1.1821702) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(-1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.5421211) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(-2.7713293) q[0];
rz(-2.8778991) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(-2.945074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5614583) q[0];
sx q[0];
rz(-2.723454) q[0];
sx q[0];
rz(2.3925376) q[0];
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
rz(pi/2) q[0];
sx q[0];
rz(-0.49911753) q[1];
sx q[1];
rz(-0.31458464) q[1];
sx q[1];
rz(2.5423074) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7408095) q[3];
sx q[3];
rz(-1.8432872) q[3];
sx q[3];
rz(-1.9201936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1281517) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(-0.91709843) q[2];
rz(0.83459485) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.6348811) q[1];
sx q[1];
rz(1.3938168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73035875) q[0];
sx q[0];
rz(-2.8011836) q[0];
sx q[0];
rz(0.22986869) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0300573) q[2];
sx q[2];
rz(-1.1984014) q[2];
sx q[2];
rz(-2.6967449) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85642568) q[1];
sx q[1];
rz(-1.835583) q[1];
sx q[1];
rz(1.1681639) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73563852) q[3];
sx q[3];
rz(-1.2218352) q[3];
sx q[3];
rz(2.5603106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1919057) q[1];
sx q[1];
rz(0.21025118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70049858) q[0];
sx q[0];
rz(-0.14963089) q[0];
sx q[0];
rz(-2.0307226) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37961752) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(-1.1773674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0381895) q[1];
sx q[1];
rz(-1.2545171) q[1];
sx q[1];
rz(-1.7722539) q[1];
rz(0.57742124) q[3];
sx q[3];
rz(-0.93772674) q[3];
sx q[3];
rz(3.0949743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82688275) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(-1.38114) q[2];
rz(-2.3794877) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(-2.7281318) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(-0.58724171) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(0.97672021) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2864838) q[0];
sx q[0];
rz(-0.47874641) q[0];
sx q[0];
rz(0.32740645) q[0];
rz(-pi) q[1];
rz(-1.4090528) q[2];
sx q[2];
rz(-1.6502893) q[2];
sx q[2];
rz(-3.1331935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8752012) q[1];
sx q[1];
rz(-0.57386639) q[1];
sx q[1];
rz(0.89504524) q[1];
rz(2.5891853) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(1.4668902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(-1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(2.2424973) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9793648) q[0];
sx q[0];
rz(-0.75408903) q[0];
sx q[0];
rz(-0.92569949) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0663509) q[2];
sx q[2];
rz(-0.62014183) q[2];
sx q[2];
rz(-1.3401741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5324459) q[1];
sx q[1];
rz(-1.4557314) q[1];
sx q[1];
rz(-1.2389517) q[1];
rz(-1.9409688) q[3];
sx q[3];
rz(-1.5709953) q[3];
sx q[3];
rz(1.1699333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(-3.0370039) q[2];
rz(-0.83834046) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33912441) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(2.4291908) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(-2.7446279) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2755651) q[0];
sx q[0];
rz(-2.4336928) q[0];
sx q[0];
rz(-2.3562354) q[0];
x q[1];
rz(2.2592779) q[2];
sx q[2];
rz(-0.46122641) q[2];
sx q[2];
rz(-0.099909401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0015928) q[1];
sx q[1];
rz(-1.4824176) q[1];
sx q[1];
rz(1.9077076) q[1];
x q[2];
rz(2.0685365) q[3];
sx q[3];
rz(-1.9314249) q[3];
sx q[3];
rz(-1.3090759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(3.0449384) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.5861355) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(1.8718406) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(-0.10791049) q[2];
sx q[2];
rz(-1.6549329) q[2];
sx q[2];
rz(-1.8555117) q[2];
rz(1.1453015) q[3];
sx q[3];
rz(-0.28278657) q[3];
sx q[3];
rz(0.84859802) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
