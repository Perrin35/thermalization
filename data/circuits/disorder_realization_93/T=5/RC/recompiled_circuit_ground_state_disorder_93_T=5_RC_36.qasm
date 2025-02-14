OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8146347) q[0];
sx q[0];
rz(3.1231413) q[0];
sx q[0];
rz(9.2508247) q[0];
rz(-1.3375018) q[1];
sx q[1];
rz(-0.87556535) q[1];
sx q[1];
rz(-2.8101885) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6397996) q[0];
sx q[0];
rz(-1.26609) q[0];
sx q[0];
rz(-1.7250389) q[0];
x q[1];
rz(0.5919257) q[2];
sx q[2];
rz(-2.3560696) q[2];
sx q[2];
rz(2.3374918) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1494721) q[1];
sx q[1];
rz(-1.3821342) q[1];
sx q[1];
rz(-2.4558398) q[1];
rz(-pi) q[2];
rz(-1.5378733) q[3];
sx q[3];
rz(-1.415807) q[3];
sx q[3];
rz(-1.8553569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6203674) q[2];
sx q[2];
rz(-1.2646893) q[2];
sx q[2];
rz(-2.550726) q[2];
rz(-1.8477731) q[3];
sx q[3];
rz(-1.7655244) q[3];
sx q[3];
rz(0.53904343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5188326) q[0];
sx q[0];
rz(-2.9157214) q[0];
sx q[0];
rz(2.2665005) q[0];
rz(3.0488293) q[1];
sx q[1];
rz(-2.0848672) q[1];
sx q[1];
rz(-2.1327532) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1160285) q[0];
sx q[0];
rz(-2.7112428) q[0];
sx q[0];
rz(-1.5034281) q[0];
x q[1];
rz(-1.5447727) q[2];
sx q[2];
rz(-0.98434005) q[2];
sx q[2];
rz(-1.9860739) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5238467) q[1];
sx q[1];
rz(-2.2154795) q[1];
sx q[1];
rz(2.9545386) q[1];
x q[2];
rz(1.6289234) q[3];
sx q[3];
rz(-1.584313) q[3];
sx q[3];
rz(1.112354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2729317) q[2];
sx q[2];
rz(-1.2596143) q[2];
sx q[2];
rz(1.7449562) q[2];
rz(0.40846387) q[3];
sx q[3];
rz(-2.4558081) q[3];
sx q[3];
rz(2.7062611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246966) q[0];
sx q[0];
rz(-2.8530687) q[0];
sx q[0];
rz(-2.1911279) q[0];
rz(1.2735927) q[1];
sx q[1];
rz(-1.344695) q[1];
sx q[1];
rz(1.3709995) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4038611) q[0];
sx q[0];
rz(-1.7086242) q[0];
sx q[0];
rz(-0.65562427) q[0];
x q[1];
rz(-2.8954178) q[2];
sx q[2];
rz(-0.44533768) q[2];
sx q[2];
rz(-0.37644437) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6505988) q[1];
sx q[1];
rz(-1.4447325) q[1];
sx q[1];
rz(-1.0647414) q[1];
rz(-2.9955735) q[3];
sx q[3];
rz(-0.64490025) q[3];
sx q[3];
rz(0.19949808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3889968) q[2];
sx q[2];
rz(-1.4138556) q[2];
sx q[2];
rz(-1.5163806) q[2];
rz(-0.73836941) q[3];
sx q[3];
rz(-1.5117896) q[3];
sx q[3];
rz(-2.4228158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7826638) q[0];
sx q[0];
rz(-1.5072701) q[0];
sx q[0];
rz(1.9184817) q[0];
rz(-0.75543985) q[1];
sx q[1];
rz(-1.1399882) q[1];
sx q[1];
rz(3.0208407) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31407315) q[0];
sx q[0];
rz(-1.4085007) q[0];
sx q[0];
rz(1.4050964) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5865929) q[2];
sx q[2];
rz(-1.3669802) q[2];
sx q[2];
rz(2.7380916) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1270406) q[1];
sx q[1];
rz(-2.2071731) q[1];
sx q[1];
rz(-2.4680952) q[1];
rz(-pi) q[2];
rz(-0.018649696) q[3];
sx q[3];
rz(-1.5232764) q[3];
sx q[3];
rz(-1.7619676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52189031) q[2];
sx q[2];
rz(-1.4554687) q[2];
sx q[2];
rz(-1.5187368) q[2];
rz(-1.8146727) q[3];
sx q[3];
rz(-0.36553317) q[3];
sx q[3];
rz(-2.0324619) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0066852) q[0];
sx q[0];
rz(-1.4445855) q[0];
sx q[0];
rz(2.4476442) q[0];
rz(-2.7130983) q[1];
sx q[1];
rz(-2.4887812) q[1];
sx q[1];
rz(1.8983715) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080414198) q[0];
sx q[0];
rz(-1.616487) q[0];
sx q[0];
rz(-1.3752328) q[0];
x q[1];
rz(2.9712293) q[2];
sx q[2];
rz(-1.3479832) q[2];
sx q[2];
rz(0.52056584) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4116333) q[1];
sx q[1];
rz(-2.6891629) q[1];
sx q[1];
rz(-1.0414473) q[1];
x q[2];
rz(-2.9443191) q[3];
sx q[3];
rz(-2.6515238) q[3];
sx q[3];
rz(-2.3207612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1474233) q[2];
sx q[2];
rz(-1.8403515) q[2];
sx q[2];
rz(-2.6541138) q[2];
rz(3.0726037) q[3];
sx q[3];
rz(-2.4853849) q[3];
sx q[3];
rz(2.5330353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92031205) q[0];
sx q[0];
rz(-0.99928904) q[0];
sx q[0];
rz(-0.47166237) q[0];
rz(-0.77700514) q[1];
sx q[1];
rz(-2.0174556) q[1];
sx q[1];
rz(1.5832925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7651833) q[0];
sx q[0];
rz(-1.5411955) q[0];
sx q[0];
rz(-1.5628003) q[0];
rz(-pi) q[1];
rz(0.42209322) q[2];
sx q[2];
rz(-1.0703841) q[2];
sx q[2];
rz(1.3304324) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.029249749) q[1];
sx q[1];
rz(-1.3549374) q[1];
sx q[1];
rz(-0.4451181) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8804017) q[3];
sx q[3];
rz(-2.3979366) q[3];
sx q[3];
rz(-0.82871515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9787057) q[2];
sx q[2];
rz(-1.6521613) q[2];
sx q[2];
rz(1.1631274) q[2];
rz(-0.85824054) q[3];
sx q[3];
rz(-1.6909928) q[3];
sx q[3];
rz(0.8391909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6100886) q[0];
sx q[0];
rz(-1.4634345) q[0];
sx q[0];
rz(2.8186744) q[0];
rz(-1.8220176) q[1];
sx q[1];
rz(-2.0009305) q[1];
sx q[1];
rz(-1.742977) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6933029) q[0];
sx q[0];
rz(-2.9041661) q[0];
sx q[0];
rz(-0.7619995) q[0];
rz(-pi) q[1];
rz(-2.7432199) q[2];
sx q[2];
rz(-2.1443416) q[2];
sx q[2];
rz(-0.15352098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7306644) q[1];
sx q[1];
rz(-2.321744) q[1];
sx q[1];
rz(1.6150722) q[1];
x q[2];
rz(-3.0121376) q[3];
sx q[3];
rz(-1.4788179) q[3];
sx q[3];
rz(-1.5935073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4601595) q[2];
sx q[2];
rz(-0.56542772) q[2];
sx q[2];
rz(-2.4088755) q[2];
rz(-0.021746246) q[3];
sx q[3];
rz(-1.6201868) q[3];
sx q[3];
rz(-2.963613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9916423) q[0];
sx q[0];
rz(-0.1217148) q[0];
sx q[0];
rz(0.44418401) q[0];
rz(1.3581879) q[1];
sx q[1];
rz(-1.5779747) q[1];
sx q[1];
rz(2.0213199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2884707) q[0];
sx q[0];
rz(-0.74672304) q[0];
sx q[0];
rz(1.4190559) q[0];
rz(0.48768576) q[2];
sx q[2];
rz(-1.926406) q[2];
sx q[2];
rz(-1.5118701) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6669338) q[1];
sx q[1];
rz(-0.28945112) q[1];
sx q[1];
rz(0.54929026) q[1];
rz(-0.46457477) q[3];
sx q[3];
rz(-1.9442055) q[3];
sx q[3];
rz(-0.68313235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84919763) q[2];
sx q[2];
rz(-2.2870543) q[2];
sx q[2];
rz(-1.9067859) q[2];
rz(-0.97709996) q[3];
sx q[3];
rz(-1.2950725) q[3];
sx q[3];
rz(2.5375514) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4567756) q[0];
sx q[0];
rz(-0.026739459) q[0];
sx q[0];
rz(-2.7944414) q[0];
rz(-0.30255643) q[1];
sx q[1];
rz(-1.1056933) q[1];
sx q[1];
rz(0.0060630719) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13841471) q[0];
sx q[0];
rz(-0.18269953) q[0];
sx q[0];
rz(-1.4783786) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89170189) q[2];
sx q[2];
rz(-1.3366586) q[2];
sx q[2];
rz(-2.8193223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8703985) q[1];
sx q[1];
rz(-2.4358303) q[1];
sx q[1];
rz(2.0466385) q[1];
rz(-pi) q[2];
rz(1.6177032) q[3];
sx q[3];
rz(-0.79104086) q[3];
sx q[3];
rz(2.6447846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89846316) q[2];
sx q[2];
rz(-2.1738238) q[2];
sx q[2];
rz(1.3904307) q[2];
rz(0.27680382) q[3];
sx q[3];
rz(-1.0172458) q[3];
sx q[3];
rz(2.0738585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7368363) q[0];
sx q[0];
rz(-2.3911349) q[0];
sx q[0];
rz(0.87338895) q[0];
rz(-0.5492754) q[1];
sx q[1];
rz(-0.6540238) q[1];
sx q[1];
rz(-2.3347847) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72708541) q[0];
sx q[0];
rz(-2.1860518) q[0];
sx q[0];
rz(1.8919957) q[0];
rz(-pi) q[1];
rz(-2.4960356) q[2];
sx q[2];
rz(-0.19536138) q[2];
sx q[2];
rz(-0.33249172) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0765105) q[1];
sx q[1];
rz(-1.8190093) q[1];
sx q[1];
rz(3.0221171) q[1];
rz(-0.30183123) q[3];
sx q[3];
rz(-0.98205295) q[3];
sx q[3];
rz(1.4422688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.4775382) q[2];
sx q[2];
rz(-1.8287903) q[2];
sx q[2];
rz(0.99267268) q[2];
rz(-2.4382675) q[3];
sx q[3];
rz(-1.1329009) q[3];
sx q[3];
rz(0.68737427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(1.0032046) q[0];
sx q[0];
rz(-0.78835543) q[0];
sx q[0];
rz(2.0302571) q[0];
rz(-2.5913024) q[1];
sx q[1];
rz(-2.7688409) q[1];
sx q[1];
rz(0.54584835) q[1];
rz(0.69725488) q[2];
sx q[2];
rz(-1.6959126) q[2];
sx q[2];
rz(-0.87903862) q[2];
rz(2.417682) q[3];
sx q[3];
rz(-0.58876287) q[3];
sx q[3];
rz(2.7947254) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
