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
rz(1.5384742) q[0];
sx q[0];
rz(4.8967946) q[0];
sx q[0];
rz(8.1660981) q[0];
rz(-1.7006682) q[1];
sx q[1];
rz(-0.52803334) q[1];
sx q[1];
rz(-0.33906403) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067056) q[0];
sx q[0];
rz(-2.7769222) q[0];
sx q[0];
rz(-1.986978) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6554675) q[2];
sx q[2];
rz(-2.4007501) q[2];
sx q[2];
rz(1.6227826) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69341371) q[1];
sx q[1];
rz(-0.95985818) q[1];
sx q[1];
rz(0.081206755) q[1];
rz(0.65646274) q[3];
sx q[3];
rz(-1.5243666) q[3];
sx q[3];
rz(-1.9551147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49860927) q[2];
sx q[2];
rz(-0.53350854) q[2];
sx q[2];
rz(1.741893) q[2];
rz(1.9881605) q[3];
sx q[3];
rz(-1.5164991) q[3];
sx q[3];
rz(1.7216518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4145849) q[0];
sx q[0];
rz(-1.4684533) q[0];
sx q[0];
rz(-0.46838316) q[0];
rz(0.98764694) q[1];
sx q[1];
rz(-1.3750319) q[1];
sx q[1];
rz(1.04029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0752234) q[0];
sx q[0];
rz(-0.98261896) q[0];
sx q[0];
rz(-2.4494314) q[0];
rz(-2.3556893) q[2];
sx q[2];
rz(-0.52009798) q[2];
sx q[2];
rz(1.4532068) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6436983) q[1];
sx q[1];
rz(-2.0308745) q[1];
sx q[1];
rz(-1.6804477) q[1];
x q[2];
rz(-2.1352102) q[3];
sx q[3];
rz(-1.4335535) q[3];
sx q[3];
rz(2.8999527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4686244) q[2];
sx q[2];
rz(-2.0277703) q[2];
sx q[2];
rz(1.4581397) q[2];
rz(0.80312076) q[3];
sx q[3];
rz(-1.2642658) q[3];
sx q[3];
rz(-1.7709812) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9010381) q[0];
sx q[0];
rz(-0.63837186) q[0];
sx q[0];
rz(-1.4580131) q[0];
rz(1.1395678) q[1];
sx q[1];
rz(-1.5007867) q[1];
sx q[1];
rz(2.4511888) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60421097) q[0];
sx q[0];
rz(-1.1539872) q[0];
sx q[0];
rz(-2.4875159) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5883303) q[2];
sx q[2];
rz(-0.25871535) q[2];
sx q[2];
rz(2.0257547) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0038759) q[1];
sx q[1];
rz(-2.5270577) q[1];
sx q[1];
rz(0.700386) q[1];
x q[2];
rz(-1.5059982) q[3];
sx q[3];
rz(-0.64842033) q[3];
sx q[3];
rz(2.5823926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0940501) q[2];
sx q[2];
rz(-1.8886781) q[2];
sx q[2];
rz(2.3599153) q[2];
rz(-2.9983669) q[3];
sx q[3];
rz(-0.6939303) q[3];
sx q[3];
rz(-3.1128939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52611065) q[0];
sx q[0];
rz(-0.019793864) q[0];
sx q[0];
rz(2.1671248) q[0];
rz(3.114585) q[1];
sx q[1];
rz(-1.4911431) q[1];
sx q[1];
rz(-0.19198051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8251569) q[0];
sx q[0];
rz(-0.46679206) q[0];
sx q[0];
rz(0.6760752) q[0];
rz(-2.9021743) q[2];
sx q[2];
rz(-1.7144616) q[2];
sx q[2];
rz(-1.3091376) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9908048) q[1];
sx q[1];
rz(-1.2876191) q[1];
sx q[1];
rz(-3.12451) q[1];
rz(2.8497898) q[3];
sx q[3];
rz(-2.553396) q[3];
sx q[3];
rz(0.10313973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30569046) q[2];
sx q[2];
rz(-1.4985936) q[2];
sx q[2];
rz(-0.72176019) q[2];
rz(0.52779683) q[3];
sx q[3];
rz(-0.58776394) q[3];
sx q[3];
rz(-1.2036948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0917621) q[0];
sx q[0];
rz(-1.993655) q[0];
sx q[0];
rz(0.21035305) q[0];
rz(-2.6166088) q[1];
sx q[1];
rz(-0.72768444) q[1];
sx q[1];
rz(0.64183372) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1686629) q[0];
sx q[0];
rz(-0.43420688) q[0];
sx q[0];
rz(0.28247873) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4125239) q[2];
sx q[2];
rz(-2.8260004) q[2];
sx q[2];
rz(-2.6305619) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8872309) q[1];
sx q[1];
rz(-1.6019643) q[1];
sx q[1];
rz(-1.7367669) q[1];
rz(3.1133223) q[3];
sx q[3];
rz(-1.8688374) q[3];
sx q[3];
rz(-0.6606797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.019023808) q[2];
sx q[2];
rz(-1.4539315) q[2];
sx q[2];
rz(0.34073487) q[2];
rz(-0.46180284) q[3];
sx q[3];
rz(-1.8739437) q[3];
sx q[3];
rz(-0.99892282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(3.0108764) q[0];
sx q[0];
rz(-2.0858481) q[0];
sx q[0];
rz(-2.6266932) q[0];
rz(1.8574832) q[1];
sx q[1];
rz(-1.037037) q[1];
sx q[1];
rz(-0.59891716) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4864082) q[0];
sx q[0];
rz(-2.1052971) q[0];
sx q[0];
rz(0.17528383) q[0];
x q[1];
rz(-1.0739051) q[2];
sx q[2];
rz(-0.23197939) q[2];
sx q[2];
rz(0.49099904) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45822847) q[1];
sx q[1];
rz(-2.1090611) q[1];
sx q[1];
rz(-0.12201536) q[1];
x q[2];
rz(-2.120047) q[3];
sx q[3];
rz(-1.7915939) q[3];
sx q[3];
rz(1.8017853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9942223) q[2];
sx q[2];
rz(-1.7662591) q[2];
sx q[2];
rz(0.71225172) q[2];
rz(-1.1143403) q[3];
sx q[3];
rz(-2.2388191) q[3];
sx q[3];
rz(-0.72079349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0459868) q[0];
sx q[0];
rz(-2.0144036) q[0];
sx q[0];
rz(1.4229232) q[0];
rz(-1.0320041) q[1];
sx q[1];
rz(-1.2675266) q[1];
sx q[1];
rz(1.9292018) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4044672) q[0];
sx q[0];
rz(-0.880768) q[0];
sx q[0];
rz(1.2002608) q[0];
rz(2.8396438) q[2];
sx q[2];
rz(-1.3808389) q[2];
sx q[2];
rz(-2.4446553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5281089) q[1];
sx q[1];
rz(-1.2917519) q[1];
sx q[1];
rz(2.319496) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3512896) q[3];
sx q[3];
rz(-2.7300082) q[3];
sx q[3];
rz(1.5615338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.568944) q[2];
sx q[2];
rz(-2.0062916) q[2];
sx q[2];
rz(-0.14986077) q[2];
rz(3.0585238) q[3];
sx q[3];
rz(-0.84851256) q[3];
sx q[3];
rz(-1.2868631) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6257553) q[0];
sx q[0];
rz(-2.5560162) q[0];
sx q[0];
rz(2.9199912) q[0];
rz(-1.340647) q[1];
sx q[1];
rz(-2.4602175) q[1];
sx q[1];
rz(0.63366205) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54468225) q[0];
sx q[0];
rz(-1.6297473) q[0];
sx q[0];
rz(0.037324104) q[0];
rz(-0.91434716) q[2];
sx q[2];
rz(-1.8546951) q[2];
sx q[2];
rz(1.0032297) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3089655) q[1];
sx q[1];
rz(-1.6436952) q[1];
sx q[1];
rz(0.64809191) q[1];
rz(-pi) q[2];
rz(-0.14505861) q[3];
sx q[3];
rz(-0.71509711) q[3];
sx q[3];
rz(1.3655658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0076912045) q[2];
sx q[2];
rz(-2.2654221) q[2];
sx q[2];
rz(-0.42789856) q[2];
rz(-0.89023542) q[3];
sx q[3];
rz(-2.4211703) q[3];
sx q[3];
rz(0.015457411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7240922) q[0];
sx q[0];
rz(-0.98560792) q[0];
sx q[0];
rz(0.86522657) q[0];
rz(-2.4100307) q[1];
sx q[1];
rz(-0.95844904) q[1];
sx q[1];
rz(2.1581214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18200066) q[0];
sx q[0];
rz(-1.0054709) q[0];
sx q[0];
rz(1.8046787) q[0];
rz(-pi) q[1];
rz(1.2484024) q[2];
sx q[2];
rz(-2.8722706) q[2];
sx q[2];
rz(1.4890738) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76170761) q[1];
sx q[1];
rz(-1.436625) q[1];
sx q[1];
rz(0.11685808) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9256225) q[3];
sx q[3];
rz(-2.4280818) q[3];
sx q[3];
rz(-0.3573561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5359155) q[2];
sx q[2];
rz(-1.3189877) q[2];
sx q[2];
rz(2.9299822) q[2];
rz(-2.0896046) q[3];
sx q[3];
rz(-2.7973723) q[3];
sx q[3];
rz(2.2229164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6720523) q[0];
sx q[0];
rz(-2.7993027) q[0];
sx q[0];
rz(-2.4218609) q[0];
rz(1.6795878) q[1];
sx q[1];
rz(-2.3639677) q[1];
sx q[1];
rz(3.0696312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3711187) q[0];
sx q[0];
rz(-1.2871871) q[0];
sx q[0];
rz(-0.78259612) q[0];
x q[1];
rz(2.9182704) q[2];
sx q[2];
rz(-2.2318948) q[2];
sx q[2];
rz(-1.2375792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5317412) q[1];
sx q[1];
rz(-2.688767) q[1];
sx q[1];
rz(-2.6792296) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8947745) q[3];
sx q[3];
rz(-1.9018963) q[3];
sx q[3];
rz(0.67321411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0422684) q[2];
sx q[2];
rz(-2.2734953) q[2];
sx q[2];
rz(2.457999) q[2];
rz(-1.4156703) q[3];
sx q[3];
rz(-0.96368027) q[3];
sx q[3];
rz(1.8583309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.100648) q[0];
sx q[0];
rz(-1.9293979) q[0];
sx q[0];
rz(1.9579493) q[0];
rz(-0.64507858) q[1];
sx q[1];
rz(-1.3007785) q[1];
sx q[1];
rz(0.29161463) q[1];
rz(-1.245247) q[2];
sx q[2];
rz(-1.3424046) q[2];
sx q[2];
rz(1.8142909) q[2];
rz(2.8937319) q[3];
sx q[3];
rz(-1.2206205) q[3];
sx q[3];
rz(-1.0222996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
