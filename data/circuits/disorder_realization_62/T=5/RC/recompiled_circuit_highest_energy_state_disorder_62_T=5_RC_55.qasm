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
rz(-1.6031185) q[0];
sx q[0];
rz(-1.7552019) q[0];
sx q[0];
rz(-1.8829128) q[0];
rz(1.4409244) q[1];
sx q[1];
rz(-2.6135593) q[1];
sx q[1];
rz(0.33906403) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10701943) q[0];
sx q[0];
rz(-1.9030754) q[0];
sx q[0];
rz(-2.9884724) q[0];
rz(-pi) q[1];
rz(-1.1669615) q[2];
sx q[2];
rz(-0.93138715) q[2];
sx q[2];
rz(-1.0013231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.69341371) q[1];
sx q[1];
rz(-2.1817345) q[1];
sx q[1];
rz(-3.0603859) q[1];
rz(-pi) q[2];
rz(-1.6293832) q[3];
sx q[3];
rz(-2.226429) q[3];
sx q[3];
rz(0.4200622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49860927) q[2];
sx q[2];
rz(-0.53350854) q[2];
sx q[2];
rz(-1.3996997) q[2];
rz(1.1534322) q[3];
sx q[3];
rz(-1.6250936) q[3];
sx q[3];
rz(-1.4199408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145849) q[0];
sx q[0];
rz(-1.4684533) q[0];
sx q[0];
rz(0.46838316) q[0];
rz(-0.98764694) q[1];
sx q[1];
rz(-1.7665607) q[1];
sx q[1];
rz(1.04029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9355277) q[0];
sx q[0];
rz(-2.1305232) q[0];
sx q[0];
rz(-2.2846643) q[0];
rz(-pi) q[1];
rz(-0.78590337) q[2];
sx q[2];
rz(-2.6214947) q[2];
sx q[2];
rz(1.4532068) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2548488) q[1];
sx q[1];
rz(-0.4720531) q[1];
sx q[1];
rz(2.9242502) q[1];
rz(-pi) q[2];
rz(1.0063824) q[3];
sx q[3];
rz(-1.4335535) q[3];
sx q[3];
rz(2.8999527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6729683) q[2];
sx q[2];
rz(-1.1138223) q[2];
sx q[2];
rz(1.4581397) q[2];
rz(0.80312076) q[3];
sx q[3];
rz(-1.8773269) q[3];
sx q[3];
rz(1.7709812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9010381) q[0];
sx q[0];
rz(-2.5032208) q[0];
sx q[0];
rz(1.6835796) q[0];
rz(1.1395678) q[1];
sx q[1];
rz(-1.640806) q[1];
sx q[1];
rz(-2.4511888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60421097) q[0];
sx q[0];
rz(-1.9876055) q[0];
sx q[0];
rz(0.65407674) q[0];
rz(1.312119) q[2];
sx q[2];
rz(-1.575282) q[2];
sx q[2];
rz(2.6696837) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1377167) q[1];
sx q[1];
rz(-0.614535) q[1];
sx q[1];
rz(-2.4412066) q[1];
rz(-pi) q[2];
rz(-3.0925678) q[3];
sx q[3];
rz(-0.92396523) q[3];
sx q[3];
rz(-2.5011592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0475426) q[2];
sx q[2];
rz(-1.8886781) q[2];
sx q[2];
rz(-2.3599153) q[2];
rz(2.9983669) q[3];
sx q[3];
rz(-2.4476624) q[3];
sx q[3];
rz(-3.1128939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52611065) q[0];
sx q[0];
rz(-3.1217988) q[0];
sx q[0];
rz(-2.1671248) q[0];
rz(-3.114585) q[1];
sx q[1];
rz(-1.6504495) q[1];
sx q[1];
rz(2.9496121) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8251569) q[0];
sx q[0];
rz(-2.6748006) q[0];
sx q[0];
rz(-2.4655174) q[0];
rz(0.23941834) q[2];
sx q[2];
rz(-1.4271311) q[2];
sx q[2];
rz(-1.832455) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0518735) q[1];
sx q[1];
rz(-0.28367821) q[1];
sx q[1];
rz(1.62943) q[1];
x q[2];
rz(0.56842808) q[3];
sx q[3];
rz(-1.7311043) q[3];
sx q[3];
rz(-1.4290546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8359022) q[2];
sx q[2];
rz(-1.6429991) q[2];
sx q[2];
rz(2.4198325) q[2];
rz(-0.52779683) q[3];
sx q[3];
rz(-0.58776394) q[3];
sx q[3];
rz(1.2036948) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0917621) q[0];
sx q[0];
rz(-1.1479377) q[0];
sx q[0];
rz(-2.9312396) q[0];
rz(-0.52498388) q[1];
sx q[1];
rz(-0.72768444) q[1];
sx q[1];
rz(-0.64183372) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1686629) q[0];
sx q[0];
rz(-2.7073858) q[0];
sx q[0];
rz(-2.8591139) q[0];
rz(-pi) q[1];
rz(1.3566229) q[2];
sx q[2];
rz(-1.3371981) q[2];
sx q[2];
rz(2.8983299) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8199379) q[1];
sx q[1];
rz(-1.4049071) q[1];
sx q[1];
rz(3.1099907) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.028270311) q[3];
sx q[3];
rz(-1.8688374) q[3];
sx q[3];
rz(2.480913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.019023808) q[2];
sx q[2];
rz(-1.6876612) q[2];
sx q[2];
rz(0.34073487) q[2];
rz(-0.46180284) q[3];
sx q[3];
rz(-1.267649) q[3];
sx q[3];
rz(0.99892282) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13071624) q[0];
sx q[0];
rz(-1.0557446) q[0];
sx q[0];
rz(2.6266932) q[0];
rz(1.2841094) q[1];
sx q[1];
rz(-1.037037) q[1];
sx q[1];
rz(-2.5426755) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8209957) q[0];
sx q[0];
rz(-0.55984523) q[0];
sx q[0];
rz(1.8572771) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3660422) q[2];
sx q[2];
rz(-1.6806112) q[2];
sx q[2];
rz(1.5653953) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22344124) q[1];
sx q[1];
rz(-2.591003) q[1];
sx q[1];
rz(-1.3697025) q[1];
rz(1.9768561) q[3];
sx q[3];
rz(-0.58771509) q[3];
sx q[3];
rz(0.57460659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9942223) q[2];
sx q[2];
rz(-1.3753336) q[2];
sx q[2];
rz(-0.71225172) q[2];
rz(2.0272523) q[3];
sx q[3];
rz(-2.2388191) q[3];
sx q[3];
rz(2.4207992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095605843) q[0];
sx q[0];
rz(-2.0144036) q[0];
sx q[0];
rz(1.4229232) q[0];
rz(1.0320041) q[1];
sx q[1];
rz(-1.874066) q[1];
sx q[1];
rz(1.9292018) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7371255) q[0];
sx q[0];
rz(-2.2608247) q[0];
sx q[0];
rz(1.9413319) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5676377) q[2];
sx q[2];
rz(-0.35517737) q[2];
sx q[2];
rz(-1.7227625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29323762) q[1];
sx q[1];
rz(-0.85746409) q[1];
sx q[1];
rz(2.7687702) q[1];
rz(0.29799975) q[3];
sx q[3];
rz(-1.8590419) q[3];
sx q[3];
rz(0.73747613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.57264868) q[2];
sx q[2];
rz(-1.135301) q[2];
sx q[2];
rz(0.14986077) q[2];
rz(-0.0830689) q[3];
sx q[3];
rz(-0.84851256) q[3];
sx q[3];
rz(1.8547295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6257553) q[0];
sx q[0];
rz(-2.5560162) q[0];
sx q[0];
rz(0.22160141) q[0];
rz(-1.8009456) q[1];
sx q[1];
rz(-2.4602175) q[1];
sx q[1];
rz(2.5079306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54468225) q[0];
sx q[0];
rz(-1.6297473) q[0];
sx q[0];
rz(-3.1042685) q[0];
rz(-pi) q[1];
rz(1.1248323) q[2];
sx q[2];
rz(-2.4348091) q[2];
sx q[2];
rz(0.91627322) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8326272) q[1];
sx q[1];
rz(-1.6436952) q[1];
sx q[1];
rz(0.64809191) q[1];
x q[2];
rz(0.70987411) q[3];
sx q[3];
rz(-1.6657192) q[3];
sx q[3];
rz(3.0462163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1339014) q[2];
sx q[2];
rz(-2.2654221) q[2];
sx q[2];
rz(2.7136941) q[2];
rz(-0.89023542) q[3];
sx q[3];
rz(-2.4211703) q[3];
sx q[3];
rz(0.015457411) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4175005) q[0];
sx q[0];
rz(-2.1559847) q[0];
sx q[0];
rz(-0.86522657) q[0];
rz(0.73156196) q[1];
sx q[1];
rz(-0.95844904) q[1];
sx q[1];
rz(2.1581214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18200066) q[0];
sx q[0];
rz(-2.1361217) q[0];
sx q[0];
rz(-1.336914) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3147361) q[2];
sx q[2];
rz(-1.4863925) q[2];
sx q[2];
rz(-2.7483482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76170761) q[1];
sx q[1];
rz(-1.436625) q[1];
sx q[1];
rz(-3.0247346) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2526475) q[3];
sx q[3];
rz(-1.3414012) q[3];
sx q[3];
rz(-1.4865686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5359155) q[2];
sx q[2];
rz(-1.8226049) q[2];
sx q[2];
rz(-2.9299822) q[2];
rz(-1.051988) q[3];
sx q[3];
rz(-0.34422031) q[3];
sx q[3];
rz(-0.91867623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4695404) q[0];
sx q[0];
rz(-2.7993027) q[0];
sx q[0];
rz(2.4218609) q[0];
rz(1.6795878) q[1];
sx q[1];
rz(-0.77762496) q[1];
sx q[1];
rz(0.07196149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6126637) q[0];
sx q[0];
rz(-0.82717973) q[0];
sx q[0];
rz(-1.1808047) q[0];
rz(-1.2934204) q[2];
sx q[2];
rz(-0.69239834) q[2];
sx q[2];
rz(-1.5497249) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60985145) q[1];
sx q[1];
rz(-2.688767) q[1];
sx q[1];
rz(2.6792296) q[1];
x q[2];
rz(-2.7937204) q[3];
sx q[3];
rz(-1.265003) q[3];
sx q[3];
rz(2.3527462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0422684) q[2];
sx q[2];
rz(-2.2734953) q[2];
sx q[2];
rz(0.68359366) q[2];
rz(1.7259224) q[3];
sx q[3];
rz(-2.1779124) q[3];
sx q[3];
rz(1.2832618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.245247) q[2];
sx q[2];
rz(-1.799188) q[2];
sx q[2];
rz(-1.3273018) q[2];
rz(1.2105) q[3];
sx q[3];
rz(-1.8033284) q[3];
sx q[3];
rz(-2.506496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
