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
rz(1.849527) q[0];
sx q[0];
rz(5.4592291) q[0];
sx q[0];
rz(8.4973314) q[0];
rz(2.0916341) q[1];
sx q[1];
rz(-0.97133049) q[1];
sx q[1];
rz(-0.8144905) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096734418) q[0];
sx q[0];
rz(-1.1219026) q[0];
sx q[0];
rz(2.7142224) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61887069) q[2];
sx q[2];
rz(-0.89779389) q[2];
sx q[2];
rz(0.25970632) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.42871745) q[1];
sx q[1];
rz(-2.4230885) q[1];
sx q[1];
rz(-0.76193049) q[1];
rz(-0.82858927) q[3];
sx q[3];
rz(-2.0999206) q[3];
sx q[3];
rz(0.55381004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3689975) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(2.7102846) q[2];
rz(1.6469693) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(-2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5876708) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(-1.824463) q[0];
rz(-1.0493086) q[1];
sx q[1];
rz(-2.6056555) q[1];
sx q[1];
rz(2.5228693) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1241708) q[0];
sx q[0];
rz(-0.96786495) q[0];
sx q[0];
rz(1.3951016) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2631386) q[2];
sx q[2];
rz(-1.8900423) q[2];
sx q[2];
rz(2.648166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2545027) q[1];
sx q[1];
rz(-1.2225593) q[1];
sx q[1];
rz(-0.13910267) q[1];
x q[2];
rz(2.5713508) q[3];
sx q[3];
rz(-0.61755731) q[3];
sx q[3];
rz(-2.9887426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.416136) q[2];
sx q[2];
rz(-1.340103) q[2];
sx q[2];
rz(0.49187342) q[2];
rz(-1.5791996) q[3];
sx q[3];
rz(-1.648936) q[3];
sx q[3];
rz(0.57945848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270776) q[0];
sx q[0];
rz(-1.0579695) q[0];
sx q[0];
rz(2.8614817) q[0];
rz(-3.0666871) q[1];
sx q[1];
rz(-0.43498755) q[1];
sx q[1];
rz(-1.7710955) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6579191) q[0];
sx q[0];
rz(-2.3020805) q[0];
sx q[0];
rz(2.8589804) q[0];
rz(-2.5016194) q[2];
sx q[2];
rz(-1.0037862) q[2];
sx q[2];
rz(-2.8550573) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.83082) q[1];
sx q[1];
rz(-1.8337063) q[1];
sx q[1];
rz(-0.57057686) q[1];
x q[2];
rz(-3.0093569) q[3];
sx q[3];
rz(-1.5894798) q[3];
sx q[3];
rz(1.4052237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28477731) q[2];
sx q[2];
rz(-1.6723526) q[2];
sx q[2];
rz(0.80186296) q[2];
rz(-2.2297468) q[3];
sx q[3];
rz(-0.82404476) q[3];
sx q[3];
rz(0.97889939) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3367679) q[0];
sx q[0];
rz(-0.66773361) q[0];
sx q[0];
rz(-0.62070745) q[0];
rz(0.2785109) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(1.0034358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8580756) q[0];
sx q[0];
rz(-1.2293613) q[0];
sx q[0];
rz(0.33591875) q[0];
x q[1];
rz(-1.7264257) q[2];
sx q[2];
rz(-0.4253201) q[2];
sx q[2];
rz(0.97753866) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.675337) q[1];
sx q[1];
rz(-1.8395443) q[1];
sx q[1];
rz(1.5142201) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5430319) q[3];
sx q[3];
rz(-1.374647) q[3];
sx q[3];
rz(2.83441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3044546) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(-3.0755074) q[2];
rz(2.0153913) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(-2.5662305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5657625) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(2.3542812) q[0];
rz(2.2832504) q[1];
sx q[1];
rz(-2.1635677) q[1];
sx q[1];
rz(2.0599005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3560548) q[0];
sx q[0];
rz(-0.83844664) q[0];
sx q[0];
rz(0.69489702) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7236347) q[2];
sx q[2];
rz(-1.1368898) q[2];
sx q[2];
rz(1.5746631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2277775) q[1];
sx q[1];
rz(-2.8568513) q[1];
sx q[1];
rz(-0.022241994) q[1];
rz(0.55629913) q[3];
sx q[3];
rz(-0.93141684) q[3];
sx q[3];
rz(-1.0368846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9786238) q[2];
sx q[2];
rz(-1.9917515) q[2];
sx q[2];
rz(-0.84257379) q[2];
rz(-0.27635559) q[3];
sx q[3];
rz(-2.1601951) q[3];
sx q[3];
rz(-1.792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.855298) q[0];
sx q[0];
rz(-1.3110302) q[0];
sx q[0];
rz(1.9516161) q[0];
rz(1.9189934) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(-0.58293265) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2806429) q[0];
sx q[0];
rz(-2.0381198) q[0];
sx q[0];
rz(-1.9625462) q[0];
x q[1];
rz(0.013986258) q[2];
sx q[2];
rz(-2.1443233) q[2];
sx q[2];
rz(1.397246) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5546023) q[1];
sx q[1];
rz(-1.6515367) q[1];
sx q[1];
rz(-0.23861532) q[1];
rz(-pi) q[2];
rz(-0.63151367) q[3];
sx q[3];
rz(-1.9296292) q[3];
sx q[3];
rz(-1.3541019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0588093) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(1.0397376) q[2];
rz(0.58935634) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(1.0065494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37650087) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(1.3366706) q[0];
rz(-2.916015) q[1];
sx q[1];
rz(-2.070319) q[1];
sx q[1];
rz(-2.371686) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8827646) q[0];
sx q[0];
rz(-2.6443091) q[0];
sx q[0];
rz(-0.33035853) q[0];
rz(-pi) q[1];
rz(2.3138397) q[2];
sx q[2];
rz(-2.0619287) q[2];
sx q[2];
rz(-2.998005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6837782) q[1];
sx q[1];
rz(-1.646161) q[1];
sx q[1];
rz(0.1981258) q[1];
x q[2];
rz(-0.60524551) q[3];
sx q[3];
rz(-1.5063933) q[3];
sx q[3];
rz(-1.4726255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62468195) q[2];
sx q[2];
rz(-1.7485488) q[2];
sx q[2];
rz(1.0153655) q[2];
rz(-2.9376302) q[3];
sx q[3];
rz(-1.5927619) q[3];
sx q[3];
rz(1.3913733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90149752) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(0.39328662) q[0];
rz(-2.6914864) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(0.030287655) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0878828) q[0];
sx q[0];
rz(-1.6018512) q[0];
sx q[0];
rz(2.8939918) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8387536) q[2];
sx q[2];
rz(-2.3620124) q[2];
sx q[2];
rz(1.5991581) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.62822484) q[1];
sx q[1];
rz(-1.0252044) q[1];
sx q[1];
rz(1.5437158) q[1];
x q[2];
rz(2.8775086) q[3];
sx q[3];
rz(-1.698016) q[3];
sx q[3];
rz(-2.7252687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1071189) q[2];
sx q[2];
rz(-1.5849042) q[2];
sx q[2];
rz(-2.0195473) q[2];
rz(-2.0864887) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9686862) q[0];
sx q[0];
rz(-2.1538669) q[0];
sx q[0];
rz(-2.3671142) q[0];
rz(2.5075746) q[1];
sx q[1];
rz(-2.1114383) q[1];
sx q[1];
rz(-1.0672306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7726583) q[0];
sx q[0];
rz(-0.83266363) q[0];
sx q[0];
rz(-1.5940109) q[0];
rz(0.048914476) q[2];
sx q[2];
rz(-1.8777784) q[2];
sx q[2];
rz(0.129667) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18548705) q[1];
sx q[1];
rz(-1.8813349) q[1];
sx q[1];
rz(0.26212543) q[1];
rz(-pi) q[2];
rz(2.9638941) q[3];
sx q[3];
rz(-1.0980933) q[3];
sx q[3];
rz(-0.79897579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5344703) q[2];
sx q[2];
rz(-3.0202713) q[2];
sx q[2];
rz(-1.2657451) q[2];
rz(3.1014465) q[3];
sx q[3];
rz(-2.7669192) q[3];
sx q[3];
rz(0.039073959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.6604615) q[0];
sx q[0];
rz(-0.88812319) q[0];
sx q[0];
rz(1.2598502) q[0];
rz(1.4866359) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(3.1390417) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5887506) q[0];
sx q[0];
rz(-1.6903983) q[0];
sx q[0];
rz(1.6970859) q[0];
rz(3.1319982) q[2];
sx q[2];
rz(-1.6348607) q[2];
sx q[2];
rz(1.715915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3348744) q[1];
sx q[1];
rz(-0.93466991) q[1];
sx q[1];
rz(-2.0680554) q[1];
x q[2];
rz(1.5137767) q[3];
sx q[3];
rz(-2.7064763) q[3];
sx q[3];
rz(2.7142937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5551787) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(-2.4671538) q[2];
rz(-2.0241578) q[3];
sx q[3];
rz(-1.0267461) q[3];
sx q[3];
rz(-0.31291541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0565223) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(0.3479192) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(-1.3580218) q[2];
sx q[2];
rz(-1.909303) q[2];
sx q[2];
rz(2.252966) q[2];
rz(0.017853768) q[3];
sx q[3];
rz(-1.717106) q[3];
sx q[3];
rz(1.794869) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
