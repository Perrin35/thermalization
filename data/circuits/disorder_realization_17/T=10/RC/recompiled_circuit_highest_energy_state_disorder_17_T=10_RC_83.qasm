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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0448582) q[0];
sx q[0];
rz(-2.0196901) q[0];
sx q[0];
rz(0.4273703) q[0];
rz(2.3454104) q[2];
sx q[2];
rz(-1.0999692) q[2];
sx q[2];
rz(-1.4126965) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.51920096) q[1];
sx q[1];
rz(-2.0425052) q[1];
sx q[1];
rz(-0.5640819) q[1];
x q[2];
rz(-0.82858927) q[3];
sx q[3];
rz(-1.0416721) q[3];
sx q[3];
rz(-0.55381004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3689975) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(-2.7102846) q[2];
rz(-1.6469693) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(-0.78540426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5539219) q[0];
sx q[0];
rz(-0.19667721) q[0];
sx q[0];
rz(-1.3171296) q[0];
rz(1.0493086) q[1];
sx q[1];
rz(-2.6056555) q[1];
sx q[1];
rz(-2.5228693) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65370377) q[0];
sx q[0];
rz(-1.426322) q[0];
sx q[0];
rz(-2.5313951) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0929865) q[2];
sx q[2];
rz(-2.3903762) q[2];
sx q[2];
rz(2.4260245) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2545027) q[1];
sx q[1];
rz(-1.2225593) q[1];
sx q[1];
rz(3.00249) q[1];
rz(2.5713508) q[3];
sx q[3];
rz(-0.61755731) q[3];
sx q[3];
rz(0.15285003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.416136) q[2];
sx q[2];
rz(-1.340103) q[2];
sx q[2];
rz(-2.6497192) q[2];
rz(1.5623931) q[3];
sx q[3];
rz(-1.648936) q[3];
sx q[3];
rz(0.57945848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270776) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(-2.8614817) q[0];
rz(-3.0666871) q[1];
sx q[1];
rz(-2.7066051) q[1];
sx q[1];
rz(-1.3704971) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460187) q[0];
sx q[0];
rz(-1.361712) q[0];
sx q[0];
rz(0.81935291) q[0];
x q[1];
rz(-0.89981516) q[2];
sx q[2];
rz(-1.0428937) q[2];
sx q[2];
rz(1.6646651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0468415) q[1];
sx q[1];
rz(-1.0221204) q[1];
sx q[1];
rz(1.8803174) q[1];
rz(-pi) q[2];
rz(0.14078056) q[3];
sx q[3];
rz(-0.13354145) q[3];
sx q[3];
rz(-2.8364766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8568153) q[2];
sx q[2];
rz(-1.6723526) q[2];
sx q[2];
rz(2.3397297) q[2];
rz(-0.91184584) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(-2.1626933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3367679) q[0];
sx q[0];
rz(-0.66773361) q[0];
sx q[0];
rz(-0.62070745) q[0];
rz(-0.2785109) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(2.1381569) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1708978) q[0];
sx q[0];
rz(-1.8866294) q[0];
sx q[0];
rz(1.210808) q[0];
x q[1];
rz(0.070095991) q[2];
sx q[2];
rz(-1.9906487) q[2];
sx q[2];
rz(-1.1481089) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0520893) q[1];
sx q[1];
rz(-1.516253) q[1];
sx q[1];
rz(2.8724345) q[1];
x q[2];
rz(-2.5430319) q[3];
sx q[3];
rz(-1.7669456) q[3];
sx q[3];
rz(2.83441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8371381) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(-0.066085286) q[2];
rz(2.0153913) q[3];
sx q[3];
rz(-0.76032138) q[3];
sx q[3];
rz(-0.57536212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5758301) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(-2.3542812) q[0];
rz(-2.2832504) q[1];
sx q[1];
rz(-2.1635677) q[1];
sx q[1];
rz(1.0816921) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72321945) q[0];
sx q[0];
rz(-2.0670509) q[0];
sx q[0];
rz(2.4346274) q[0];
rz(-pi) q[1];
rz(0.85140793) q[2];
sx q[2];
rz(-2.5485196) q[2];
sx q[2];
rz(-0.75411686) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8199205) q[1];
sx q[1];
rz(-1.5645488) q[1];
sx q[1];
rz(-2.856918) q[1];
x q[2];
rz(0.55629913) q[3];
sx q[3];
rz(-0.93141684) q[3];
sx q[3];
rz(2.104708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9786238) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(-0.84257379) q[2];
rz(2.8652371) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.855298) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(1.9516161) q[0];
rz(-1.2225993) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(2.55866) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8609498) q[0];
sx q[0];
rz(-2.0381198) q[0];
sx q[0];
rz(1.1790465) q[0];
rz(-pi) q[1];
x q[1];
rz(0.013986258) q[2];
sx q[2];
rz(-2.1443233) q[2];
sx q[2];
rz(1.397246) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.035810245) q[1];
sx q[1];
rz(-1.3329734) q[1];
sx q[1];
rz(1.4877122) q[1];
rz(-0.56598466) q[3];
sx q[3];
rz(-0.71403394) q[3];
sx q[3];
rz(-0.66431724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0827834) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(-1.0397376) q[2];
rz(-0.58935634) q[3];
sx q[3];
rz(-1.6683234) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650918) q[0];
sx q[0];
rz(-0.14987513) q[0];
sx q[0];
rz(1.804922) q[0];
rz(-0.22557766) q[1];
sx q[1];
rz(-2.070319) q[1];
sx q[1];
rz(-0.76990661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8827646) q[0];
sx q[0];
rz(-0.49728359) q[0];
sx q[0];
rz(0.33035853) q[0];
rz(-pi) q[1];
rz(-0.628148) q[2];
sx q[2];
rz(-0.9315812) q[2];
sx q[2];
rz(-2.1232427) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47178956) q[1];
sx q[1];
rz(-2.9297929) q[1];
sx q[1];
rz(-2.7752911) q[1];
rz(-pi) q[2];
rz(-3.0287241) q[3];
sx q[3];
rz(-0.60823554) q[3];
sx q[3];
rz(3.1361767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(-2.1262271) q[2];
rz(-0.20396248) q[3];
sx q[3];
rz(-1.5927619) q[3];
sx q[3];
rz(1.7502194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2400951) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(2.748306) q[0];
rz(-0.45010629) q[1];
sx q[1];
rz(-1.422912) q[1];
sx q[1];
rz(0.030287655) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5092376) q[0];
sx q[0];
rz(-1.3233174) q[0];
sx q[0];
rz(1.6028274) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3322163) q[2];
sx q[2];
rz(-1.7580108) q[2];
sx q[2];
rz(2.9771114) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4612164) q[1];
sx q[1];
rz(-2.595397) q[1];
sx q[1];
rz(3.097018) q[1];
x q[2];
rz(2.6859517) q[3];
sx q[3];
rz(-0.29248387) q[3];
sx q[3];
rz(0.71575056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0344737) q[2];
sx q[2];
rz(-1.5849042) q[2];
sx q[2];
rz(2.0195473) q[2];
rz(2.0864887) q[3];
sx q[3];
rz(-0.3873581) q[3];
sx q[3];
rz(1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9686862) q[0];
sx q[0];
rz(-2.1538669) q[0];
sx q[0];
rz(-2.3671142) q[0];
rz(-2.5075746) q[1];
sx q[1];
rz(-2.1114383) q[1];
sx q[1];
rz(-2.074362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2174847) q[0];
sx q[0];
rz(-1.5536246) q[0];
sx q[0];
rz(2.4033258) q[0];
rz(-pi) q[1];
rz(1.7238317) q[2];
sx q[2];
rz(-2.8308597) q[2];
sx q[2];
rz(2.8513214) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9561056) q[1];
sx q[1];
rz(-1.8813349) q[1];
sx q[1];
rz(0.26212543) q[1];
x q[2];
rz(-1.2379856) q[3];
sx q[3];
rz(-0.5026256) q[3];
sx q[3];
rz(-0.42325936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5344703) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(-1.8758476) q[2];
rz(-0.040146116) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-0.039073959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6604615) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(-1.2598502) q[0];
rz(-1.4866359) q[1];
sx q[1];
rz(-1.0340374) q[1];
sx q[1];
rz(3.1390417) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.369285) q[0];
sx q[0];
rz(-0.17371674) q[0];
sx q[0];
rz(2.332721) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6348636) q[2];
sx q[2];
rz(-1.5803711) q[2];
sx q[2];
rz(0.14450444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5469909) q[1];
sx q[1];
rz(-2.3560215) q[1];
sx q[1];
rz(-2.5681096) q[1];
rz(-pi) q[2];
rz(-1.5137767) q[3];
sx q[3];
rz(-2.7064763) q[3];
sx q[3];
rz(0.42729898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.586414) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(0.67443887) q[2];
rz(-1.1174348) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(-0.31291541) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085070327) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(0.3479192) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(2.7958776) q[2];
sx q[2];
rz(-1.7713265) q[2];
sx q[2];
rz(-2.3878018) q[2];
rz(-1.717129) q[3];
sx q[3];
rz(-1.5884593) q[3];
sx q[3];
rz(-2.9149169) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
