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
rz(-0.82395616) q[0];
sx q[0];
rz(-0.92744654) q[0];
rz(2.0916341) q[1];
sx q[1];
rz(-0.97133049) q[1];
sx q[1];
rz(-0.8144905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626635) q[0];
sx q[0];
rz(-1.9534847) q[0];
sx q[0];
rz(-2.0576059) q[0];
rz(-pi) q[1];
rz(-0.79618228) q[2];
sx q[2];
rz(-1.0999692) q[2];
sx q[2];
rz(1.7288961) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3315288) q[1];
sx q[1];
rz(-2.0671856) q[1];
sx q[1];
rz(2.1138825) q[1];
x q[2];
rz(0.82858927) q[3];
sx q[3];
rz(-2.0999206) q[3];
sx q[3];
rz(-0.55381004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77259511) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(-2.7102846) q[2];
rz(1.4946233) q[3];
sx q[3];
rz(-1.5458115) q[3];
sx q[3];
rz(0.78540426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539219) q[0];
sx q[0];
rz(-0.19667721) q[0];
sx q[0];
rz(-1.824463) q[0];
rz(1.0493086) q[1];
sx q[1];
rz(-0.53593719) q[1];
sx q[1];
rz(2.5228693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71402446) q[0];
sx q[0];
rz(-2.51665) q[0];
sx q[0];
rz(-0.24863909) q[0];
rz(0.87845408) q[2];
sx q[2];
rz(-1.8900423) q[2];
sx q[2];
rz(0.4934267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6438457) q[1];
sx q[1];
rz(-0.37393716) q[1];
sx q[1];
rz(-1.2059597) q[1];
rz(-2.5713508) q[3];
sx q[3];
rz(-0.61755731) q[3];
sx q[3];
rz(2.9887426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72545663) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(-2.6497192) q[2];
rz(-1.5623931) q[3];
sx q[3];
rz(-1.648936) q[3];
sx q[3];
rz(-0.57945848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270776) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(-2.8614817) q[0];
rz(3.0666871) q[1];
sx q[1];
rz(-2.7066051) q[1];
sx q[1];
rz(1.3704971) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48367354) q[0];
sx q[0];
rz(-2.3020805) q[0];
sx q[0];
rz(0.28261225) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2417775) q[2];
sx q[2];
rz(-1.0428937) q[2];
sx q[2];
rz(1.6646651) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31077267) q[1];
sx q[1];
rz(-1.8337063) q[1];
sx q[1];
rz(-0.57057686) q[1];
rz(-pi) q[2];
rz(-3.0093569) q[3];
sx q[3];
rz(-1.5521129) q[3];
sx q[3];
rz(-1.4052237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28477731) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(0.80186296) q[2];
rz(2.2297468) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(0.97889939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80482471) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(-0.62070745) q[0];
rz(0.2785109) q[1];
sx q[1];
rz(-1.4418437) q[1];
sx q[1];
rz(-1.0034358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28351706) q[0];
sx q[0];
rz(-1.9122313) q[0];
sx q[0];
rz(0.33591875) q[0];
x q[1];
rz(-0.070095991) q[2];
sx q[2];
rz(-1.9906487) q[2];
sx q[2];
rz(1.1481089) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.089503376) q[1];
sx q[1];
rz(-1.6253396) q[1];
sx q[1];
rz(-0.2691582) q[1];
rz(-pi) q[2];
rz(-1.8068321) q[3];
sx q[3];
rz(-0.98525648) q[3];
sx q[3];
rz(-1.3957617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3044546) q[2];
sx q[2];
rz(-2.0292323) q[2];
sx q[2];
rz(3.0755074) q[2];
rz(-2.0153913) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(2.5662305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5758301) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(2.3542812) q[0];
rz(0.85834223) q[1];
sx q[1];
rz(-0.97802496) q[1];
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
rz(-2.4466956) q[0];
x q[1];
rz(-2.0400286) q[2];
sx q[2];
rz(-1.9479556) q[2];
sx q[2];
rz(-2.9531329) q[2];
rz(-pi) q[3];
x q[3];
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
rz(-pi) q[2];
rz(-0.85162195) q[3];
sx q[3];
rz(-1.1332261) q[3];
sx q[3];
rz(2.9629666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9786238) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(-0.84257379) q[2];
rz(-2.8652371) q[3];
sx q[3];
rz(-0.98139757) q[3];
sx q[3];
rz(-1.792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2862947) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(-1.1899765) q[0];
rz(-1.2225993) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(2.55866) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8609498) q[0];
sx q[0];
rz(-1.1034729) q[0];
sx q[0];
rz(-1.1790465) q[0];
rz(-pi) q[1];
rz(3.1276064) q[2];
sx q[2];
rz(-2.1443233) q[2];
sx q[2];
rz(-1.397246) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1057824) q[1];
sx q[1];
rz(-1.3329734) q[1];
sx q[1];
rz(1.4877122) q[1];
rz(-pi) q[2];
rz(1.1357901) q[3];
sx q[3];
rz(-0.9851176) q[3];
sx q[3];
rz(-0.034736573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0588093) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(-1.0397376) q[2];
rz(0.58935634) q[3];
sx q[3];
rz(-1.6683234) q[3];
sx q[3];
rz(-1.0065494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.0712737) q[1];
sx q[1];
rz(0.76990661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8827646) q[0];
sx q[0];
rz(-2.6443091) q[0];
sx q[0];
rz(0.33035853) q[0];
x q[1];
rz(-2.3138397) q[2];
sx q[2];
rz(-2.0619287) q[2];
sx q[2];
rz(-0.14358768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47178956) q[1];
sx q[1];
rz(-0.21179971) q[1];
sx q[1];
rz(2.7752911) q[1];
rz(-1.4925333) q[3];
sx q[3];
rz(-2.1746082) q[3];
sx q[3];
rz(-0.1426689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.7485488) q[2];
sx q[2];
rz(-1.0153655) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2400951) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(-2.748306) q[0];
rz(-2.6914864) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(-3.111305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05370985) q[0];
sx q[0];
rz(-1.6018512) q[0];
sx q[0];
rz(2.8939918) q[0];
rz(-0.25595902) q[2];
sx q[2];
rz(-2.3156907) q[2];
sx q[2];
rz(1.9108552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4612164) q[1];
sx q[1];
rz(-2.595397) q[1];
sx q[1];
rz(3.097018) q[1];
rz(-1.7025331) q[3];
sx q[3];
rz(-1.3088969) q[3];
sx q[3];
rz(-1.1887664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1071189) q[2];
sx q[2];
rz(-1.5849042) q[2];
sx q[2];
rz(1.1220453) q[2];
rz(-2.0864887) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(-1.25157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.0301544) q[1];
sx q[1];
rz(-1.0672306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8071496) q[0];
sx q[0];
rz(-2.4031638) q[0];
sx q[0];
rz(3.1160808) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8781233) q[2];
sx q[2];
rz(-1.5241703) q[2];
sx q[2];
rz(1.7152552) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9561056) q[1];
sx q[1];
rz(-1.8813349) q[1];
sx q[1];
rz(2.8794672) q[1];
x q[2];
rz(-2.9638941) q[3];
sx q[3];
rz(-2.0434994) q[3];
sx q[3];
rz(2.3426169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5344703) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(1.8758476) q[2];
rz(-0.040146116) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-0.039073959) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48113111) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(-1.2598502) q[0];
rz(-1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(3.1390417) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77230763) q[0];
sx q[0];
rz(-2.9678759) q[0];
sx q[0];
rz(0.80887167) q[0];
x q[1];
rz(3.1319982) q[2];
sx q[2];
rz(-1.506732) q[2];
sx q[2];
rz(1.4256776) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5469909) q[1];
sx q[1];
rz(-0.78557116) q[1];
sx q[1];
rz(0.57348302) q[1];
rz(-pi) q[2];
rz(1.5137767) q[3];
sx q[3];
rz(-0.43511639) q[3];
sx q[3];
rz(0.42729898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5551787) q[2];
sx q[2];
rz(-1.5191583) q[2];
sx q[2];
rz(0.67443887) q[2];
rz(2.0241578) q[3];
sx q[3];
rz(-1.0267461) q[3];
sx q[3];
rz(0.31291541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085070327) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(-2.7936735) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(1.3580218) q[2];
sx q[2];
rz(-1.2322896) q[2];
sx q[2];
rz(-0.88862669) q[2];
rz(-0.017853768) q[3];
sx q[3];
rz(-1.4244867) q[3];
sx q[3];
rz(-1.3467237) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
