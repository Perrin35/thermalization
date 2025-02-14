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
rz(-1.0499586) q[1];
sx q[1];
rz(-2.1702622) q[1];
sx q[1];
rz(0.8144905) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0448582) q[0];
sx q[0];
rz(-2.0196901) q[0];
sx q[0];
rz(2.7142224) q[0];
rz(2.1998946) q[2];
sx q[2];
rz(-0.88028958) q[2];
sx q[2];
rz(-2.5494573) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42871745) q[1];
sx q[1];
rz(-0.71850417) q[1];
sx q[1];
rz(-0.76193049) q[1];
rz(-pi) q[2];
rz(-0.82858927) q[3];
sx q[3];
rz(-2.0999206) q[3];
sx q[3];
rz(-2.5877826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3689975) q[2];
sx q[2];
rz(-1.3694222) q[2];
sx q[2];
rz(-0.43130809) q[2];
rz(1.6469693) q[3];
sx q[3];
rz(-1.5458115) q[3];
sx q[3];
rz(2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539219) q[0];
sx q[0];
rz(-0.19667721) q[0];
sx q[0];
rz(1.3171296) q[0];
rz(2.092284) q[1];
sx q[1];
rz(-0.53593719) q[1];
sx q[1];
rz(-2.5228693) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174219) q[0];
sx q[0];
rz(-0.96786495) q[0];
sx q[0];
rz(-1.3951016) q[0];
rz(-2.7359782) q[2];
sx q[2];
rz(-2.2219293) q[2];
sx q[2];
rz(1.8096015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8730388) q[1];
sx q[1];
rz(-1.4400926) q[1];
sx q[1];
rz(1.2194344) q[1];
rz(-pi) q[2];
rz(2.6027565) q[3];
sx q[3];
rz(-1.2528786) q[3];
sx q[3];
rz(-0.93618083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.72545663) q[2];
sx q[2];
rz(-1.340103) q[2];
sx q[2];
rz(-2.6497192) q[2];
rz(-1.5791996) q[3];
sx q[3];
rz(-1.4926566) q[3];
sx q[3];
rz(2.5621342) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48367354) q[0];
sx q[0];
rz(-0.83951211) q[0];
sx q[0];
rz(0.28261225) q[0];
rz(-pi) q[1];
rz(0.81746705) q[2];
sx q[2];
rz(-0.82767476) q[2];
sx q[2];
rz(2.4823014) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4968003) q[1];
sx q[1];
rz(-0.62207568) q[1];
sx q[1];
rz(-2.679307) q[1];
rz(-pi) q[2];
rz(-3.0008121) q[3];
sx q[3];
rz(-0.13354145) q[3];
sx q[3];
rz(0.30511607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8568153) q[2];
sx q[2];
rz(-1.6723526) q[2];
sx q[2];
rz(0.80186296) q[2];
rz(-0.91184584) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(0.97889939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80482471) q[0];
sx q[0];
rz(-0.66773361) q[0];
sx q[0];
rz(0.62070745) q[0];
rz(-2.8630818) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(-2.1381569) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8580756) q[0];
sx q[0];
rz(-1.2293613) q[0];
sx q[0];
rz(2.8056739) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.070095991) q[2];
sx q[2];
rz(-1.150944) q[2];
sx q[2];
rz(1.9934838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6764073) q[1];
sx q[1];
rz(-2.8670951) q[1];
sx q[1];
rz(2.9390915) q[1];
x q[2];
rz(-1.3347606) q[3];
sx q[3];
rz(-2.1563362) q[3];
sx q[3];
rz(-1.3957617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8371381) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(-3.0755074) q[2];
rz(1.1262013) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(-0.57536212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5657625) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(2.3542812) q[0];
rz(-0.85834223) q[1];
sx q[1];
rz(-2.1635677) q[1];
sx q[1];
rz(-1.0816921) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6803241) q[0];
sx q[0];
rz(-0.96285614) q[0];
sx q[0];
rz(0.95198385) q[0];
x q[1];
rz(-2.2901847) q[2];
sx q[2];
rz(-2.5485196) q[2];
sx q[2];
rz(2.3874758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9138152) q[1];
sx q[1];
rz(-0.28474131) q[1];
sx q[1];
rz(-3.1193507) q[1];
x q[2];
rz(-2.2899707) q[3];
sx q[3];
rz(-1.1332261) q[3];
sx q[3];
rz(-2.9629666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9786238) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(2.2990189) q[2];
rz(-0.27635559) q[3];
sx q[3];
rz(-2.1601951) q[3];
sx q[3];
rz(1.3493376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.855298) q[0];
sx q[0];
rz(-1.3110302) q[0];
sx q[0];
rz(-1.1899765) q[0];
rz(1.2225993) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(-2.55866) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.615743) q[0];
sx q[0];
rz(-1.2229563) q[0];
sx q[0];
rz(-2.6418532) q[0];
rz(-pi) q[1];
rz(1.5491484) q[2];
sx q[2];
rz(-2.5679143) q[2];
sx q[2];
rz(1.4230185) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8376298) q[1];
sx q[1];
rz(-2.8899341) q[1];
sx q[1];
rz(2.8117517) q[1];
rz(-pi) q[2];
rz(-1.1357901) q[3];
sx q[3];
rz(-0.9851176) q[3];
sx q[3];
rz(0.034736573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0588093) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(2.1018551) q[2];
rz(2.5522363) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(2.1350433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650918) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(1.3366706) q[0];
rz(-0.22557766) q[1];
sx q[1];
rz(-1.0712737) q[1];
sx q[1];
rz(0.76990661) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192359) q[0];
sx q[0];
rz(-1.7261639) q[0];
sx q[0];
rz(-0.47433406) q[0];
x q[1];
rz(-0.628148) q[2];
sx q[2];
rz(-0.9315812) q[2];
sx q[2];
rz(-2.1232427) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4578145) q[1];
sx q[1];
rz(-1.646161) q[1];
sx q[1];
rz(0.1981258) q[1];
rz(-pi) q[2];
rz(0.1128686) q[3];
sx q[3];
rz(-2.5333571) q[3];
sx q[3];
rz(0.0054159482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5169107) q[2];
sx q[2];
rz(-1.7485488) q[2];
sx q[2];
rz(-2.1262271) q[2];
rz(0.20396248) q[3];
sx q[3];
rz(-1.5927619) q[3];
sx q[3];
rz(1.3913733) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90149752) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(-2.748306) q[0];
rz(0.45010629) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(-3.111305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6323551) q[0];
sx q[0];
rz(-1.8182753) q[0];
sx q[0];
rz(1.6028274) q[0];
rz(0.80937635) q[2];
sx q[2];
rz(-1.7580108) q[2];
sx q[2];
rz(2.9771114) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4612164) q[1];
sx q[1];
rz(-0.54619563) q[1];
sx q[1];
rz(-0.044574634) q[1];
x q[2];
rz(-2.6859517) q[3];
sx q[3];
rz(-2.8491088) q[3];
sx q[3];
rz(0.71575056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1071189) q[2];
sx q[2];
rz(-1.5566885) q[2];
sx q[2];
rz(-2.0195473) q[2];
rz(-1.0551039) q[3];
sx q[3];
rz(-0.3873581) q[3];
sx q[3];
rz(-1.25157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.17290641) q[0];
sx q[0];
rz(-0.98772573) q[0];
sx q[0];
rz(0.77447844) q[0];
rz(2.5075746) q[1];
sx q[1];
rz(-2.1114383) q[1];
sx q[1];
rz(-1.0672306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9241079) q[0];
sx q[0];
rz(-1.5536246) q[0];
sx q[0];
rz(-2.4033258) q[0];
rz(-pi) q[1];
rz(-1.7238317) q[2];
sx q[2];
rz(-0.31073292) q[2];
sx q[2];
rz(2.8513214) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6065403) q[1];
sx q[1];
rz(-2.7379704) q[1];
sx q[1];
rz(2.2500747) q[1];
x q[2];
rz(-1.0916294) q[3];
sx q[3];
rz(-1.7288343) q[3];
sx q[3];
rz(2.288186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60712236) q[2];
sx q[2];
rz(-3.0202713) q[2];
sx q[2];
rz(1.2657451) q[2];
rz(0.040146116) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-3.1025187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.48113111) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(-1.8817425) q[0];
rz(-1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(-0.0025509603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0028062971) q[0];
sx q[0];
rz(-1.696179) q[0];
sx q[0];
rz(3.0210397) q[0];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(1.8067183) q[1];
sx q[1];
rz(-2.2069227) q[1];
sx q[1];
rz(-1.0735372) q[1];
rz(-pi) q[2];
x q[2];
rz(0.02648377) q[3];
sx q[3];
rz(-2.0051574) q[3];
sx q[3];
rz(2.6514298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.586414) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(2.4671538) q[2];
rz(1.1174348) q[3];
sx q[3];
rz(-1.0267461) q[3];
sx q[3];
rz(-0.31291541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0565223) q[0];
sx q[0];
rz(-1.1310348) q[0];
sx q[0];
rz(-0.52666589) q[0];
rz(-2.7936735) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(2.6013026) q[2];
sx q[2];
rz(-2.7439596) q[2];
sx q[2];
rz(-0.31184218) q[2];
rz(-1.6913577) q[3];
sx q[3];
rz(-2.9942054) q[3];
sx q[3];
rz(-1.2248539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
