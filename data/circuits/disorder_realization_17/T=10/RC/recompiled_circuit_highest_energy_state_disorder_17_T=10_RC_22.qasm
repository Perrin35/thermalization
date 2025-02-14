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
rz(-1.2920657) q[0];
sx q[0];
rz(-2.3176365) q[0];
sx q[0];
rz(-2.2141461) q[0];
rz(2.0916341) q[1];
sx q[1];
rz(-0.97133049) q[1];
sx q[1];
rz(-0.8144905) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2789291) q[0];
sx q[0];
rz(-1.188108) q[0];
sx q[0];
rz(1.0839868) q[0];
rz(-pi) q[1];
rz(0.94169803) q[2];
sx q[2];
rz(-0.88028958) q[2];
sx q[2];
rz(2.5494573) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3315288) q[1];
sx q[1];
rz(-2.0671856) q[1];
sx q[1];
rz(1.0277102) q[1];
rz(-0.82858927) q[3];
sx q[3];
rz(-2.0999206) q[3];
sx q[3];
rz(0.55381004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77259511) q[2];
sx q[2];
rz(-1.3694222) q[2];
sx q[2];
rz(0.43130809) q[2];
rz(1.6469693) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(0.78540426) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876708) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(1.824463) q[0];
rz(-2.092284) q[1];
sx q[1];
rz(-2.6056555) q[1];
sx q[1];
rz(-2.5228693) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0174219) q[0];
sx q[0];
rz(-0.96786495) q[0];
sx q[0];
rz(-1.746491) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0929865) q[2];
sx q[2];
rz(-2.3903762) q[2];
sx q[2];
rz(-0.71556811) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2545027) q[1];
sx q[1];
rz(-1.2225593) q[1];
sx q[1];
rz(-0.13910267) q[1];
rz(-1.9369164) q[3];
sx q[3];
rz(-2.0799326) q[3];
sx q[3];
rz(-2.3222271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.416136) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(2.6497192) q[2];
rz(-1.5791996) q[3];
sx q[3];
rz(-1.4926566) q[3];
sx q[3];
rz(-0.57945848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270776) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(0.28011093) q[0];
rz(-3.0666871) q[1];
sx q[1];
rz(-2.7066051) q[1];
sx q[1];
rz(1.7710955) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89557392) q[0];
sx q[0];
rz(-1.361712) q[0];
sx q[0];
rz(0.81935291) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5016194) q[2];
sx q[2];
rz(-1.0037862) q[2];
sx q[2];
rz(0.28653539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6447923) q[1];
sx q[1];
rz(-0.62207568) q[1];
sx q[1];
rz(0.46228564) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5896443) q[3];
sx q[3];
rz(-1.7030088) q[3];
sx q[3];
rz(2.978505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8568153) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(-0.80186296) q[2];
rz(0.91184584) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(-0.97889939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80482471) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(-2.5208852) q[0];
rz(-2.8630818) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(1.0034358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1708978) q[0];
sx q[0];
rz(-1.2549632) q[0];
sx q[0];
rz(1.9307846) q[0];
rz(1.150028) q[2];
sx q[2];
rz(-1.506797) q[2];
sx q[2];
rz(-0.45129946) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.089503376) q[1];
sx q[1];
rz(-1.6253396) q[1];
sx q[1];
rz(2.8724345) q[1];
rz(-2.8025556) q[3];
sx q[3];
rz(-0.62612306) q[3];
sx q[3];
rz(2.1563185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8371381) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(-3.0755074) q[2];
rz(1.1262013) q[3];
sx q[3];
rz(-0.76032138) q[3];
sx q[3];
rz(0.57536212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5657625) q[0];
sx q[0];
rz(-3.059721) q[0];
sx q[0];
rz(0.78731147) q[0];
rz(0.85834223) q[1];
sx q[1];
rz(-0.97802496) q[1];
sx q[1];
rz(2.0599005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72321945) q[0];
sx q[0];
rz(-1.0745418) q[0];
sx q[0];
rz(-0.70696522) q[0];
rz(-pi) q[1];
rz(1.101564) q[2];
sx q[2];
rz(-1.9479556) q[2];
sx q[2];
rz(-2.9531329) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.25095233) q[1];
sx q[1];
rz(-1.8554652) q[1];
sx q[1];
rz(-1.5773058) q[1];
x q[2];
rz(2.2899707) q[3];
sx q[3];
rz(-1.1332261) q[3];
sx q[3];
rz(2.9629666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9786238) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(0.84257379) q[2];
rz(0.27635559) q[3];
sx q[3];
rz(-0.98139757) q[3];
sx q[3];
rz(1.3493376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2862947) q[0];
sx q[0];
rz(-1.3110302) q[0];
sx q[0];
rz(-1.1899765) q[0];
rz(1.2225993) q[1];
sx q[1];
rz(-1.2346376) q[1];
sx q[1];
rz(2.55866) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5385206) q[0];
sx q[0];
rz(-2.5412509) q[0];
sx q[0];
rz(-2.4938514) q[0];
rz(0.013986258) q[2];
sx q[2];
rz(-0.99726935) q[2];
sx q[2];
rz(-1.397246) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5869904) q[1];
sx q[1];
rz(-1.6515367) q[1];
sx q[1];
rz(2.9029773) q[1];
rz(-pi) q[2];
rz(1.1357901) q[3];
sx q[3];
rz(-0.9851176) q[3];
sx q[3];
rz(-0.034736573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0588093) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(2.1018551) q[2];
rz(-2.5522363) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(1.0065494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650918) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(1.3366706) q[0];
rz(-2.916015) q[1];
sx q[1];
rz(-2.070319) q[1];
sx q[1];
rz(-2.371686) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63084376) q[0];
sx q[0];
rz(-2.0389557) q[0];
sx q[0];
rz(-1.3965142) q[0];
x q[1];
rz(0.82775292) q[2];
sx q[2];
rz(-1.079664) q[2];
sx q[2];
rz(0.14358768) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.47178956) q[1];
sx q[1];
rz(-0.21179971) q[1];
sx q[1];
rz(0.36630156) q[1];
rz(-0.60524551) q[3];
sx q[3];
rz(-1.5063933) q[3];
sx q[3];
rz(1.6689672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(1.0153655) q[2];
rz(-2.9376302) q[3];
sx q[3];
rz(-1.5927619) q[3];
sx q[3];
rz(-1.7502194) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90149752) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(2.748306) q[0];
rz(0.45010629) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(-3.111305) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6393042) q[0];
sx q[0];
rz(-2.8920916) q[0];
sx q[0];
rz(-3.0155103) q[0];
x q[1];
rz(0.80937635) q[2];
sx q[2];
rz(-1.7580108) q[2];
sx q[2];
rz(-0.16448122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4612164) q[1];
sx q[1];
rz(-0.54619563) q[1];
sx q[1];
rz(0.044574634) q[1];
rz(-1.7025331) q[3];
sx q[3];
rz(-1.3088969) q[3];
sx q[3];
rz(1.9528263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0344737) q[2];
sx q[2];
rz(-1.5849042) q[2];
sx q[2];
rz(2.0195473) q[2];
rz(1.0551039) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(1.8900227) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17290641) q[0];
sx q[0];
rz(-0.98772573) q[0];
sx q[0];
rz(0.77447844) q[0];
rz(0.6340181) q[1];
sx q[1];
rz(-1.0301544) q[1];
sx q[1];
rz(2.074362) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8071496) q[0];
sx q[0];
rz(-2.4031638) q[0];
sx q[0];
rz(-0.025511857) q[0];
rz(-pi) q[1];
rz(3.0926782) q[2];
sx q[2];
rz(-1.2638143) q[2];
sx q[2];
rz(0.129667) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6744819) q[1];
sx q[1];
rz(-1.8201105) q[1];
sx q[1];
rz(-1.8915908) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9036071) q[3];
sx q[3];
rz(-0.5026256) q[3];
sx q[3];
rz(0.42325936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5344703) q[2];
sx q[2];
rz(-3.0202713) q[2];
sx q[2];
rz(-1.2657451) q[2];
rz(0.040146116) q[3];
sx q[3];
rz(-2.7669192) q[3];
sx q[3];
rz(-0.039073959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6604615) q[0];
sx q[0];
rz(-0.88812319) q[0];
sx q[0];
rz(1.2598502) q[0];
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
rz(3.1387864) q[0];
sx q[0];
rz(-1.4454137) q[0];
sx q[0];
rz(0.12055293) q[0];
rz(-3.1319982) q[2];
sx q[2];
rz(-1.506732) q[2];
sx q[2];
rz(1.715915) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8067183) q[1];
sx q[1];
rz(-2.2069227) q[1];
sx q[1];
rz(2.0680554) q[1];
x q[2];
rz(1.1363013) q[3];
sx q[3];
rz(-1.5467724) q[3];
sx q[3];
rz(2.0498118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.586414) q[2];
sx q[2];
rz(-1.5191583) q[2];
sx q[2];
rz(0.67443887) q[2];
rz(2.0241578) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(-0.31291541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085070327) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(-0.3479192) q[1];
sx q[1];
rz(-1.2384474) q[1];
sx q[1];
rz(1.7927982) q[1];
rz(-1.3580218) q[2];
sx q[2];
rz(-1.909303) q[2];
sx q[2];
rz(2.252966) q[2];
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
