OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(2.252993) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(-1.5421901) q[1];
sx q[1];
rz(1.3808274) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5262573) q[0];
sx q[0];
rz(-1.3852296) q[0];
sx q[0];
rz(-0.19677563) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1031453) q[2];
sx q[2];
rz(-0.84538904) q[2];
sx q[2];
rz(-0.40276819) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0128768) q[1];
sx q[1];
rz(-0.33564645) q[1];
sx q[1];
rz(3.1124093) q[1];
rz(-pi) q[2];
rz(-2.9167261) q[3];
sx q[3];
rz(-1.9266085) q[3];
sx q[3];
rz(2.0840621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8218653) q[2];
sx q[2];
rz(-1.3250019) q[2];
sx q[2];
rz(0.74903178) q[2];
rz(2.9876515) q[3];
sx q[3];
rz(-2.1059683) q[3];
sx q[3];
rz(3.0416987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.446949) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(1.9248167) q[0];
rz(-1.0617537) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(-2.7144576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0674853) q[0];
sx q[0];
rz(-2.4541313) q[0];
sx q[0];
rz(-0.56337728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92556503) q[2];
sx q[2];
rz(-1.7665909) q[2];
sx q[2];
rz(-2.2044971) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6500351) q[1];
sx q[1];
rz(-1.2833529) q[1];
sx q[1];
rz(-0.14658714) q[1];
rz(3.0892761) q[3];
sx q[3];
rz(-1.3671095) q[3];
sx q[3];
rz(1.2869204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.079166807) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(1.3986577) q[2];
rz(0.22377293) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1200714) q[0];
sx q[0];
rz(-1.7912309) q[0];
sx q[0];
rz(3.0960826) q[0];
rz(-1.1687763) q[1];
sx q[1];
rz(-1.903542) q[1];
sx q[1];
rz(-2.5659836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94852266) q[0];
sx q[0];
rz(-0.93174495) q[0];
sx q[0];
rz(1.5847413) q[0];
rz(3.124442) q[2];
sx q[2];
rz(-0.91020012) q[2];
sx q[2];
rz(-1.2397546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82390416) q[1];
sx q[1];
rz(-2.1145193) q[1];
sx q[1];
rz(-2.8871817) q[1];
rz(-pi) q[2];
rz(1.5381673) q[3];
sx q[3];
rz(-0.65308076) q[3];
sx q[3];
rz(-3.0663222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0098972926) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(0.95727813) q[2];
rz(2.5668528) q[3];
sx q[3];
rz(-1.6114019) q[3];
sx q[3];
rz(-1.3055698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15605536) q[0];
sx q[0];
rz(-0.95273459) q[0];
sx q[0];
rz(-1.2611457) q[0];
rz(3.0592697) q[1];
sx q[1];
rz(-2.0539093) q[1];
sx q[1];
rz(-1.7877158) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1607504) q[0];
sx q[0];
rz(-0.35355648) q[0];
sx q[0];
rz(-0.73676957) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8026695) q[2];
sx q[2];
rz(-2.6312575) q[2];
sx q[2];
rz(-1.5938544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30198797) q[1];
sx q[1];
rz(-1.0763554) q[1];
sx q[1];
rz(1.5093263) q[1];
rz(-pi) q[2];
rz(2.6526101) q[3];
sx q[3];
rz(-0.95150286) q[3];
sx q[3];
rz(1.796738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42133078) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(-1.2565695) q[2];
rz(0.71150696) q[3];
sx q[3];
rz(-1.3308176) q[3];
sx q[3];
rz(-2.8594657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.310815) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(2.7984483) q[0];
rz(-0.061773069) q[1];
sx q[1];
rz(-2.1715178) q[1];
sx q[1];
rz(1.4168581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23487906) q[0];
sx q[0];
rz(-0.32787927) q[0];
sx q[0];
rz(-1.0176786) q[0];
x q[1];
rz(0.70059641) q[2];
sx q[2];
rz(-1.3178409) q[2];
sx q[2];
rz(-2.7299936) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2572391) q[1];
sx q[1];
rz(-1.7063008) q[1];
sx q[1];
rz(2.0550904) q[1];
rz(2.9067578) q[3];
sx q[3];
rz(-1.1560625) q[3];
sx q[3];
rz(-2.107634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5477649) q[2];
sx q[2];
rz(-1.4676899) q[2];
sx q[2];
rz(-1.0859547) q[2];
rz(-0.91935277) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(2.5277188) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74349657) q[0];
sx q[0];
rz(-1.9323823) q[0];
sx q[0];
rz(-2.3486775) q[0];
rz(-0.74053699) q[1];
sx q[1];
rz(-0.99895993) q[1];
sx q[1];
rz(2.3103796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8342469) q[0];
sx q[0];
rz(-1.9458658) q[0];
sx q[0];
rz(-2.9177279) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0074499091) q[2];
sx q[2];
rz(-1.8478571) q[2];
sx q[2];
rz(-1.3598833) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38922849) q[1];
sx q[1];
rz(-0.78654754) q[1];
sx q[1];
rz(-2.0081372) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5300203) q[3];
sx q[3];
rz(-0.47104657) q[3];
sx q[3];
rz(-0.99179964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.071216019) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(-2.0224723) q[2];
rz(1.8708771) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(-1.7601815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374461) q[0];
sx q[0];
rz(-0.16682145) q[0];
sx q[0];
rz(-1.5995837) q[0];
rz(-2.1265325) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(-0.17280811) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8625582) q[0];
sx q[0];
rz(-2.1365215) q[0];
sx q[0];
rz(-0.59086694) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4577499) q[2];
sx q[2];
rz(-1.2798066) q[2];
sx q[2];
rz(0.061372193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.648293) q[1];
sx q[1];
rz(-1.9956335) q[1];
sx q[1];
rz(0.27314911) q[1];
x q[2];
rz(0.5676078) q[3];
sx q[3];
rz(-2.2877573) q[3];
sx q[3];
rz(0.53138083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66199866) q[2];
sx q[2];
rz(-2.2981503) q[2];
sx q[2];
rz(-2.1638828) q[2];
rz(-2.5665723) q[3];
sx q[3];
rz(-1.9366555) q[3];
sx q[3];
rz(-1.9112816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32790312) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(-0.015722474) q[0];
rz(-1.188259) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(-1.9421008) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5003842) q[0];
sx q[0];
rz(-0.50725021) q[0];
sx q[0];
rz(1.070606) q[0];
x q[1];
rz(-2.0616777) q[2];
sx q[2];
rz(-1.5245617) q[2];
sx q[2];
rz(-0.85109988) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6576865) q[1];
sx q[1];
rz(-2.4281341) q[1];
sx q[1];
rz(-0.76303996) q[1];
rz(2.0301129) q[3];
sx q[3];
rz(-2.5190341) q[3];
sx q[3];
rz(-2.7186269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99379313) q[2];
sx q[2];
rz(-1.9027998) q[2];
sx q[2];
rz(1.3851059) q[2];
rz(0.6238474) q[3];
sx q[3];
rz(-1.8892989) q[3];
sx q[3];
rz(2.0417716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1856336) q[0];
sx q[0];
rz(-0.82056844) q[0];
sx q[0];
rz(-0.69325915) q[0];
rz(0.26501003) q[1];
sx q[1];
rz(-2.3131504) q[1];
sx q[1];
rz(1.5230491) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023848195) q[0];
sx q[0];
rz(-1.5818412) q[0];
sx q[0];
rz(2.3680229) q[0];
rz(-pi) q[1];
rz(-2.393894) q[2];
sx q[2];
rz(-0.77538632) q[2];
sx q[2];
rz(-0.86624399) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56227389) q[1];
sx q[1];
rz(-2.8011311) q[1];
sx q[1];
rz(3.0698983) q[1];
rz(-pi) q[2];
rz(-1.6783488) q[3];
sx q[3];
rz(-0.58379025) q[3];
sx q[3];
rz(-1.3046169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30274621) q[2];
sx q[2];
rz(-1.5134209) q[2];
sx q[2];
rz(1.4576853) q[2];
rz(0.95528209) q[3];
sx q[3];
rz(-2.6421319) q[3];
sx q[3];
rz(2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38462287) q[0];
sx q[0];
rz(-0.56461016) q[0];
sx q[0];
rz(0.48267522) q[0];
rz(-0.92974281) q[1];
sx q[1];
rz(-2.1371806) q[1];
sx q[1];
rz(0.39628705) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7571054) q[0];
sx q[0];
rz(-1.9845909) q[0];
sx q[0];
rz(1.1767469) q[0];
rz(2.9626289) q[2];
sx q[2];
rz(-1.71993) q[2];
sx q[2];
rz(2.3863132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0490108) q[1];
sx q[1];
rz(-1.7515469) q[1];
sx q[1];
rz(2.3011219) q[1];
rz(-pi) q[2];
rz(-1.2108874) q[3];
sx q[3];
rz(-1.9511838) q[3];
sx q[3];
rz(1.4402657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79260176) q[2];
sx q[2];
rz(-1.7020117) q[2];
sx q[2];
rz(-0.3271884) q[2];
rz(2.3606825) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(-1.6646632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0384211) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(-0.36956638) q[1];
sx q[1];
rz(-0.8820487) q[1];
sx q[1];
rz(-0.65912156) q[1];
rz(-2.4482881) q[2];
sx q[2];
rz(-1.605576) q[2];
sx q[2];
rz(-1.3645542) q[2];
rz(-0.9624858) q[3];
sx q[3];
rz(-1.8909834) q[3];
sx q[3];
rz(-3.0743619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
