OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6586128) q[0];
sx q[0];
rz(-0.40402544) q[0];
sx q[0];
rz(-2.7513096) q[0];
rz(3.1080988) q[1];
sx q[1];
rz(-0.53116763) q[1];
sx q[1];
rz(2.9386428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4056978) q[0];
sx q[0];
rz(-0.66758388) q[0];
sx q[0];
rz(-2.1632458) q[0];
rz(-pi) q[1];
rz(-3.13596) q[2];
sx q[2];
rz(-1.4818824) q[2];
sx q[2];
rz(-2.2064457) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81357274) q[1];
sx q[1];
rz(-1.7576694) q[1];
sx q[1];
rz(0.55208167) q[1];
x q[2];
rz(-3.1139939) q[3];
sx q[3];
rz(-2.2970133) q[3];
sx q[3];
rz(0.82311326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94886327) q[2];
sx q[2];
rz(-2.296083) q[2];
sx q[2];
rz(-1.753099) q[2];
rz(0.071831547) q[3];
sx q[3];
rz(-2.6303232) q[3];
sx q[3];
rz(0.82740074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097505957) q[0];
sx q[0];
rz(-2.9768017) q[0];
sx q[0];
rz(0.49758115) q[0];
rz(-1.3961821) q[1];
sx q[1];
rz(-2.1131682) q[1];
sx q[1];
rz(-0.57026774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1193863) q[0];
sx q[0];
rz(-2.019068) q[0];
sx q[0];
rz(-1.717685) q[0];
rz(-pi) q[1];
rz(1.360434) q[2];
sx q[2];
rz(-1.3704925) q[2];
sx q[2];
rz(0.64678538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1062938) q[1];
sx q[1];
rz(-2.4199644) q[1];
sx q[1];
rz(-0.053573805) q[1];
rz(-pi) q[2];
rz(-0.87941283) q[3];
sx q[3];
rz(-2.0447404) q[3];
sx q[3];
rz(0.64521257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1121062) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(0.71375978) q[2];
rz(-1.3826238) q[3];
sx q[3];
rz(-0.52240038) q[3];
sx q[3];
rz(2.6944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5722028) q[0];
sx q[0];
rz(-2.1694006) q[0];
sx q[0];
rz(-2.8045281) q[0];
rz(0.49346787) q[1];
sx q[1];
rz(-2.4512873) q[1];
sx q[1];
rz(0.92672551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2703646) q[0];
sx q[0];
rz(-2.7634826) q[0];
sx q[0];
rz(-2.3367995) q[0];
rz(2.3078467) q[2];
sx q[2];
rz(-1.0930702) q[2];
sx q[2];
rz(-1.1857978) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.265343) q[1];
sx q[1];
rz(-2.690965) q[1];
sx q[1];
rz(-1.9540811) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51982359) q[3];
sx q[3];
rz(-0.9614203) q[3];
sx q[3];
rz(1.6345616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9530764) q[2];
sx q[2];
rz(-2.2761554) q[2];
sx q[2];
rz(2.4600929) q[2];
rz(-1.6279047) q[3];
sx q[3];
rz(-1.3368006) q[3];
sx q[3];
rz(2.9708235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3156768) q[0];
sx q[0];
rz(-3.0382394) q[0];
sx q[0];
rz(-2.090825) q[0];
rz(0.61327618) q[1];
sx q[1];
rz(-0.79137099) q[1];
sx q[1];
rz(-1.223986) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67658778) q[0];
sx q[0];
rz(-0.95116827) q[0];
sx q[0];
rz(-1.1628435) q[0];
rz(2.1992219) q[2];
sx q[2];
rz(-2.0258459) q[2];
sx q[2];
rz(-1.1466743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.115009) q[1];
sx q[1];
rz(-1.5582651) q[1];
sx q[1];
rz(2.5471155) q[1];
rz(-pi) q[2];
rz(0.23777407) q[3];
sx q[3];
rz(-0.60543767) q[3];
sx q[3];
rz(2.0914222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0918538) q[2];
sx q[2];
rz(-0.78511304) q[2];
sx q[2];
rz(-2.4741057) q[2];
rz(-1.5664172) q[3];
sx q[3];
rz(-0.19078855) q[3];
sx q[3];
rz(-0.13458399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5176373) q[0];
sx q[0];
rz(-2.1006382) q[0];
sx q[0];
rz(-0.55476302) q[0];
rz(-1.293921) q[1];
sx q[1];
rz(-2.7298253) q[1];
sx q[1];
rz(-2.256934) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8644476) q[0];
sx q[0];
rz(-0.53098035) q[0];
sx q[0];
rz(1.619521) q[0];
x q[1];
rz(-2.4403202) q[2];
sx q[2];
rz(-1.9375357) q[2];
sx q[2];
rz(-1.203095) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.64120871) q[1];
sx q[1];
rz(-0.51189089) q[1];
sx q[1];
rz(-1.1275395) q[1];
rz(-pi) q[2];
rz(1.8507666) q[3];
sx q[3];
rz(-1.1445224) q[3];
sx q[3];
rz(0.72457641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.027792949) q[2];
sx q[2];
rz(-2.5774559) q[2];
sx q[2];
rz(-1.0998868) q[2];
rz(0.22282985) q[3];
sx q[3];
rz(-1.9375216) q[3];
sx q[3];
rz(-3.0075464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851538) q[0];
sx q[0];
rz(-2.8414861) q[0];
sx q[0];
rz(-1.9578178) q[0];
rz(-1.3166332) q[1];
sx q[1];
rz(-0.96375179) q[1];
sx q[1];
rz(1.0354985) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1610704) q[0];
sx q[0];
rz(-0.59605205) q[0];
sx q[0];
rz(-0.97575326) q[0];
rz(-pi) q[1];
rz(2.1450317) q[2];
sx q[2];
rz(-2.445526) q[2];
sx q[2];
rz(2.0526759) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.07044) q[1];
sx q[1];
rz(-2.8912918) q[1];
sx q[1];
rz(0.76257272) q[1];
rz(3.0665056) q[3];
sx q[3];
rz(-2.3737566) q[3];
sx q[3];
rz(-1.7164643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3731132) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(2.7601472) q[2];
rz(-1.2162195) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(-0.60429627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4878047) q[0];
sx q[0];
rz(-0.31823802) q[0];
sx q[0];
rz(2.8741264) q[0];
rz(-1.5221315) q[1];
sx q[1];
rz(-2.5289502) q[1];
sx q[1];
rz(2.6681275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8657314) q[0];
sx q[0];
rz(-1.9124219) q[0];
sx q[0];
rz(0.39535687) q[0];
rz(-pi) q[1];
rz(-1.7554088) q[2];
sx q[2];
rz(-1.0981596) q[2];
sx q[2];
rz(-0.71834823) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5468041) q[1];
sx q[1];
rz(-1.4546397) q[1];
sx q[1];
rz(-0.22211566) q[1];
rz(-pi) q[2];
rz(-2.7678185) q[3];
sx q[3];
rz(-1.7229986) q[3];
sx q[3];
rz(-2.0836326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21753103) q[2];
sx q[2];
rz(-1.4950098) q[2];
sx q[2];
rz(-1.3727429) q[2];
rz(-2.8352906) q[3];
sx q[3];
rz(-2.114571) q[3];
sx q[3];
rz(-0.33941227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.821625) q[0];
sx q[0];
rz(-2.8447633) q[0];
sx q[0];
rz(-0.4739652) q[0];
rz(-2.3560246) q[1];
sx q[1];
rz(-2.5626917) q[1];
sx q[1];
rz(-2.2483291) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1635904) q[0];
sx q[0];
rz(-1.8390391) q[0];
sx q[0];
rz(0.017649529) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39126663) q[2];
sx q[2];
rz(-0.80384582) q[2];
sx q[2];
rz(0.93144691) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.230669) q[1];
sx q[1];
rz(-1.5381457) q[1];
sx q[1];
rz(2.9327675) q[1];
rz(3.0419218) q[3];
sx q[3];
rz(-1.1766889) q[3];
sx q[3];
rz(-1.3083713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6713509) q[2];
sx q[2];
rz(-1.525815) q[2];
sx q[2];
rz(-0.5980171) q[2];
rz(0.20448576) q[3];
sx q[3];
rz(-0.11071591) q[3];
sx q[3];
rz(2.428875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04190271) q[0];
sx q[0];
rz(-2.1429017) q[0];
sx q[0];
rz(-2.4424851) q[0];
rz(2.7514669) q[1];
sx q[1];
rz(-2.4603619) q[1];
sx q[1];
rz(2.1652538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028653305) q[0];
sx q[0];
rz(-1.8277455) q[0];
sx q[0];
rz(-2.9963958) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62841793) q[2];
sx q[2];
rz(-2.6172514) q[2];
sx q[2];
rz(0.4195329) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0864549) q[1];
sx q[1];
rz(-1.6443776) q[1];
sx q[1];
rz(-1.0615968) q[1];
rz(-pi) q[2];
rz(2.7224119) q[3];
sx q[3];
rz(-2.5435735) q[3];
sx q[3];
rz(0.49367762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98328996) q[2];
sx q[2];
rz(-0.96902865) q[2];
sx q[2];
rz(0.14652531) q[2];
rz(-2.8807785) q[3];
sx q[3];
rz(-2.0050037) q[3];
sx q[3];
rz(-0.28723106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89676595) q[0];
sx q[0];
rz(-0.34938669) q[0];
sx q[0];
rz(0.31325999) q[0];
rz(0.80097711) q[1];
sx q[1];
rz(-1.6311092) q[1];
sx q[1];
rz(-0.40447485) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0718719) q[0];
sx q[0];
rz(-1.7554465) q[0];
sx q[0];
rz(1.429256) q[0];
x q[1];
rz(1.7687665) q[2];
sx q[2];
rz(-2.1083197) q[2];
sx q[2];
rz(-0.24453577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69533379) q[1];
sx q[1];
rz(-2.0768169) q[1];
sx q[1];
rz(2.395588) q[1];
rz(1.9451109) q[3];
sx q[3];
rz(-0.35983837) q[3];
sx q[3];
rz(-1.8051749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6910088) q[2];
sx q[2];
rz(-2.6477224) q[2];
sx q[2];
rz(2.3502926) q[2];
rz(-0.50296909) q[3];
sx q[3];
rz(-2.0937604) q[3];
sx q[3];
rz(-2.3414229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41225152) q[0];
sx q[0];
rz(-1.3164192) q[0];
sx q[0];
rz(-1.3715716) q[0];
rz(0.62660632) q[1];
sx q[1];
rz(-1.659844) q[1];
sx q[1];
rz(-1.0214092) q[1];
rz(1.4087497) q[2];
sx q[2];
rz(-1.9377943) q[2];
sx q[2];
rz(1.6429907) q[2];
rz(1.5833686) q[3];
sx q[3];
rz(-1.8042121) q[3];
sx q[3];
rz(1.369759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
