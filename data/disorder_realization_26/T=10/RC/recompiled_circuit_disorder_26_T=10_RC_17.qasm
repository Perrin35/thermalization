OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(-1.4533071) q[0];
sx q[0];
rz(0.31153554) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(4.9649927) q[1];
sx q[1];
rz(8.8658219) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6214949) q[0];
sx q[0];
rz(-1.9858452) q[0];
sx q[0];
rz(0.15226224) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55732255) q[2];
sx q[2];
rz(-1.5814591) q[2];
sx q[2];
rz(2.9023841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4999057) q[1];
sx q[1];
rz(-0.99398621) q[1];
sx q[1];
rz(0.66024248) q[1];
x q[2];
rz(-1.8564838) q[3];
sx q[3];
rz(-1.9938855) q[3];
sx q[3];
rz(2.3959514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7493593) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(-0.63670811) q[2];
rz(-2.2926245) q[3];
sx q[3];
rz(-0.62148062) q[3];
sx q[3];
rz(-2.9076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37671509) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(-0.15287457) q[0];
rz(2.3846467) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(-0.98639948) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2043641) q[0];
sx q[0];
rz(-1.5110656) q[0];
sx q[0];
rz(-1.5304969) q[0];
rz(-pi) q[1];
rz(1.0436922) q[2];
sx q[2];
rz(-2.2615848) q[2];
sx q[2];
rz(1.3259128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4437372) q[1];
sx q[1];
rz(-1.9471696) q[1];
sx q[1];
rz(1.8259551) q[1];
rz(-pi) q[2];
rz(0.13568474) q[3];
sx q[3];
rz(-1.8968582) q[3];
sx q[3];
rz(-0.56860926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5622921) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(2.3584649) q[2];
rz(3.1230208) q[3];
sx q[3];
rz(-1.5037856) q[3];
sx q[3];
rz(2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0531533) q[0];
sx q[0];
rz(-0.29215559) q[0];
sx q[0];
rz(0.96167481) q[0];
rz(0.36034521) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(-0.12869421) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0592244) q[0];
sx q[0];
rz(-1.6186065) q[0];
sx q[0];
rz(1.6533018) q[0];
rz(-pi) q[1];
rz(-2.6016597) q[2];
sx q[2];
rz(-2.2579102) q[2];
sx q[2];
rz(-0.5772669) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5161908) q[1];
sx q[1];
rz(-2.0149724) q[1];
sx q[1];
rz(-1.8810012) q[1];
x q[2];
rz(-0.13462984) q[3];
sx q[3];
rz(-0.80727808) q[3];
sx q[3];
rz(1.5869629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0814357) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(1.8998247) q[2];
rz(-2.5545819) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(-2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509165) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(-2.0571016) q[0];
rz(-1.4831316) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(-0.09253563) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4400892) q[0];
sx q[0];
rz(-1.2883696) q[0];
sx q[0];
rz(-2.861172) q[0];
rz(-pi) q[1];
rz(0.34543583) q[2];
sx q[2];
rz(-1.1159117) q[2];
sx q[2];
rz(-2.5330184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.98574084) q[1];
sx q[1];
rz(-2.7936613) q[1];
sx q[1];
rz(0.17133979) q[1];
rz(-pi) q[2];
rz(-1.7792286) q[3];
sx q[3];
rz(-0.66550335) q[3];
sx q[3];
rz(0.044737577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3080421) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(-0.33205024) q[2];
rz(-1.0559233) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8191391) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(2.9300368) q[0];
rz(1.8353204) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(0.64770118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6015198) q[0];
sx q[0];
rz(-1.411502) q[0];
sx q[0];
rz(1.6838616) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7251882) q[2];
sx q[2];
rz(-1.9152181) q[2];
sx q[2];
rz(-1.7248578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.065121) q[1];
sx q[1];
rz(-1.3741125) q[1];
sx q[1];
rz(-2.0842488) q[1];
rz(1.5557489) q[3];
sx q[3];
rz(-1.589937) q[3];
sx q[3];
rz(2.0862938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5806879) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(0.69331759) q[2];
rz(2.4723315) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(-0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(-1.2325226) q[0];
rz(-2.0690074) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(0.17428621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78290735) q[0];
sx q[0];
rz(-1.2619962) q[0];
sx q[0];
rz(2.3147644) q[0];
x q[1];
rz(0.86485483) q[2];
sx q[2];
rz(-1.0360498) q[2];
sx q[2];
rz(-2.6851482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.41387687) q[1];
sx q[1];
rz(-1.9503647) q[1];
sx q[1];
rz(2.3351401) q[1];
x q[2];
rz(-1.2234736) q[3];
sx q[3];
rz(-2.8109549) q[3];
sx q[3];
rz(-0.91352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8217414) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(-2.9439587) q[2];
rz(-0.28891426) q[3];
sx q[3];
rz(-2.2646326) q[3];
sx q[3];
rz(1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.063868) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(2.9329964) q[0];
rz(0.96616191) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(-1.5055515) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063232139) q[0];
sx q[0];
rz(-0.97023836) q[0];
sx q[0];
rz(-2.2230704) q[0];
rz(-2.7111972) q[2];
sx q[2];
rz(-1.6039404) q[2];
sx q[2];
rz(2.7660649) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.437285) q[1];
sx q[1];
rz(-1.1204801) q[1];
sx q[1];
rz(0.98547658) q[1];
x q[2];
rz(0.286245) q[3];
sx q[3];
rz(-1.6702594) q[3];
sx q[3];
rz(-3.0486097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.85764) q[2];
sx q[2];
rz(-2.3985034) q[2];
sx q[2];
rz(-0.097578438) q[2];
rz(-1.3939259) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(2.424749) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39032787) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(0.51399291) q[0];
rz(-0.12318525) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(0.93200144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1876353) q[0];
sx q[0];
rz(-1.4353416) q[0];
sx q[0];
rz(2.016469) q[0];
rz(-0.26884218) q[2];
sx q[2];
rz(-2.1890867) q[2];
sx q[2];
rz(2.6593069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92703687) q[1];
sx q[1];
rz(-2.3408457) q[1];
sx q[1];
rz(-2.3826249) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85429116) q[3];
sx q[3];
rz(-1.4201418) q[3];
sx q[3];
rz(0.52779576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0294068) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(2.712148) q[2];
rz(-1.2094234) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76686239) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(-1.3388348) q[0];
rz(0.70612899) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(2.1910117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50981748) q[0];
sx q[0];
rz(-2.0777367) q[0];
sx q[0];
rz(-0.99960534) q[0];
rz(3.0663475) q[2];
sx q[2];
rz(-1.9855472) q[2];
sx q[2];
rz(-3.0441949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95987684) q[1];
sx q[1];
rz(-1.0392337) q[1];
sx q[1];
rz(-2.9243484) q[1];
x q[2];
rz(-0.25456984) q[3];
sx q[3];
rz(-1.1708784) q[3];
sx q[3];
rz(-3.1018156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(2.2144923) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1442239) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(-2.9123059) q[0];
rz(-0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(0.71892175) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7867891) q[0];
sx q[0];
rz(-2.1956586) q[0];
sx q[0];
rz(3.025269) q[0];
rz(1.3942137) q[2];
sx q[2];
rz(-0.59497661) q[2];
sx q[2];
rz(-0.31756155) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.194866) q[1];
sx q[1];
rz(-0.39992878) q[1];
sx q[1];
rz(-2.7547794) q[1];
rz(-pi) q[2];
rz(-2.1925681) q[3];
sx q[3];
rz(-1.1488631) q[3];
sx q[3];
rz(-0.66872795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7583313) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(0.74679217) q[2];
rz(-0.87219277) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(2.1993568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.223021) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(0.37721286) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(2.626426) q[2];
sx q[2];
rz(-0.58044051) q[2];
sx q[2];
rz(2.6945111) q[2];
rz(-1.4217581) q[3];
sx q[3];
rz(-2.2825713) q[3];
sx q[3];
rz(1.5406516) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
