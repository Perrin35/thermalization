OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31792274) q[0];
sx q[0];
rz(-1.5225141) q[0];
sx q[0];
rz(-0.063152753) q[0];
rz(-1.3178648) q[1];
sx q[1];
rz(3.5409034) q[1];
sx q[1];
rz(10.200722) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0948333) q[0];
sx q[0];
rz(-1.6164808) q[0];
sx q[0];
rz(1.5976357) q[0];
rz(-pi) q[1];
rz(-2.6444998) q[2];
sx q[2];
rz(-0.99537206) q[2];
sx q[2];
rz(-2.6445553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5835442) q[1];
sx q[1];
rz(-0.87178245) q[1];
sx q[1];
rz(-2.9093161) q[1];
rz(-pi) q[2];
rz(1.8100061) q[3];
sx q[3];
rz(-2.0612217) q[3];
sx q[3];
rz(-1.3077206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45106384) q[2];
sx q[2];
rz(-1.7454001) q[2];
sx q[2];
rz(0.53831354) q[2];
rz(2.8484143) q[3];
sx q[3];
rz(-1.5978975) q[3];
sx q[3];
rz(-1.258491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20741589) q[0];
sx q[0];
rz(-0.84686142) q[0];
sx q[0];
rz(-2.9127981) q[0];
rz(-3.0711807) q[1];
sx q[1];
rz(-1.8682624) q[1];
sx q[1];
rz(2.0796897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7477729) q[0];
sx q[0];
rz(-2.5508587) q[0];
sx q[0];
rz(-0.48818199) q[0];
x q[1];
rz(1.1106165) q[2];
sx q[2];
rz(-2.0041582) q[2];
sx q[2];
rz(0.037511911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8613597) q[1];
sx q[1];
rz(-1.2462843) q[1];
sx q[1];
rz(1.8779511) q[1];
rz(-pi) q[2];
rz(-1.358374) q[3];
sx q[3];
rz(-1.1512892) q[3];
sx q[3];
rz(2.7537997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.215302) q[2];
sx q[2];
rz(-0.44507241) q[2];
sx q[2];
rz(0.11623795) q[2];
rz(-0.91533533) q[3];
sx q[3];
rz(-1.4696308) q[3];
sx q[3];
rz(0.62469971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4794469) q[0];
sx q[0];
rz(-1.5597458) q[0];
sx q[0];
rz(2.8060655) q[0];
rz(1.4115931) q[1];
sx q[1];
rz(-0.84452191) q[1];
sx q[1];
rz(-1.9427049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39349651) q[0];
sx q[0];
rz(-0.61433799) q[0];
sx q[0];
rz(-2.0779646) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9556276) q[2];
sx q[2];
rz(-2.129938) q[2];
sx q[2];
rz(2.2988479) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4863805) q[1];
sx q[1];
rz(-0.81334121) q[1];
sx q[1];
rz(1.7310623) q[1];
x q[2];
rz(-2.059778) q[3];
sx q[3];
rz(-1.7918402) q[3];
sx q[3];
rz(-2.3810015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9592643) q[2];
sx q[2];
rz(-2.3458643) q[2];
sx q[2];
rz(-2.0286593) q[2];
rz(2.6304604) q[3];
sx q[3];
rz(-1.3730349) q[3];
sx q[3];
rz(0.50890499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48753259) q[0];
sx q[0];
rz(-0.075981058) q[0];
sx q[0];
rz(2.7677166) q[0];
rz(-1.182425) q[1];
sx q[1];
rz(-1.9805209) q[1];
sx q[1];
rz(-0.16071308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041695853) q[0];
sx q[0];
rz(-1.6536435) q[0];
sx q[0];
rz(-0.13819466) q[0];
rz(0.89779519) q[2];
sx q[2];
rz(-2.2429969) q[2];
sx q[2];
rz(-2.3703645) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7294893) q[1];
sx q[1];
rz(-1.0588875) q[1];
sx q[1];
rz(-2.6964534) q[1];
rz(-pi) q[2];
rz(-2.8237958) q[3];
sx q[3];
rz(-1.8565912) q[3];
sx q[3];
rz(-2.5366549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.60761991) q[2];
sx q[2];
rz(-2.0325568) q[2];
sx q[2];
rz(-2.4268761) q[2];
rz(-0.25501069) q[3];
sx q[3];
rz(-1.3304354) q[3];
sx q[3];
rz(2.3431006) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0890927) q[0];
sx q[0];
rz(-2.8369501) q[0];
sx q[0];
rz(-2.3783045) q[0];
rz(2.8870562) q[1];
sx q[1];
rz(-1.7691879) q[1];
sx q[1];
rz(1.4483784) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7970552) q[0];
sx q[0];
rz(-1.2303924) q[0];
sx q[0];
rz(-1.1974687) q[0];
rz(-1.4035475) q[2];
sx q[2];
rz(-0.75772253) q[2];
sx q[2];
rz(3.0191985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9824286) q[1];
sx q[1];
rz(-2.4953173) q[1];
sx q[1];
rz(-0.29781945) q[1];
x q[2];
rz(2.6877787) q[3];
sx q[3];
rz(-2.3771913) q[3];
sx q[3];
rz(-1.3168471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86108565) q[2];
sx q[2];
rz(-2.0528767) q[2];
sx q[2];
rz(1.0293993) q[2];
rz(-1.6086802) q[3];
sx q[3];
rz(-0.89919388) q[3];
sx q[3];
rz(-2.5077584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3935881) q[0];
sx q[0];
rz(-2.2279976) q[0];
sx q[0];
rz(0.21011259) q[0];
rz(0.70136079) q[1];
sx q[1];
rz(-1.7751834) q[1];
sx q[1];
rz(-2.0393541) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2909246) q[0];
sx q[0];
rz(-1.0449261) q[0];
sx q[0];
rz(2.6024502) q[0];
rz(2.4793825) q[2];
sx q[2];
rz(-1.5689611) q[2];
sx q[2];
rz(2.3037132) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6893856) q[1];
sx q[1];
rz(-1.2579134) q[1];
sx q[1];
rz(0.72231267) q[1];
rz(-pi) q[2];
rz(-2.5999448) q[3];
sx q[3];
rz(-1.4337501) q[3];
sx q[3];
rz(-0.26462072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8474951) q[2];
sx q[2];
rz(-1.5832486) q[2];
sx q[2];
rz(-1.6451277) q[2];
rz(1.7485113) q[3];
sx q[3];
rz(-0.93904177) q[3];
sx q[3];
rz(2.4699672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74712718) q[0];
sx q[0];
rz(-1.324993) q[0];
sx q[0];
rz(2.6649244) q[0];
rz(1.3564302) q[1];
sx q[1];
rz(-1.2973659) q[1];
sx q[1];
rz(2.2043601) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7362226) q[0];
sx q[0];
rz(-2.4511466) q[0];
sx q[0];
rz(1.5140945) q[0];
rz(-pi) q[1];
x q[1];
rz(1.632948) q[2];
sx q[2];
rz(-2.2878433) q[2];
sx q[2];
rz(-0.99822068) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7644219) q[1];
sx q[1];
rz(-0.57767361) q[1];
sx q[1];
rz(2.38297) q[1];
rz(-pi) q[2];
rz(-0.44853278) q[3];
sx q[3];
rz(-1.483102) q[3];
sx q[3];
rz(1.0486856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6095907) q[2];
sx q[2];
rz(-1.3401745) q[2];
sx q[2];
rz(1.0756005) q[2];
rz(0.39014751) q[3];
sx q[3];
rz(-1.8307999) q[3];
sx q[3];
rz(0.80316097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53668642) q[0];
sx q[0];
rz(-1.0746047) q[0];
sx q[0];
rz(0.54501504) q[0];
rz(2.1489428) q[1];
sx q[1];
rz(-2.1558709) q[1];
sx q[1];
rz(-1.4535905) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75460583) q[0];
sx q[0];
rz(-2.0029066) q[0];
sx q[0];
rz(-2.291392) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1686859) q[2];
sx q[2];
rz(-2.997749) q[2];
sx q[2];
rz(-2.3013702) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0664898) q[1];
sx q[1];
rz(-2.3716092) q[1];
sx q[1];
rz(-1.7335197) q[1];
x q[2];
rz(-2.4798315) q[3];
sx q[3];
rz(-0.076120928) q[3];
sx q[3];
rz(-0.95835987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7327026) q[2];
sx q[2];
rz(-1.7444892) q[2];
sx q[2];
rz(-0.2962386) q[2];
rz(-2.564751) q[3];
sx q[3];
rz(-0.39671612) q[3];
sx q[3];
rz(1.2805987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.3532362) q[0];
sx q[0];
rz(-0.78859538) q[0];
sx q[0];
rz(0.22407918) q[0];
rz(2.5375986) q[1];
sx q[1];
rz(-2.150034) q[1];
sx q[1];
rz(-1.6557065) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7753405) q[0];
sx q[0];
rz(-1.9320017) q[0];
sx q[0];
rz(-1.3068178) q[0];
rz(-pi) q[1];
rz(-2.2182216) q[2];
sx q[2];
rz(-0.24590547) q[2];
sx q[2];
rz(2.3594423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60843492) q[1];
sx q[1];
rz(-0.67109334) q[1];
sx q[1];
rz(-2.72481) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2717683) q[3];
sx q[3];
rz(-1.097603) q[3];
sx q[3];
rz(0.92834559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6622582) q[2];
sx q[2];
rz(-2.0749638) q[2];
sx q[2];
rz(-0.38702854) q[2];
rz(2.8582063) q[3];
sx q[3];
rz(-2.3895388) q[3];
sx q[3];
rz(-0.40248218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1166444) q[0];
sx q[0];
rz(-2.8706757) q[0];
sx q[0];
rz(1.1234294) q[0];
rz(-1.4943538) q[1];
sx q[1];
rz(-1.6554183) q[1];
sx q[1];
rz(2.2656238) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6982684) q[0];
sx q[0];
rz(-0.7900376) q[0];
sx q[0];
rz(-2.7774071) q[0];
rz(-pi) q[1];
rz(-2.1112305) q[2];
sx q[2];
rz(-1.9649995) q[2];
sx q[2];
rz(-2.9351007) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3800602) q[1];
sx q[1];
rz(-2.426154) q[1];
sx q[1];
rz(-1.9081294) q[1];
rz(-pi) q[2];
rz(-1.749529) q[3];
sx q[3];
rz(-1.1471841) q[3];
sx q[3];
rz(-0.42631794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4103849) q[2];
sx q[2];
rz(-1.3599675) q[2];
sx q[2];
rz(3.0044921) q[2];
rz(1.3827263) q[3];
sx q[3];
rz(-2.2189249) q[3];
sx q[3];
rz(-1.9130798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41505861) q[0];
sx q[0];
rz(-0.4701604) q[0];
sx q[0];
rz(-0.62248019) q[0];
rz(2.7525735) q[1];
sx q[1];
rz(-0.35094378) q[1];
sx q[1];
rz(2.6019179) q[1];
rz(-0.48609235) q[2];
sx q[2];
rz(-1.9584502) q[2];
sx q[2];
rz(2.9050585) q[2];
rz(-2.9403654) q[3];
sx q[3];
rz(-1.9126127) q[3];
sx q[3];
rz(0.9780608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
