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
rz(1.6871356) q[0];
sx q[0];
rz(1.9782826) q[0];
sx q[0];
rz(10.692698) q[0];
rz(-1.4866225) q[1];
sx q[1];
rz(-1.2868737) q[1];
sx q[1];
rz(-0.16965228) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9244512) q[0];
sx q[0];
rz(-2.0471153) q[0];
sx q[0];
rz(2.7360271) q[0];
rz(1.0526161) q[2];
sx q[2];
rz(-1.7917601) q[2];
sx q[2];
rz(-2.882304) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4625068) q[1];
sx q[1];
rz(-0.39606491) q[1];
sx q[1];
rz(1.5157265) q[1];
x q[2];
rz(2.9291487) q[3];
sx q[3];
rz(-1.8768969) q[3];
sx q[3];
rz(2.3116632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9156645) q[2];
sx q[2];
rz(-1.516284) q[2];
sx q[2];
rz(3.1212659) q[2];
rz(-2.9018371) q[3];
sx q[3];
rz(-0.25529796) q[3];
sx q[3];
rz(-2.3026626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13424419) q[0];
sx q[0];
rz(-2.5237995) q[0];
sx q[0];
rz(2.0452621) q[0];
rz(-1.4758551) q[1];
sx q[1];
rz(-1.9225537) q[1];
sx q[1];
rz(2.347167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4969032) q[0];
sx q[0];
rz(-2.6711426) q[0];
sx q[0];
rz(0.75017922) q[0];
rz(-0.6068318) q[2];
sx q[2];
rz(-0.85222679) q[2];
sx q[2];
rz(-2.5963714) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8960171) q[1];
sx q[1];
rz(-1.4652677) q[1];
sx q[1];
rz(-0.42550931) q[1];
x q[2];
rz(2.9555213) q[3];
sx q[3];
rz(-2.3181462) q[3];
sx q[3];
rz(-0.88445437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.085792556) q[2];
sx q[2];
rz(-1.5726568) q[2];
sx q[2];
rz(2.9158578) q[2];
rz(2.8977532) q[3];
sx q[3];
rz(-0.95832458) q[3];
sx q[3];
rz(0.17914151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8826411) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(0.42695811) q[0];
rz(0.34307617) q[1];
sx q[1];
rz(-1.9648569) q[1];
sx q[1];
rz(2.8899736) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8247525) q[0];
sx q[0];
rz(-1.3469633) q[0];
sx q[0];
rz(0.59455183) q[0];
rz(-pi) q[1];
rz(0.80406739) q[2];
sx q[2];
rz(-1.7323426) q[2];
sx q[2];
rz(1.3093914) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15089825) q[1];
sx q[1];
rz(-2.7643235) q[1];
sx q[1];
rz(0.47028415) q[1];
rz(-pi) q[2];
rz(1.512647) q[3];
sx q[3];
rz(-0.71134252) q[3];
sx q[3];
rz(-2.1429246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0106657) q[2];
sx q[2];
rz(-1.6341354) q[2];
sx q[2];
rz(-0.44535401) q[2];
rz(1.917786) q[3];
sx q[3];
rz(-1.2659975) q[3];
sx q[3];
rz(0.38770097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7447516) q[0];
sx q[0];
rz(-0.80533177) q[0];
sx q[0];
rz(-1.9824363) q[0];
rz(1.9388439) q[1];
sx q[1];
rz(-1.1801327) q[1];
sx q[1];
rz(0.73701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7173968) q[0];
sx q[0];
rz(-1.5079466) q[0];
sx q[0];
rz(-0.24817384) q[0];
rz(-2.0194172) q[2];
sx q[2];
rz(-2.4146842) q[2];
sx q[2];
rz(0.87269937) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9056978) q[1];
sx q[1];
rz(-1.7445843) q[1];
sx q[1];
rz(0.64069616) q[1];
rz(-pi) q[2];
rz(2.9748671) q[3];
sx q[3];
rz(-1.6279711) q[3];
sx q[3];
rz(-2.056332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1272588) q[2];
sx q[2];
rz(-0.88982439) q[2];
sx q[2];
rz(-0.8117525) q[2];
rz(2.3938866) q[3];
sx q[3];
rz(-2.2816608) q[3];
sx q[3];
rz(-3.0809793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6478445) q[0];
sx q[0];
rz(-1.0890549) q[0];
sx q[0];
rz(-1.3519979) q[0];
rz(-0.79212517) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(-1.0101213) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3528004) q[0];
sx q[0];
rz(-1.9051587) q[0];
sx q[0];
rz(-1.0092606) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0094658) q[2];
sx q[2];
rz(-2.7458522) q[2];
sx q[2];
rz(0.49116116) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7923342) q[1];
sx q[1];
rz(-2.5790303) q[1];
sx q[1];
rz(2.9922561) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.185818) q[3];
sx q[3];
rz(-0.20032756) q[3];
sx q[3];
rz(-2.1708787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1048364) q[2];
sx q[2];
rz(-0.75883055) q[2];
sx q[2];
rz(-1.7581615) q[2];
rz(0.029959921) q[3];
sx q[3];
rz(-0.99830097) q[3];
sx q[3];
rz(-1.2768607) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6303915) q[0];
sx q[0];
rz(-0.41655219) q[0];
sx q[0];
rz(3.0101486) q[0];
rz(0.99880544) q[1];
sx q[1];
rz(-1.9824332) q[1];
sx q[1];
rz(-1.3712032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5906487) q[0];
sx q[0];
rz(-1.9134132) q[0];
sx q[0];
rz(0.37775535) q[0];
rz(-2.8949685) q[2];
sx q[2];
rz(-2.4402251) q[2];
sx q[2];
rz(-1.3878617) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.882521) q[1];
sx q[1];
rz(-1.4132199) q[1];
sx q[1];
rz(-1.2941918) q[1];
x q[2];
rz(1.5742947) q[3];
sx q[3];
rz(-2.281684) q[3];
sx q[3];
rz(1.9128127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3580631) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(1.0895458) q[2];
rz(-2.774636) q[3];
sx q[3];
rz(-2.7373382) q[3];
sx q[3];
rz(-0.75016108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51339665) q[0];
sx q[0];
rz(-0.80444002) q[0];
sx q[0];
rz(-0.87271571) q[0];
rz(-1.5645507) q[1];
sx q[1];
rz(-1.6460452) q[1];
sx q[1];
rz(0.73211342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.064678) q[0];
sx q[0];
rz(-0.87192649) q[0];
sx q[0];
rz(1.8346321) q[0];
rz(-pi) q[1];
rz(-0.68552981) q[2];
sx q[2];
rz(-1.4741802) q[2];
sx q[2];
rz(-1.8308507) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9804025) q[1];
sx q[1];
rz(-1.4136317) q[1];
sx q[1];
rz(1.0970626) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1752582) q[3];
sx q[3];
rz(-2.8962781) q[3];
sx q[3];
rz(-2.3504864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6096036) q[2];
sx q[2];
rz(-1.0010109) q[2];
sx q[2];
rz(-2.5606489) q[2];
rz(1.1254492) q[3];
sx q[3];
rz(-1.9476798) q[3];
sx q[3];
rz(-1.9862991) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67227143) q[0];
sx q[0];
rz(-3.0348365) q[0];
sx q[0];
rz(-0.17053764) q[0];
rz(-1.8960309) q[1];
sx q[1];
rz(-1.5252557) q[1];
sx q[1];
rz(2.4571498) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.754622) q[0];
sx q[0];
rz(-0.59387654) q[0];
sx q[0];
rz(-2.519849) q[0];
x q[1];
rz(1.4438875) q[2];
sx q[2];
rz(-2.0503245) q[2];
sx q[2];
rz(-0.084464262) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0741006) q[1];
sx q[1];
rz(-1.8163741) q[1];
sx q[1];
rz(-0.26439338) q[1];
rz(1.6309647) q[3];
sx q[3];
rz(-2.2188713) q[3];
sx q[3];
rz(-1.5435406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.838852) q[2];
sx q[2];
rz(-1.3417696) q[2];
sx q[2];
rz(2.0212685) q[2];
rz(-0.66425792) q[3];
sx q[3];
rz(-2.7393326) q[3];
sx q[3];
rz(-2.6334488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.53617322) q[0];
sx q[0];
rz(-0.37261951) q[0];
sx q[0];
rz(0.63825178) q[0];
rz(-3.1328746) q[1];
sx q[1];
rz(-1.1161085) q[1];
sx q[1];
rz(0.4062103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063651872) q[0];
sx q[0];
rz(-1.6839875) q[0];
sx q[0];
rz(-1.8454396) q[0];
rz(-1.4582602) q[2];
sx q[2];
rz(-0.43607682) q[2];
sx q[2];
rz(-1.9364408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6510424) q[1];
sx q[1];
rz(-0.79572751) q[1];
sx q[1];
rz(-2.9814238) q[1];
rz(-pi) q[2];
rz(1.5453734) q[3];
sx q[3];
rz(-1.6442199) q[3];
sx q[3];
rz(-1.5591476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.093988769) q[2];
sx q[2];
rz(-2.0145907) q[2];
sx q[2];
rz(-0.36063933) q[2];
rz(2.4141198) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(0.31807652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3280846) q[0];
sx q[0];
rz(-0.94539517) q[0];
sx q[0];
rz(2.9720921) q[0];
rz(2.4318579) q[1];
sx q[1];
rz(-1.4163481) q[1];
sx q[1];
rz(-1.08606) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57993488) q[0];
sx q[0];
rz(-0.57889639) q[0];
sx q[0];
rz(2.7282342) q[0];
rz(2.8302221) q[2];
sx q[2];
rz(-1.1527921) q[2];
sx q[2];
rz(3.0504464) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5063613) q[1];
sx q[1];
rz(-1.6690134) q[1];
sx q[1];
rz(-2.1327095) q[1];
rz(-0.38137718) q[3];
sx q[3];
rz(-1.3998195) q[3];
sx q[3];
rz(-1.387763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3544932) q[2];
sx q[2];
rz(-3.0566065) q[2];
sx q[2];
rz(-2.7717915) q[2];
rz(-1.9434816) q[3];
sx q[3];
rz(-1.5503649) q[3];
sx q[3];
rz(-0.91355356) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13854606) q[0];
sx q[0];
rz(-1.1644762) q[0];
sx q[0];
rz(-1.2559011) q[0];
rz(2.1239602) q[1];
sx q[1];
rz(-1.7435278) q[1];
sx q[1];
rz(1.7494038) q[1];
rz(-1.6412777) q[2];
sx q[2];
rz(-1.4983836) q[2];
sx q[2];
rz(-0.56624779) q[2];
rz(-0.66987517) q[3];
sx q[3];
rz(-2.9333339) q[3];
sx q[3];
rz(0.40706709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
