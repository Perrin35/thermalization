OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1922167) q[0];
sx q[0];
rz(-2.0944216) q[0];
sx q[0];
rz(3.0728683) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(-1.9332164) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2657304) q[0];
sx q[0];
rz(-3.0348572) q[0];
sx q[0];
rz(-2.021832) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0157929) q[2];
sx q[2];
rz(-1.558779) q[2];
sx q[2];
rz(-1.2596241) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5282643) q[1];
sx q[1];
rz(-1.4600888) q[1];
sx q[1];
rz(-2.5221591) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73322202) q[3];
sx q[3];
rz(-2.4823722) q[3];
sx q[3];
rz(-0.99430195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3258813) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(-1.4665843) q[2];
rz(-2.4438434) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(0.74716032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0242457) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(1.9673989) q[0];
rz(-0.17114561) q[1];
sx q[1];
rz(-1.0447964) q[1];
sx q[1];
rz(0.29719621) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575465) q[0];
sx q[0];
rz(-1.6539126) q[0];
sx q[0];
rz(-3.0931285) q[0];
rz(0.71248033) q[2];
sx q[2];
rz(-1.9409632) q[2];
sx q[2];
rz(-0.21288255) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9651523) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(-2.7494207) q[1];
rz(0.76241775) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(-0.03014119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.979636) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(-0.87614122) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0959594) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(2.8339548) q[0];
rz(0.74854198) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(-0.83980733) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1727985) q[0];
sx q[0];
rz(-1.8546687) q[0];
sx q[0];
rz(0.42885212) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33939056) q[2];
sx q[2];
rz(-1.4815785) q[2];
sx q[2];
rz(2.3238376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3629344) q[1];
sx q[1];
rz(-2.0211453) q[1];
sx q[1];
rz(1.2308916) q[1];
rz(-pi) q[2];
x q[2];
rz(2.967756) q[3];
sx q[3];
rz(-0.94508119) q[3];
sx q[3];
rz(-3.0415149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.7904736) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(-1.2157724) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(2.3777305) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(0.44556251) q[0];
rz(0.62082779) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-2.1760118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8863106) q[0];
sx q[0];
rz(-2.3112486) q[0];
sx q[0];
rz(-1.1925973) q[0];
x q[1];
rz(0.84237174) q[2];
sx q[2];
rz(-1.5475376) q[2];
sx q[2];
rz(2.1881441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6009778) q[1];
sx q[1];
rz(-1.8243196) q[1];
sx q[1];
rz(-0.032932245) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2772917) q[3];
sx q[3];
rz(-2.9923277) q[3];
sx q[3];
rz(-2.5316558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(-2.2122673) q[2];
rz(0.6435414) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(1.8575645) q[0];
rz(-2.8517826) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(1.0481542) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6402733) q[0];
sx q[0];
rz(-0.87448705) q[0];
sx q[0];
rz(2.2107844) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58446144) q[2];
sx q[2];
rz(-2.2197154) q[2];
sx q[2];
rz(-1.4635758) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2939799) q[1];
sx q[1];
rz(-1.8925397) q[1];
sx q[1];
rz(-1.4732248) q[1];
rz(-pi) q[2];
rz(-0.9976451) q[3];
sx q[3];
rz(-0.76612681) q[3];
sx q[3];
rz(1.0539953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8706878) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-2.2407545) q[2];
rz(-1.0926931) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(-1.2341011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822405) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(2.545488) q[0];
rz(1.6456564) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.8966819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67958528) q[0];
sx q[0];
rz(-1.8903268) q[0];
sx q[0];
rz(1.0380448) q[0];
rz(2.5858324) q[2];
sx q[2];
rz(-0.55609497) q[2];
sx q[2];
rz(0.13242002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7024755) q[1];
sx q[1];
rz(-0.87279746) q[1];
sx q[1];
rz(0.67081397) q[1];
x q[2];
rz(-1.8025493) q[3];
sx q[3];
rz(-1.8125121) q[3];
sx q[3];
rz(-1.4178604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.231455) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(-0.088430017) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(2.6627873) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46463075) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(-0.27994573) q[0];
rz(-1.4631368) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(0.25269145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8425927) q[0];
sx q[0];
rz(-2.8447066) q[0];
sx q[0];
rz(1.7345558) q[0];
rz(0.93637034) q[2];
sx q[2];
rz(-1.5555824) q[2];
sx q[2];
rz(-2.0903367) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6934769) q[1];
sx q[1];
rz(-2.2709393) q[1];
sx q[1];
rz(-2.4561988) q[1];
x q[2];
rz(0.85315506) q[3];
sx q[3];
rz(-0.95082885) q[3];
sx q[3];
rz(-0.48228797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8213886) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(1.5318711) q[2];
rz(-1.948471) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(-2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10483345) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(-1.4861134) q[0];
rz(-0.2688109) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(0.2789467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3725961) q[0];
sx q[0];
rz(-2.0645752) q[0];
sx q[0];
rz(0.84817024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0686764) q[2];
sx q[2];
rz(-1.4197822) q[2];
sx q[2];
rz(-0.47889027) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3153509) q[1];
sx q[1];
rz(-2.3126174) q[1];
sx q[1];
rz(-0.51893236) q[1];
rz(-pi) q[2];
rz(0.99854462) q[3];
sx q[3];
rz(-1.6958941) q[3];
sx q[3];
rz(2.3541114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3387317) q[2];
sx q[2];
rz(-1.4638476) q[2];
sx q[2];
rz(-1.5489244) q[2];
rz(2.0907949) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.46362296) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(-1.2532225) q[0];
rz(-0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.1463096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63908731) q[0];
sx q[0];
rz(-0.82477942) q[0];
sx q[0];
rz(1.4750925) q[0];
rz(-0.65593221) q[2];
sx q[2];
rz(-2.1428875) q[2];
sx q[2];
rz(-0.42665542) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.031530023) q[1];
sx q[1];
rz(-1.4974125) q[1];
sx q[1];
rz(2.258638) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9911489) q[3];
sx q[3];
rz(-0.97878362) q[3];
sx q[3];
rz(-2.1136485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9877732) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(2.1257607) q[2];
rz(-2.231797) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(0.01424271) q[0];
rz(-2.2968538) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(-1.3815809) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81603564) q[0];
sx q[0];
rz(-2.9986458) q[0];
sx q[0];
rz(1.4965034) q[0];
rz(-pi) q[1];
rz(0.99190418) q[2];
sx q[2];
rz(-1.5169946) q[2];
sx q[2];
rz(-2.3827041) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33243079) q[1];
sx q[1];
rz(-2.5659487) q[1];
sx q[1];
rz(1.8368506) q[1];
x q[2];
rz(1.3515477) q[3];
sx q[3];
rz(-1.0478813) q[3];
sx q[3];
rz(3.074309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0766729) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(-0.98999611) q[2];
rz(2.2475217) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(-2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0513231) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(1.8016215) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(2.0409394) q[2];
sx q[2];
rz(-2.2307776) q[2];
sx q[2];
rz(-2.270436) q[2];
rz(-0.99466956) q[3];
sx q[3];
rz(-0.70636402) q[3];
sx q[3];
rz(-0.93887226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];