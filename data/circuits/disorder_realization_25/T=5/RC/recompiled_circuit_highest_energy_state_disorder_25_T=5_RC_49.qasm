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
rz(-0.93936062) q[0];
sx q[0];
rz(-2.0225749) q[0];
sx q[0];
rz(0.033666704) q[0];
rz(-1.4576003) q[1];
sx q[1];
rz(-1.8063318) q[1];
sx q[1];
rz(-1.6330947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7608838) q[0];
sx q[0];
rz(-0.11822001) q[0];
sx q[0];
rz(1.8272754) q[0];
x q[1];
rz(-1.2766507) q[2];
sx q[2];
rz(-1.642881) q[2];
sx q[2];
rz(-2.8800137) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37225393) q[1];
sx q[1];
rz(-2.4318691) q[1];
sx q[1];
rz(0.22270568) q[1];
rz(0.90773031) q[3];
sx q[3];
rz(-2.2624348) q[3];
sx q[3];
rz(1.7691106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9713251) q[2];
sx q[2];
rz(-0.89189947) q[2];
sx q[2];
rz(-0.83667052) q[2];
rz(2.3352046) q[3];
sx q[3];
rz(-0.92525768) q[3];
sx q[3];
rz(-0.18536082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0760913) q[0];
sx q[0];
rz(-1.0668904) q[0];
sx q[0];
rz(-1.2170894) q[0];
rz(-2.941046) q[1];
sx q[1];
rz(-2.3417818) q[1];
sx q[1];
rz(-2.5977871) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4763308) q[0];
sx q[0];
rz(-1.7996739) q[0];
sx q[0];
rz(-0.32461597) q[0];
x q[1];
rz(-0.5949623) q[2];
sx q[2];
rz(-1.5820855) q[2];
sx q[2];
rz(-1.4608698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75419352) q[1];
sx q[1];
rz(-2.6541944) q[1];
sx q[1];
rz(0.91895603) q[1];
rz(-pi) q[2];
rz(1.614566) q[3];
sx q[3];
rz(-0.67568103) q[3];
sx q[3];
rz(-1.1202433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4726938) q[2];
sx q[2];
rz(-0.74287477) q[2];
sx q[2];
rz(-1.6061858) q[2];
rz(-1.499048) q[3];
sx q[3];
rz(-2.5583772) q[3];
sx q[3];
rz(-1.2523874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.1010308) q[0];
sx q[0];
rz(-0.56483785) q[0];
sx q[0];
rz(-3.132013) q[0];
rz(2.1122232) q[1];
sx q[1];
rz(-1.5734438) q[1];
sx q[1];
rz(-1.8052489) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1642417) q[0];
sx q[0];
rz(-2.079043) q[0];
sx q[0];
rz(-2.4524816) q[0];
rz(-0.32149489) q[2];
sx q[2];
rz(-1.090637) q[2];
sx q[2];
rz(0.26115044) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6402898) q[1];
sx q[1];
rz(-2.7510018) q[1];
sx q[1];
rz(-1.0949446) q[1];
rz(-pi) q[2];
rz(-1.9591145) q[3];
sx q[3];
rz(-0.42436436) q[3];
sx q[3];
rz(-0.96411112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10793992) q[2];
sx q[2];
rz(-0.5113217) q[2];
sx q[2];
rz(1.8281724) q[2];
rz(-0.039947346) q[3];
sx q[3];
rz(-0.91244709) q[3];
sx q[3];
rz(1.1279172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0082598) q[0];
sx q[0];
rz(-1.2902211) q[0];
sx q[0];
rz(-2.7449352) q[0];
rz(2.2465514) q[1];
sx q[1];
rz(-1.7170693) q[1];
sx q[1];
rz(2.3779714) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9485003) q[0];
sx q[0];
rz(-1.4903127) q[0];
sx q[0];
rz(-1.3250687) q[0];
x q[1];
rz(-3.0210546) q[2];
sx q[2];
rz(-1.5876895) q[2];
sx q[2];
rz(2.8617045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3862003) q[1];
sx q[1];
rz(-1.5441276) q[1];
sx q[1];
rz(1.3332572) q[1];
rz(-pi) q[2];
rz(-0.29080363) q[3];
sx q[3];
rz(-1.606308) q[3];
sx q[3];
rz(-2.5574656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.414173) q[2];
sx q[2];
rz(-1.5199993) q[2];
sx q[2];
rz(-1.1302036) q[2];
rz(1.1388418) q[3];
sx q[3];
rz(-1.9500705) q[3];
sx q[3];
rz(1.9738919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30911699) q[0];
sx q[0];
rz(-2.583355) q[0];
sx q[0];
rz(0.48156893) q[0];
rz(-2.8140977) q[1];
sx q[1];
rz(-1.4188674) q[1];
sx q[1];
rz(-0.96424261) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1314581) q[0];
sx q[0];
rz(-1.1842964) q[0];
sx q[0];
rz(-0.23614998) q[0];
rz(-pi) q[1];
rz(0.19179253) q[2];
sx q[2];
rz(-0.99770498) q[2];
sx q[2];
rz(1.7157784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22589265) q[1];
sx q[1];
rz(-1.4488765) q[1];
sx q[1];
rz(-1.8173161) q[1];
x q[2];
rz(0.39059095) q[3];
sx q[3];
rz(-1.7601305) q[3];
sx q[3];
rz(-2.6238497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97011906) q[2];
sx q[2];
rz(-1.118266) q[2];
sx q[2];
rz(1.6469275) q[2];
rz(-1.0129048) q[3];
sx q[3];
rz(-2.2154112) q[3];
sx q[3];
rz(-2.6265898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2038912) q[0];
sx q[0];
rz(-1.0045811) q[0];
sx q[0];
rz(-2.0830233) q[0];
rz(2.2189498) q[1];
sx q[1];
rz(-1.0170931) q[1];
sx q[1];
rz(0.10291084) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424868) q[0];
sx q[0];
rz(-2.1916336) q[0];
sx q[0];
rz(2.5607361) q[0];
x q[1];
rz(2.5526664) q[2];
sx q[2];
rz(-1.803033) q[2];
sx q[2];
rz(3.1265772) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2051112) q[1];
sx q[1];
rz(-1.0432341) q[1];
sx q[1];
rz(-0.21326201) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2872171) q[3];
sx q[3];
rz(-1.5660785) q[3];
sx q[3];
rz(-1.9172052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4395478) q[2];
sx q[2];
rz(-0.44946686) q[2];
sx q[2];
rz(2.8257545) q[2];
rz(-0.99509197) q[3];
sx q[3];
rz(-1.8683878) q[3];
sx q[3];
rz(-2.4163213) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.134326) q[0];
sx q[0];
rz(-0.26679978) q[0];
sx q[0];
rz(-2.5079978) q[0];
rz(2.7569356) q[1];
sx q[1];
rz(-1.9622012) q[1];
sx q[1];
rz(-2.7124009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10776183) q[0];
sx q[0];
rz(-2.0380479) q[0];
sx q[0];
rz(-0.28470566) q[0];
x q[1];
rz(-0.68558461) q[2];
sx q[2];
rz(-1.597817) q[2];
sx q[2];
rz(2.9179058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0269484) q[1];
sx q[1];
rz(-0.99851313) q[1];
sx q[1];
rz(-1.9932822) q[1];
rz(-pi) q[2];
rz(2.1294534) q[3];
sx q[3];
rz(-1.1691165) q[3];
sx q[3];
rz(-2.5334849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.089243285) q[2];
sx q[2];
rz(-0.77755916) q[2];
sx q[2];
rz(0.26697978) q[2];
rz(-3.0698245) q[3];
sx q[3];
rz(-1.0177344) q[3];
sx q[3];
rz(1.7223541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778022) q[0];
sx q[0];
rz(-3.0402122) q[0];
sx q[0];
rz(-1.4005533) q[0];
rz(-1.1446704) q[1];
sx q[1];
rz(-2.0796227) q[1];
sx q[1];
rz(-0.58715075) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82281821) q[0];
sx q[0];
rz(-2.1923034) q[0];
sx q[0];
rz(1.0448635) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7233974) q[2];
sx q[2];
rz(-1.4640239) q[2];
sx q[2];
rz(-2.622294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59881951) q[1];
sx q[1];
rz(-0.89869281) q[1];
sx q[1];
rz(-1.2818976) q[1];
x q[2];
rz(0.31717698) q[3];
sx q[3];
rz(-0.97517386) q[3];
sx q[3];
rz(0.45444876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8553541) q[2];
sx q[2];
rz(-1.9619532) q[2];
sx q[2];
rz(0.41909763) q[2];
rz(-1.5239483) q[3];
sx q[3];
rz(-2.2231299) q[3];
sx q[3];
rz(-1.6379448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54290259) q[0];
sx q[0];
rz(-1.9200696) q[0];
sx q[0];
rz(0.30938095) q[0];
rz(1.7912553) q[1];
sx q[1];
rz(-2.5534936) q[1];
sx q[1];
rz(2.5877171) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079289201) q[0];
sx q[0];
rz(-1.252018) q[0];
sx q[0];
rz(-0.43894025) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1872841) q[2];
sx q[2];
rz(-1.4330136) q[2];
sx q[2];
rz(2.1509508) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83912056) q[1];
sx q[1];
rz(-2.1697147) q[1];
sx q[1];
rz(0.37521404) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7782544) q[3];
sx q[3];
rz(-2.4060892) q[3];
sx q[3];
rz(2.9309762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.89173335) q[2];
sx q[2];
rz(-1.5313818) q[2];
sx q[2];
rz(-1.7943133) q[2];
rz(-0.76752082) q[3];
sx q[3];
rz(-1.9556421) q[3];
sx q[3];
rz(0.84848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42221853) q[0];
sx q[0];
rz(-2.8415871) q[0];
sx q[0];
rz(1.2147709) q[0];
rz(2.1396554) q[1];
sx q[1];
rz(-1.6286214) q[1];
sx q[1];
rz(-2.217206) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0262605) q[0];
sx q[0];
rz(-0.9965653) q[0];
sx q[0];
rz(2.6747245) q[0];
rz(-0.027410313) q[2];
sx q[2];
rz(-1.0263342) q[2];
sx q[2];
rz(1.4278864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0707911) q[1];
sx q[1];
rz(-1.3475349) q[1];
sx q[1];
rz(-1.5469458) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2460292) q[3];
sx q[3];
rz(-2.319284) q[3];
sx q[3];
rz(0.12982097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34710106) q[2];
sx q[2];
rz(-0.93239409) q[2];
sx q[2];
rz(-0.36716983) q[2];
rz(1.1651039) q[3];
sx q[3];
rz(-0.96787435) q[3];
sx q[3];
rz(-0.86618209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544871) q[0];
sx q[0];
rz(-0.73910537) q[0];
sx q[0];
rz(0.11866971) q[0];
rz(-0.26609303) q[1];
sx q[1];
rz(-2.2292744) q[1];
sx q[1];
rz(1.6899065) q[1];
rz(-2.8546147) q[2];
sx q[2];
rz(-1.8331883) q[2];
sx q[2];
rz(2.2856648) q[2];
rz(-0.96919555) q[3];
sx q[3];
rz(-1.7432913) q[3];
sx q[3];
rz(-0.023719214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
