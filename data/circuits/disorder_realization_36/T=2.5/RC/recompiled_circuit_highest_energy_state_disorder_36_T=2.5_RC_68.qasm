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
rz(3.1004768) q[0];
sx q[0];
rz(-1.2007204) q[0];
sx q[0];
rz(-1.6706985) q[0];
rz(-0.38493758) q[1];
sx q[1];
rz(-0.5783143) q[1];
sx q[1];
rz(-0.90775031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0484412) q[0];
sx q[0];
rz(-1.0142097) q[0];
sx q[0];
rz(2.0970035) q[0];
x q[1];
rz(1.8706246) q[2];
sx q[2];
rz(-1.289164) q[2];
sx q[2];
rz(2.9099438) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5564087) q[1];
sx q[1];
rz(-0.5250465) q[1];
sx q[1];
rz(-0.37918143) q[1];
x q[2];
rz(-1.5159881) q[3];
sx q[3];
rz(-1.6920838) q[3];
sx q[3];
rz(0.19586043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3140728) q[2];
sx q[2];
rz(-1.4853518) q[2];
sx q[2];
rz(2.8208288) q[2];
rz(3.1406506) q[3];
sx q[3];
rz(-0.92206803) q[3];
sx q[3];
rz(2.6607386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8568521) q[0];
sx q[0];
rz(-2.4212403) q[0];
sx q[0];
rz(0.62069935) q[0];
rz(1.8547828) q[1];
sx q[1];
rz(-1.4871253) q[1];
sx q[1];
rz(-0.99498814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5704203) q[0];
sx q[0];
rz(-1.3745527) q[0];
sx q[0];
rz(-1.061383) q[0];
rz(2.1844145) q[2];
sx q[2];
rz(-2.6189096) q[2];
sx q[2];
rz(-1.7651287) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2590005) q[1];
sx q[1];
rz(-1.126312) q[1];
sx q[1];
rz(-2.9372272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6334559) q[3];
sx q[3];
rz(-2.2656815) q[3];
sx q[3];
rz(2.9949403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2113125) q[2];
sx q[2];
rz(-1.2418396) q[2];
sx q[2];
rz(2.9105183) q[2];
rz(1.6127582) q[3];
sx q[3];
rz(-1.684618) q[3];
sx q[3];
rz(0.55135623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7174317) q[0];
sx q[0];
rz(-1.1671952) q[0];
sx q[0];
rz(-0.396808) q[0];
rz(1.1948168) q[1];
sx q[1];
rz(-1.8794941) q[1];
sx q[1];
rz(1.6669115) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68434381) q[0];
sx q[0];
rz(-1.5822344) q[0];
sx q[0];
rz(1.6586763) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0662768) q[2];
sx q[2];
rz(-1.0766509) q[2];
sx q[2];
rz(-2.0574428) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.35016216) q[1];
sx q[1];
rz(-1.9099978) q[1];
sx q[1];
rz(-2.2712525) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5422882) q[3];
sx q[3];
rz(-0.39428821) q[3];
sx q[3];
rz(-0.35343227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.84612238) q[2];
sx q[2];
rz(-2.2961871) q[2];
sx q[2];
rz(1.6579312) q[2];
rz(2.1814003) q[3];
sx q[3];
rz(-0.63513297) q[3];
sx q[3];
rz(0.6944164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33261499) q[0];
sx q[0];
rz(-1.806145) q[0];
sx q[0];
rz(-0.72840011) q[0];
rz(-1.867846) q[1];
sx q[1];
rz(-1.961901) q[1];
sx q[1];
rz(1.5789998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1611164) q[0];
sx q[0];
rz(-3.118782) q[0];
sx q[0];
rz(-1.7177607) q[0];
x q[1];
rz(2.118763) q[2];
sx q[2];
rz(-1.5336516) q[2];
sx q[2];
rz(1.9859344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82484805) q[1];
sx q[1];
rz(-1.4364752) q[1];
sx q[1];
rz(2.9243241) q[1];
rz(-1.9779491) q[3];
sx q[3];
rz(-2.4302357) q[3];
sx q[3];
rz(-2.9748981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1100715) q[2];
sx q[2];
rz(-1.4036125) q[2];
sx q[2];
rz(2.9384379) q[2];
rz(1.7848232) q[3];
sx q[3];
rz(-1.9316831) q[3];
sx q[3];
rz(1.8206966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.2585231) q[0];
sx q[0];
rz(-1.3354381) q[0];
sx q[0];
rz(1.6254599) q[0];
rz(2.7698703) q[1];
sx q[1];
rz(-0.55551353) q[1];
sx q[1];
rz(-2.1518167) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55732841) q[0];
sx q[0];
rz(-0.66270486) q[0];
sx q[0];
rz(1.8813547) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6838275) q[2];
sx q[2];
rz(-0.80146061) q[2];
sx q[2];
rz(-2.46704) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6367148) q[1];
sx q[1];
rz(-1.2413919) q[1];
sx q[1];
rz(-0.18784951) q[1];
x q[2];
rz(1.8797713) q[3];
sx q[3];
rz(-1.1011964) q[3];
sx q[3];
rz(-2.0275786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61701361) q[2];
sx q[2];
rz(-2.9144574) q[2];
sx q[2];
rz(-0.058509286) q[2];
rz(-1.2837563) q[3];
sx q[3];
rz(-2.3144898) q[3];
sx q[3];
rz(-1.1684928) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44130317) q[0];
sx q[0];
rz(-2.0952201) q[0];
sx q[0];
rz(-2.6562279) q[0];
rz(-0.46733388) q[1];
sx q[1];
rz(-2.554437) q[1];
sx q[1];
rz(2.4553518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29083911) q[0];
sx q[0];
rz(-2.3510482) q[0];
sx q[0];
rz(-2.683862) q[0];
x q[1];
rz(-1.2755308) q[2];
sx q[2];
rz(-0.64744189) q[2];
sx q[2];
rz(1.9290627) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8279449) q[1];
sx q[1];
rz(-2.0454496) q[1];
sx q[1];
rz(-0.56874911) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3090906) q[3];
sx q[3];
rz(-1.0847391) q[3];
sx q[3];
rz(1.6532236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83963362) q[2];
sx q[2];
rz(-1.0884103) q[2];
sx q[2];
rz(0.23183091) q[2];
rz(0.92136446) q[3];
sx q[3];
rz(-1.6145584) q[3];
sx q[3];
rz(-0.025402633) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63474083) q[0];
sx q[0];
rz(-1.0036108) q[0];
sx q[0];
rz(-2.4941709) q[0];
rz(-2.1000775) q[1];
sx q[1];
rz(-1.2728649) q[1];
sx q[1];
rz(-1.1127068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0395687) q[0];
sx q[0];
rz(-1.0728076) q[0];
sx q[0];
rz(0.51135333) q[0];
x q[1];
rz(-1.9693807) q[2];
sx q[2];
rz(-2.9867134) q[2];
sx q[2];
rz(-1.1501555) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0354516) q[1];
sx q[1];
rz(-1.3961432) q[1];
sx q[1];
rz(-2.5779252) q[1];
rz(0.010257234) q[3];
sx q[3];
rz(-2.0358277) q[3];
sx q[3];
rz(-0.51562515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5645494) q[2];
sx q[2];
rz(-0.81114045) q[2];
sx q[2];
rz(0.88456336) q[2];
rz(1.0386284) q[3];
sx q[3];
rz(-1.1218772) q[3];
sx q[3];
rz(-1.60892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3603947) q[0];
sx q[0];
rz(-2.9304657) q[0];
sx q[0];
rz(1.1398844) q[0];
rz(-0.79849157) q[1];
sx q[1];
rz(-1.4356828) q[1];
sx q[1];
rz(3.0671157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7388106) q[0];
sx q[0];
rz(-2.180679) q[0];
sx q[0];
rz(-0.29891157) q[0];
x q[1];
rz(-0.42718446) q[2];
sx q[2];
rz(-2.3318034) q[2];
sx q[2];
rz(2.0860738) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18596126) q[1];
sx q[1];
rz(-2.4096386) q[1];
sx q[1];
rz(2.6630528) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8387855) q[3];
sx q[3];
rz(-1.9085371) q[3];
sx q[3];
rz(-0.95888058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8188339) q[2];
sx q[2];
rz(-2.5311311) q[2];
sx q[2];
rz(-0.16767137) q[2];
rz(-0.75336114) q[3];
sx q[3];
rz(-1.2246917) q[3];
sx q[3];
rz(-1.0838449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3566256) q[0];
sx q[0];
rz(-1.4545119) q[0];
sx q[0];
rz(1.0474569) q[0];
rz(0.52182237) q[1];
sx q[1];
rz(-1.1818886) q[1];
sx q[1];
rz(-0.80400115) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76011953) q[0];
sx q[0];
rz(-1.9895041) q[0];
sx q[0];
rz(-2.9497854) q[0];
rz(-0.1754403) q[2];
sx q[2];
rz(-2.1198303) q[2];
sx q[2];
rz(-0.4859449) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9404133) q[1];
sx q[1];
rz(-2.6087954) q[1];
sx q[1];
rz(-2.6032107) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6860572) q[3];
sx q[3];
rz(-0.75596327) q[3];
sx q[3];
rz(-0.29827932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1890586) q[2];
sx q[2];
rz(-1.4905832) q[2];
sx q[2];
rz(2.789759) q[2];
rz(0.072988836) q[3];
sx q[3];
rz(-1.9609541) q[3];
sx q[3];
rz(-0.73968691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22013448) q[0];
sx q[0];
rz(-2.298482) q[0];
sx q[0];
rz(1.5699009) q[0];
rz(2.4193071) q[1];
sx q[1];
rz(-2.0000358) q[1];
sx q[1];
rz(1.7438181) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0972526) q[0];
sx q[0];
rz(-1.3201137) q[0];
sx q[0];
rz(-0.41608019) q[0];
x q[1];
rz(2.4358814) q[2];
sx q[2];
rz(-1.7448261) q[2];
sx q[2];
rz(1.1104012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0666982) q[1];
sx q[1];
rz(-0.5579307) q[1];
sx q[1];
rz(2.9900223) q[1];
x q[2];
rz(-0.13171036) q[3];
sx q[3];
rz(-1.9870821) q[3];
sx q[3];
rz(0.10271969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7950644) q[2];
sx q[2];
rz(-1.3115735) q[2];
sx q[2];
rz(1.2947003) q[2];
rz(2.0857701) q[3];
sx q[3];
rz(-2.081223) q[3];
sx q[3];
rz(1.0052217) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9697363) q[0];
sx q[0];
rz(-2.2315401) q[0];
sx q[0];
rz(1.6730614) q[0];
rz(-0.83364529) q[1];
sx q[1];
rz(-1.2794762) q[1];
sx q[1];
rz(1.5247482) q[1];
rz(2.0019309) q[2];
sx q[2];
rz(-2.3464602) q[2];
sx q[2];
rz(-2.4577557) q[2];
rz(0.41224319) q[3];
sx q[3];
rz(-1.1664248) q[3];
sx q[3];
rz(2.5965367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
