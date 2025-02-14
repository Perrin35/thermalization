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
rz(-1.3462525) q[0];
sx q[0];
rz(-1.4486382) q[0];
sx q[0];
rz(-1.9899272) q[0];
rz(1.0898074) q[1];
sx q[1];
rz(-1.4767708) q[1];
sx q[1];
rz(-2.9906315) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7825086) q[0];
sx q[0];
rz(-1.623296) q[0];
sx q[0];
rz(1.2760602) q[0];
rz(-pi) q[1];
rz(1.6577166) q[2];
sx q[2];
rz(-1.2582694) q[2];
sx q[2];
rz(0.0258044) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6358984) q[1];
sx q[1];
rz(-1.555996) q[1];
sx q[1];
rz(-3.1379051) q[1];
rz(-pi) q[2];
rz(-0.096681194) q[3];
sx q[3];
rz(-1.5301609) q[3];
sx q[3];
rz(1.7698897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1208518) q[2];
sx q[2];
rz(-3.1181702) q[2];
sx q[2];
rz(-0.49129301) q[2];
rz(1.8190207) q[3];
sx q[3];
rz(-1.5343851) q[3];
sx q[3];
rz(1.3278495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78111929) q[0];
sx q[0];
rz(-2.5871215) q[0];
sx q[0];
rz(2.4419899) q[0];
rz(1.5910925) q[1];
sx q[1];
rz(-2.6467549) q[1];
sx q[1];
rz(-2.9717305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79029523) q[0];
sx q[0];
rz(-1.6093495) q[0];
sx q[0];
rz(3.0638873) q[0];
x q[1];
rz(-1.4361977) q[2];
sx q[2];
rz(-1.2685809) q[2];
sx q[2];
rz(-0.58537358) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58386034) q[1];
sx q[1];
rz(-0.61994821) q[1];
sx q[1];
rz(0.9264733) q[1];
rz(-3.0624495) q[3];
sx q[3];
rz(-1.1149661) q[3];
sx q[3];
rz(-1.7439708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26213172) q[2];
sx q[2];
rz(-2.7792271) q[2];
sx q[2];
rz(-1.4169089) q[2];
rz(-1.0828241) q[3];
sx q[3];
rz(-0.88734752) q[3];
sx q[3];
rz(0.82720238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5812434) q[0];
sx q[0];
rz(-2.2214948) q[0];
sx q[0];
rz(1.5823407) q[0];
rz(2.4371367) q[1];
sx q[1];
rz(-1.9943941) q[1];
sx q[1];
rz(-2.6047883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0999609) q[0];
sx q[0];
rz(-0.88425626) q[0];
sx q[0];
rz(-1.6306535) q[0];
x q[1];
rz(-0.38449826) q[2];
sx q[2];
rz(-3.0442723) q[2];
sx q[2];
rz(-2.26869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0743064) q[1];
sx q[1];
rz(-0.18504772) q[1];
sx q[1];
rz(-0.61850278) q[1];
rz(-pi) q[2];
rz(2.267089) q[3];
sx q[3];
rz(-1.3811967) q[3];
sx q[3];
rz(0.43280552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74865666) q[2];
sx q[2];
rz(-1.8042678) q[2];
sx q[2];
rz(-0.7788457) q[2];
rz(-0.084176453) q[3];
sx q[3];
rz(-0.86936969) q[3];
sx q[3];
rz(1.4028153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2445143) q[0];
sx q[0];
rz(-2.2205413) q[0];
sx q[0];
rz(0.47065863) q[0];
rz(-1.6828407) q[1];
sx q[1];
rz(-0.0048480821) q[1];
sx q[1];
rz(-0.86476129) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03173657) q[0];
sx q[0];
rz(-1.6288647) q[0];
sx q[0];
rz(2.6266897) q[0];
rz(-pi) q[1];
rz(-1.7030348) q[2];
sx q[2];
rz(-1.7900428) q[2];
sx q[2];
rz(-1.2546778) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4196849) q[1];
sx q[1];
rz(-3.0583755) q[1];
sx q[1];
rz(-0.91099847) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49089499) q[3];
sx q[3];
rz(-2.5449736) q[3];
sx q[3];
rz(0.20363775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1137769) q[2];
sx q[2];
rz(-2.2769589) q[2];
sx q[2];
rz(1.9846385) q[2];
rz(-0.35465869) q[3];
sx q[3];
rz(-1.2361453) q[3];
sx q[3];
rz(-3.0310596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.329634) q[0];
sx q[0];
rz(-2.2060153) q[0];
sx q[0];
rz(-2.6079566) q[0];
rz(1.2136906) q[1];
sx q[1];
rz(-0.023093725) q[1];
sx q[1];
rz(-1.1812706) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2642519) q[0];
sx q[0];
rz(-0.16942693) q[0];
sx q[0];
rz(-2.8657783) q[0];
x q[1];
rz(1.1540765) q[2];
sx q[2];
rz(-0.96229711) q[2];
sx q[2];
rz(1.5173943) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5983728) q[1];
sx q[1];
rz(-2.4717165) q[1];
sx q[1];
rz(1.3175439) q[1];
rz(-0.13792928) q[3];
sx q[3];
rz(-1.0702881) q[3];
sx q[3];
rz(-0.59964666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.10826762) q[2];
sx q[2];
rz(-0.51563087) q[2];
sx q[2];
rz(-2.2132614) q[2];
rz(0.27836529) q[3];
sx q[3];
rz(-1.8230349) q[3];
sx q[3];
rz(0.068647169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5334897) q[0];
sx q[0];
rz(-2.6187496) q[0];
sx q[0];
rz(0.19304092) q[0];
rz(-0.67735425) q[1];
sx q[1];
rz(-3.1132128) q[1];
sx q[1];
rz(-1.631564) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2264435) q[0];
sx q[0];
rz(-0.14227501) q[0];
sx q[0];
rz(2.9998542) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4807499) q[2];
sx q[2];
rz(-2.6777732) q[2];
sx q[2];
rz(0.12741379) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.90551299) q[1];
sx q[1];
rz(-1.326099) q[1];
sx q[1];
rz(-2.3062745) q[1];
x q[2];
rz(0.55431788) q[3];
sx q[3];
rz(-0.80762562) q[3];
sx q[3];
rz(-2.9452711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3626927) q[2];
sx q[2];
rz(-2.429246) q[2];
sx q[2];
rz(2.1544797) q[2];
rz(-2.0569862) q[3];
sx q[3];
rz(-0.54607138) q[3];
sx q[3];
rz(0.62091056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5622332) q[0];
sx q[0];
rz(-0.99246228) q[0];
sx q[0];
rz(1.5480504) q[0];
rz(-2.4515737) q[1];
sx q[1];
rz(-0.023612173) q[1];
sx q[1];
rz(1.4366368) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0853135) q[0];
sx q[0];
rz(-1.5271957) q[0];
sx q[0];
rz(0.52012238) q[0];
rz(3.1146554) q[2];
sx q[2];
rz(-2.4566894) q[2];
sx q[2];
rz(-1.8933344) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1093398) q[1];
sx q[1];
rz(-2.5962862) q[1];
sx q[1];
rz(-1.6857573) q[1];
x q[2];
rz(-0.5593188) q[3];
sx q[3];
rz(-2.6169861) q[3];
sx q[3];
rz(0.23422262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37316608) q[2];
sx q[2];
rz(-1.2370141) q[2];
sx q[2];
rz(-2.5567059) q[2];
rz(-1.54555) q[3];
sx q[3];
rz(-1.4228179) q[3];
sx q[3];
rz(-2.2271631) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5698513) q[0];
sx q[0];
rz(-0.91747704) q[0];
sx q[0];
rz(-1.5745987) q[0];
rz(0.51708108) q[1];
sx q[1];
rz(-2.2061429) q[1];
sx q[1];
rz(1.8749974) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0027813214) q[0];
sx q[0];
rz(-1.5673593) q[0];
sx q[0];
rz(0.0050188662) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6041853) q[2];
sx q[2];
rz(-0.88055842) q[2];
sx q[2];
rz(1.3956525) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4609485) q[1];
sx q[1];
rz(-2.1105745) q[1];
sx q[1];
rz(2.7311027) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.35584) q[3];
sx q[3];
rz(-0.17221299) q[3];
sx q[3];
rz(2.7453051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7685585) q[2];
sx q[2];
rz(-0.80058241) q[2];
sx q[2];
rz(-1.8689092) q[2];
rz(0.4396762) q[3];
sx q[3];
rz(-1.1594073) q[3];
sx q[3];
rz(-0.064597733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8398297) q[0];
sx q[0];
rz(-3.1253212) q[0];
sx q[0];
rz(-0.29626688) q[0];
rz(-3.0773194) q[1];
sx q[1];
rz(-1.6769033) q[1];
sx q[1];
rz(1.5110663) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.481547) q[0];
sx q[0];
rz(-1.440795) q[0];
sx q[0];
rz(-1.2370716) q[0];
rz(1.3675646) q[2];
sx q[2];
rz(-2.2335238) q[2];
sx q[2];
rz(1.4858147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.293773) q[1];
sx q[1];
rz(-0.34983954) q[1];
sx q[1];
rz(1.185525) q[1];
rz(1.6919208) q[3];
sx q[3];
rz(-1.559404) q[3];
sx q[3];
rz(2.4081972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54603798) q[2];
sx q[2];
rz(-1.7067355) q[2];
sx q[2];
rz(-2.0175374) q[2];
rz(-1.3042287) q[3];
sx q[3];
rz(-0.29574695) q[3];
sx q[3];
rz(-0.058163253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27125636) q[0];
sx q[0];
rz(-0.79126343) q[0];
sx q[0];
rz(1.1668209) q[0];
rz(-1.600949) q[1];
sx q[1];
rz(-0.29778844) q[1];
sx q[1];
rz(1.3265532) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013915283) q[0];
sx q[0];
rz(-1.3260889) q[0];
sx q[0];
rz(0.6689594) q[0];
rz(0.4076414) q[2];
sx q[2];
rz(-1.2349814) q[2];
sx q[2];
rz(-1.4996355) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1455166) q[1];
sx q[1];
rz(-2.4156486) q[1];
sx q[1];
rz(-0.32246253) q[1];
rz(-1.79779) q[3];
sx q[3];
rz(-2.7309037) q[3];
sx q[3];
rz(-0.80055492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2549609) q[2];
sx q[2];
rz(-0.16356629) q[2];
sx q[2];
rz(1.6939885) q[2];
rz(3.0149031) q[3];
sx q[3];
rz(-1.4796175) q[3];
sx q[3];
rz(1.1160342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3348677) q[0];
sx q[0];
rz(-1.8341478) q[0];
sx q[0];
rz(-1.3679416) q[0];
rz(-1.5927636) q[1];
sx q[1];
rz(-0.82850414) q[1];
sx q[1];
rz(-2.9864476) q[1];
rz(-2.6947179) q[2];
sx q[2];
rz(-1.4021654) q[2];
sx q[2];
rz(-1.5709102) q[2];
rz(2.1550989) q[3];
sx q[3];
rz(-1.7135847) q[3];
sx q[3];
rz(1.8691487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
