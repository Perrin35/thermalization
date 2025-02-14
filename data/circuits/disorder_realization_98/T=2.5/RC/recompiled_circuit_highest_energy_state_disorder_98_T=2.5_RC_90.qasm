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
rz(3.0124445) q[0];
sx q[0];
rz(-1.4718066) q[0];
sx q[0];
rz(-0.9653402) q[0];
rz(0.020429285) q[1];
sx q[1];
rz(-1.483622) q[1];
sx q[1];
rz(1.7641915) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55104461) q[0];
sx q[0];
rz(-2.2529896) q[0];
sx q[0];
rz(2.5709573) q[0];
rz(-0.92490478) q[2];
sx q[2];
rz(-1.5450302) q[2];
sx q[2];
rz(0.58810593) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6377351) q[1];
sx q[1];
rz(-1.4252932) q[1];
sx q[1];
rz(2.9521431) q[1];
rz(-pi) q[2];
rz(0.018947424) q[3];
sx q[3];
rz(-2.3403882) q[3];
sx q[3];
rz(-2.612118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2237902) q[2];
sx q[2];
rz(-1.4897572) q[2];
sx q[2];
rz(-1.5000783) q[2];
rz(2.3598059) q[3];
sx q[3];
rz(-0.89038554) q[3];
sx q[3];
rz(-0.75209832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1246474) q[0];
sx q[0];
rz(-2.3638159) q[0];
sx q[0];
rz(2.1710904) q[0];
rz(1.3264725) q[1];
sx q[1];
rz(-1.3423723) q[1];
sx q[1];
rz(-0.91189799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058033179) q[0];
sx q[0];
rz(-2.5290997) q[0];
sx q[0];
rz(2.8508194) q[0];
x q[1];
rz(0.17120338) q[2];
sx q[2];
rz(-1.3774301) q[2];
sx q[2];
rz(2.4457576) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7378547) q[1];
sx q[1];
rz(-0.7016088) q[1];
sx q[1];
rz(-1.5075633) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6025869) q[3];
sx q[3];
rz(-0.45511757) q[3];
sx q[3];
rz(-1.2869175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2328925) q[2];
sx q[2];
rz(-2.0286655) q[2];
sx q[2];
rz(2.1523037) q[2];
rz(-2.7214637) q[3];
sx q[3];
rz(-2.1412854) q[3];
sx q[3];
rz(0.72107983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99783889) q[0];
sx q[0];
rz(-1.9793352) q[0];
sx q[0];
rz(-2.1036527) q[0];
rz(2.2833917) q[1];
sx q[1];
rz(-0.52640262) q[1];
sx q[1];
rz(-0.16955489) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9586048) q[0];
sx q[0];
rz(-0.14359328) q[0];
sx q[0];
rz(1.5426226) q[0];
rz(-pi) q[1];
rz(-0.900678) q[2];
sx q[2];
rz(-2.6478516) q[2];
sx q[2];
rz(0.37767071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3225786) q[1];
sx q[1];
rz(-0.59016229) q[1];
sx q[1];
rz(-1.0250837) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1812594) q[3];
sx q[3];
rz(-0.3892322) q[3];
sx q[3];
rz(-1.3078944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1797336) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(-0.40599424) q[2];
rz(2.3640442) q[3];
sx q[3];
rz(-2.8113139) q[3];
sx q[3];
rz(2.3269261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6430214) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(-2.2136069) q[0];
rz(-0.058874933) q[1];
sx q[1];
rz(-1.6935655) q[1];
sx q[1];
rz(1.4307129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64739156) q[0];
sx q[0];
rz(-1.6483161) q[0];
sx q[0];
rz(-1.1679542) q[0];
rz(-pi) q[1];
rz(0.36305444) q[2];
sx q[2];
rz(-1.7573234) q[2];
sx q[2];
rz(2.8721953) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9560066) q[1];
sx q[1];
rz(-1.7934467) q[1];
sx q[1];
rz(-1.0658479) q[1];
rz(0.30117463) q[3];
sx q[3];
rz(-0.96350017) q[3];
sx q[3];
rz(0.09610304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4579939) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(1.1452453) q[2];
rz(-0.84732071) q[3];
sx q[3];
rz(-1.8073795) q[3];
sx q[3];
rz(-1.5020717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1640846) q[0];
sx q[0];
rz(-1.1825528) q[0];
sx q[0];
rz(-2.7567647) q[0];
rz(0.79967868) q[1];
sx q[1];
rz(-0.5189907) q[1];
sx q[1];
rz(-1.5934561) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.528397) q[0];
sx q[0];
rz(-1.8146744) q[0];
sx q[0];
rz(-0.24922483) q[0];
rz(1.9673011) q[2];
sx q[2];
rz(-1.7616211) q[2];
sx q[2];
rz(2.3421974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0210798) q[1];
sx q[1];
rz(-1.9200293) q[1];
sx q[1];
rz(1.2362739) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0209072) q[3];
sx q[3];
rz(-1.8916777) q[3];
sx q[3];
rz(-2.4060017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3758214) q[2];
sx q[2];
rz(-2.4479726) q[2];
sx q[2];
rz(-3.0613464) q[2];
rz(2.5806228) q[3];
sx q[3];
rz(-1.4813981) q[3];
sx q[3];
rz(-2.5188353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65492594) q[0];
sx q[0];
rz(-0.66513649) q[0];
sx q[0];
rz(-1.2605865) q[0];
rz(0.27443019) q[1];
sx q[1];
rz(-1.6703037) q[1];
sx q[1];
rz(-0.71896499) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711766) q[0];
sx q[0];
rz(-2.377802) q[0];
sx q[0];
rz(1.0538573) q[0];
rz(-pi) q[1];
rz(-1.5814797) q[2];
sx q[2];
rz(-2.3950393) q[2];
sx q[2];
rz(0.03554666) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58108854) q[1];
sx q[1];
rz(-1.6916654) q[1];
sx q[1];
rz(2.493119) q[1];
rz(1.1996693) q[3];
sx q[3];
rz(-1.7748482) q[3];
sx q[3];
rz(0.95641092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4869953) q[2];
sx q[2];
rz(-1.8751112) q[2];
sx q[2];
rz(2.5385762) q[2];
rz(-1.4740137) q[3];
sx q[3];
rz(-1.0242198) q[3];
sx q[3];
rz(2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6239768) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(0.18727592) q[0];
rz(-2.1693443) q[1];
sx q[1];
rz(-2.9863803) q[1];
sx q[1];
rz(0.42047277) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12026006) q[0];
sx q[0];
rz(-1.1545404) q[0];
sx q[0];
rz(0.72159213) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7626129) q[2];
sx q[2];
rz(-1.0595269) q[2];
sx q[2];
rz(-3.1261217) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.17294614) q[1];
sx q[1];
rz(-1.9526281) q[1];
sx q[1];
rz(2.8416425) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22962946) q[3];
sx q[3];
rz(-2.1329704) q[3];
sx q[3];
rz(-2.0571518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.04574) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(-0.25071684) q[2];
rz(-1.8935253) q[3];
sx q[3];
rz(-1.3919132) q[3];
sx q[3];
rz(2.5417476) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2774169) q[0];
sx q[0];
rz(-2.5499948) q[0];
sx q[0];
rz(-2.199882) q[0];
rz(-0.038453728) q[1];
sx q[1];
rz(-1.7105303) q[1];
sx q[1];
rz(0.11631913) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8909047) q[0];
sx q[0];
rz(-2.2330992) q[0];
sx q[0];
rz(0.57470365) q[0];
x q[1];
rz(1.3384678) q[2];
sx q[2];
rz(-2.4470377) q[2];
sx q[2];
rz(1.6603927) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1585576) q[1];
sx q[1];
rz(-1.8280892) q[1];
sx q[1];
rz(1.757181) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6809471) q[3];
sx q[3];
rz(-1.9578736) q[3];
sx q[3];
rz(-0.36110669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2945127) q[2];
sx q[2];
rz(-0.81081644) q[2];
sx q[2];
rz(-1.026356) q[2];
rz(-1.2096842) q[3];
sx q[3];
rz(-1.0299725) q[3];
sx q[3];
rz(2.1799555) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.472979) q[0];
sx q[0];
rz(-2.0675779) q[0];
sx q[0];
rz(-2.7630254) q[0];
rz(-2.0319132) q[1];
sx q[1];
rz(-0.30866426) q[1];
sx q[1];
rz(0.027677061) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6768764) q[0];
sx q[0];
rz(-2.4188571) q[0];
sx q[0];
rz(-2.5739772) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1849298) q[2];
sx q[2];
rz(-1.7835296) q[2];
sx q[2];
rz(-2.8367868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81144801) q[1];
sx q[1];
rz(-1.7551577) q[1];
sx q[1];
rz(1.2069279) q[1];
rz(-pi) q[2];
rz(2.8227771) q[3];
sx q[3];
rz(-1.1622815) q[3];
sx q[3];
rz(0.2948979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3045584) q[2];
sx q[2];
rz(-0.95418945) q[2];
sx q[2];
rz(-0.41998106) q[2];
rz(-0.84152451) q[3];
sx q[3];
rz(-1.0171112) q[3];
sx q[3];
rz(2.7628472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751752) q[0];
sx q[0];
rz(-0.11688047) q[0];
sx q[0];
rz(0.28512678) q[0];
rz(2.9410703) q[1];
sx q[1];
rz(-1.2031809) q[1];
sx q[1];
rz(2.3060422) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1211981) q[0];
sx q[0];
rz(-0.9011974) q[0];
sx q[0];
rz(-2.5479937) q[0];
x q[1];
rz(1.5642688) q[2];
sx q[2];
rz(-2.2178429) q[2];
sx q[2];
rz(-1.4359635) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0802059) q[1];
sx q[1];
rz(-2.5701984) q[1];
sx q[1];
rz(2.4953043) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9799913) q[3];
sx q[3];
rz(-2.9029663) q[3];
sx q[3];
rz(2.0493226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7150813) q[2];
sx q[2];
rz(-1.0914404) q[2];
sx q[2];
rz(0.69768989) q[2];
rz(2.6155124) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(1.6349767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1220916) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(1.9922235) q[1];
sx q[1];
rz(-0.72723564) q[1];
sx q[1];
rz(-0.74692187) q[1];
rz(2.2483038) q[2];
sx q[2];
rz(-1.8468241) q[2];
sx q[2];
rz(2.1459864) q[2];
rz(2.5935843) q[3];
sx q[3];
rz(-1.3153793) q[3];
sx q[3];
rz(1.8438189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
