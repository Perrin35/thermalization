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
rz(2.1762525) q[0];
rz(0.020429285) q[1];
sx q[1];
rz(-1.483622) q[1];
sx q[1];
rz(-1.3774011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7962153) q[0];
sx q[0];
rz(-0.85897972) q[0];
sx q[0];
rz(-2.1576361) q[0];
x q[1];
rz(0.032261511) q[2];
sx q[2];
rz(-0.92515495) q[2];
sx q[2];
rz(-2.1394859) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6377351) q[1];
sx q[1];
rz(-1.4252932) q[1];
sx q[1];
rz(2.9521431) q[1];
x q[2];
rz(-0.018947424) q[3];
sx q[3];
rz(-2.3403882) q[3];
sx q[3];
rz(-0.5294746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2237902) q[2];
sx q[2];
rz(-1.4897572) q[2];
sx q[2];
rz(-1.6415143) q[2];
rz(0.7817868) q[3];
sx q[3];
rz(-2.2512071) q[3];
sx q[3];
rz(-0.75209832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1246474) q[0];
sx q[0];
rz(-2.3638159) q[0];
sx q[0];
rz(-0.97050226) q[0];
rz(1.3264725) q[1];
sx q[1];
rz(-1.3423723) q[1];
sx q[1];
rz(2.2296947) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7528943) q[0];
sx q[0];
rz(-1.4052183) q[0];
sx q[0];
rz(0.59247156) q[0];
x q[1];
rz(0.85477306) q[2];
sx q[2];
rz(-0.25755421) q[2];
sx q[2];
rz(1.7130898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8205631) q[1];
sx q[1];
rz(-2.2707175) q[1];
sx q[1];
rz(0.053348347) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5390057) q[3];
sx q[3];
rz(-0.45511757) q[3];
sx q[3];
rz(-1.8546752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9087002) q[2];
sx q[2];
rz(-1.1129271) q[2];
sx q[2];
rz(2.1523037) q[2];
rz(2.7214637) q[3];
sx q[3];
rz(-2.1412854) q[3];
sx q[3];
rz(2.4205128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1437538) q[0];
sx q[0];
rz(-1.9793352) q[0];
sx q[0];
rz(2.1036527) q[0];
rz(-0.85820091) q[1];
sx q[1];
rz(-2.61519) q[1];
sx q[1];
rz(-2.9720378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9301382) q[0];
sx q[0];
rz(-1.4272604) q[0];
sx q[0];
rz(-3.1375196) q[0];
x q[1];
rz(2.8190024) q[2];
sx q[2];
rz(-1.9513522) q[2];
sx q[2];
rz(-2.0311462) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3225786) q[1];
sx q[1];
rz(-2.5514304) q[1];
sx q[1];
rz(1.0250837) q[1];
x q[2];
rz(1.9334458) q[3];
sx q[3];
rz(-1.4261822) q[3];
sx q[3];
rz(-0.1000769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1797336) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(-0.40599424) q[2];
rz(-0.77754846) q[3];
sx q[3];
rz(-2.8113139) q[3];
sx q[3];
rz(-0.8146666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6430214) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(-0.92798573) q[0];
rz(-0.058874933) q[1];
sx q[1];
rz(-1.4480271) q[1];
sx q[1];
rz(-1.4307129) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95639455) q[0];
sx q[0];
rz(-1.972359) q[0];
sx q[0];
rz(3.0573581) q[0];
rz(-pi) q[1];
rz(-2.6531327) q[2];
sx q[2];
rz(-2.7353228) q[2];
sx q[2];
rz(0.84727188) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7652755) q[1];
sx q[1];
rz(-2.5936454) q[1];
sx q[1];
rz(-2.0085232) q[1];
rz(-2.840418) q[3];
sx q[3];
rz(-2.1780925) q[3];
sx q[3];
rz(-0.09610304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6835988) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(1.9963473) q[2];
rz(-2.2942719) q[3];
sx q[3];
rz(-1.8073795) q[3];
sx q[3];
rz(-1.639521) q[3];
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
rz(0.97750807) q[0];
sx q[0];
rz(-1.1825528) q[0];
sx q[0];
rz(-0.384828) q[0];
rz(2.341914) q[1];
sx q[1];
rz(-0.5189907) q[1];
sx q[1];
rz(1.5934561) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.528397) q[0];
sx q[0];
rz(-1.8146744) q[0];
sx q[0];
rz(2.8923678) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9673011) q[2];
sx q[2];
rz(-1.3799715) q[2];
sx q[2];
rz(2.3421974) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0210798) q[1];
sx q[1];
rz(-1.9200293) q[1];
sx q[1];
rz(-1.9053188) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1206854) q[3];
sx q[3];
rz(-1.249915) q[3];
sx q[3];
rz(2.4060017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7657713) q[2];
sx q[2];
rz(-0.69362005) q[2];
sx q[2];
rz(-3.0613464) q[2];
rz(-2.5806228) q[3];
sx q[3];
rz(-1.4813981) q[3];
sx q[3];
rz(2.5188353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65492594) q[0];
sx q[0];
rz(-2.4764562) q[0];
sx q[0];
rz(1.8810062) q[0];
rz(-2.8671625) q[1];
sx q[1];
rz(-1.471289) q[1];
sx q[1];
rz(-2.4226277) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9372805) q[0];
sx q[0];
rz(-0.92568362) q[0];
sx q[0];
rz(0.44207032) q[0];
x q[1];
rz(1.5601129) q[2];
sx q[2];
rz(-2.3950393) q[2];
sx q[2];
rz(-3.106046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9940954) q[1];
sx q[1];
rz(-0.65804099) q[1];
sx q[1];
rz(-2.9431353) q[1];
rz(2.9230887) q[3];
sx q[3];
rz(-1.2077304) q[3];
sx q[3];
rz(2.4485112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6545973) q[2];
sx q[2];
rz(-1.8751112) q[2];
sx q[2];
rz(2.5385762) q[2];
rz(-1.4740137) q[3];
sx q[3];
rz(-2.1173729) q[3];
sx q[3];
rz(-2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6239768) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(2.9543167) q[0];
rz(2.1693443) q[1];
sx q[1];
rz(-2.9863803) q[1];
sx q[1];
rz(2.7211199) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.792345) q[0];
sx q[0];
rz(-0.92206832) q[0];
sx q[0];
rz(1.0386085) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37897972) q[2];
sx q[2];
rz(-2.0820658) q[2];
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
rz(-pi/2) q[0];
rz(-1.6290038) q[1];
sx q[1];
rz(-1.8485475) q[1];
sx q[1];
rz(-1.968683) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1450348) q[3];
sx q[3];
rz(-1.3769994) q[3];
sx q[3];
rz(0.36239788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.04574) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(-2.8908758) q[2];
rz(-1.2480674) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(2.5417476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8641758) q[0];
sx q[0];
rz(-2.5499948) q[0];
sx q[0];
rz(-0.94171062) q[0];
rz(-3.1031389) q[1];
sx q[1];
rz(-1.4310623) q[1];
sx q[1];
rz(-3.0252735) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4394035) q[0];
sx q[0];
rz(-0.84745126) q[0];
sx q[0];
rz(0.96203104) q[0];
rz(1.8031249) q[2];
sx q[2];
rz(-0.69455494) q[2];
sx q[2];
rz(1.6603927) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47919905) q[1];
sx q[1];
rz(-2.8251007) q[1];
sx q[1];
rz(-2.5280158) q[1];
rz(-0.26340254) q[3];
sx q[3];
rz(-0.40168328) q[3];
sx q[3];
rz(2.4954737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2945127) q[2];
sx q[2];
rz(-2.3307762) q[2];
sx q[2];
rz(2.1152367) q[2];
rz(1.9319084) q[3];
sx q[3];
rz(-1.0299725) q[3];
sx q[3];
rz(2.1799555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.472979) q[0];
sx q[0];
rz(-1.0740148) q[0];
sx q[0];
rz(0.37856722) q[0];
rz(-1.1096795) q[1];
sx q[1];
rz(-2.8329284) q[1];
sx q[1];
rz(0.027677061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8016114) q[0];
sx q[0];
rz(-1.2072354) q[0];
sx q[0];
rz(2.5021863) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1849298) q[2];
sx q[2];
rz(-1.3580631) q[2];
sx q[2];
rz(0.30480584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31071404) q[1];
sx q[1];
rz(-2.7355477) q[1];
sx q[1];
rz(-2.0534404) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94407336) q[3];
sx q[3];
rz(-0.51261307) q[3];
sx q[3];
rz(2.1533898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3045584) q[2];
sx q[2];
rz(-2.1874032) q[2];
sx q[2];
rz(0.41998106) q[2];
rz(-2.3000681) q[3];
sx q[3];
rz(-1.0171112) q[3];
sx q[3];
rz(-2.7628472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-2.8564659) q[0];
rz(0.20052234) q[1];
sx q[1];
rz(-1.2031809) q[1];
sx q[1];
rz(-2.3060422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0203945) q[0];
sx q[0];
rz(-2.2403952) q[0];
sx q[0];
rz(-2.5479937) q[0];
rz(-pi) q[1];
rz(-0.64705683) q[2];
sx q[2];
rz(-1.5655883) q[2];
sx q[2];
rz(-3.0106949) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7924462) q[1];
sx q[1];
rz(-1.1243774) q[1];
sx q[1];
rz(-1.940215) q[1];
rz(2.9799913) q[3];
sx q[3];
rz(-0.23862632) q[3];
sx q[3];
rz(2.0493226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7150813) q[2];
sx q[2];
rz(-1.0914404) q[2];
sx q[2];
rz(-0.69768989) q[2];
rz(-0.52608025) q[3];
sx q[3];
rz(-2.3487921) q[3];
sx q[3];
rz(-1.6349767) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1220916) q[0];
sx q[0];
rz(-0.87845907) q[0];
sx q[0];
rz(-0.39977951) q[0];
rz(1.9922235) q[1];
sx q[1];
rz(-0.72723564) q[1];
sx q[1];
rz(-0.74692187) q[1];
rz(1.1463852) q[2];
sx q[2];
rz(-0.72327264) q[2];
sx q[2];
rz(-2.8930152) q[2];
rz(-2.6769842) q[3];
sx q[3];
rz(-2.5425442) q[3];
sx q[3];
rz(3.0221593) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
