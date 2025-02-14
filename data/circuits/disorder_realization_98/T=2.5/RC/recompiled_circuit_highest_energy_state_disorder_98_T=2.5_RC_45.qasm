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
rz(-1.3774011) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55104461) q[0];
sx q[0];
rz(-2.2529896) q[0];
sx q[0];
rz(2.5709573) q[0];
rz(-1.6135869) q[2];
sx q[2];
rz(-2.4952609) q[2];
sx q[2];
rz(2.1930694) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5038576) q[1];
sx q[1];
rz(-1.4252932) q[1];
sx q[1];
rz(-0.18944959) q[1];
rz(-pi) q[2];
rz(0.018947424) q[3];
sx q[3];
rz(-0.80120443) q[3];
sx q[3];
rz(2.612118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2237902) q[2];
sx q[2];
rz(-1.6518355) q[2];
sx q[2];
rz(1.6415143) q[2];
rz(-2.3598059) q[3];
sx q[3];
rz(-2.2512071) q[3];
sx q[3];
rz(-0.75209832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1246474) q[0];
sx q[0];
rz(-0.77777672) q[0];
sx q[0];
rz(-0.97050226) q[0];
rz(1.8151201) q[1];
sx q[1];
rz(-1.7992203) q[1];
sx q[1];
rz(2.2296947) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3886984) q[0];
sx q[0];
rz(-1.4052183) q[0];
sx q[0];
rz(0.59247156) q[0];
rz(-0.17120338) q[2];
sx q[2];
rz(-1.7641626) q[2];
sx q[2];
rz(-0.69583508) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32102958) q[1];
sx q[1];
rz(-2.2707175) q[1];
sx q[1];
rz(3.0882443) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6025869) q[3];
sx q[3];
rz(-0.45511757) q[3];
sx q[3];
rz(1.2869175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9087002) q[2];
sx q[2];
rz(-2.0286655) q[2];
sx q[2];
rz(0.98928893) q[2];
rz(-0.42012897) q[3];
sx q[3];
rz(-1.0003072) q[3];
sx q[3];
rz(0.72107983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99783889) q[0];
sx q[0];
rz(-1.9793352) q[0];
sx q[0];
rz(1.0379399) q[0];
rz(0.85820091) q[1];
sx q[1];
rz(-2.61519) q[1];
sx q[1];
rz(-0.16955489) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9586048) q[0];
sx q[0];
rz(-2.9979994) q[0];
sx q[0];
rz(1.5989701) q[0];
rz(-pi) q[1];
rz(2.8190024) q[2];
sx q[2];
rz(-1.1902404) q[2];
sx q[2];
rz(-1.1104465) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3225786) q[1];
sx q[1];
rz(-2.5514304) q[1];
sx q[1];
rz(2.1165089) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15451984) q[3];
sx q[3];
rz(-1.9294881) q[3];
sx q[3];
rz(1.4160938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1797336) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(-2.7355984) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6430214) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(-2.2136069) q[0];
rz(0.058874933) q[1];
sx q[1];
rz(-1.6935655) q[1];
sx q[1];
rz(-1.4307129) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1851981) q[0];
sx q[0];
rz(-1.972359) q[0];
sx q[0];
rz(-3.0573581) q[0];
rz(0.48845993) q[2];
sx q[2];
rz(-0.40626981) q[2];
sx q[2];
rz(-0.84727188) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3763172) q[1];
sx q[1];
rz(-2.5936454) q[1];
sx q[1];
rz(2.0085232) q[1];
rz(-pi) q[2];
rz(-1.9742613) q[3];
sx q[3];
rz(-2.472252) q[3];
sx q[3];
rz(-2.739213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6835988) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(1.1452453) q[2];
rz(0.84732071) q[3];
sx q[3];
rz(-1.8073795) q[3];
sx q[3];
rz(1.5020717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1640846) q[0];
sx q[0];
rz(-1.9590398) q[0];
sx q[0];
rz(2.7567647) q[0];
rz(-0.79967868) q[1];
sx q[1];
rz(-0.5189907) q[1];
sx q[1];
rz(-1.5481366) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.858289) q[0];
sx q[0];
rz(-2.7946804) q[0];
sx q[0];
rz(-0.78979413) q[0];
rz(-pi) q[1];
rz(1.1069894) q[2];
sx q[2];
rz(-2.7037604) q[2];
sx q[2];
rz(-0.34604117) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1205129) q[1];
sx q[1];
rz(-1.9200293) q[1];
sx q[1];
rz(1.2362739) q[1];
x q[2];
rz(-2.1372651) q[3];
sx q[3];
rz(-0.62823717) q[3];
sx q[3];
rz(-2.7816176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3758214) q[2];
sx q[2];
rz(-0.69362005) q[2];
sx q[2];
rz(3.0613464) q[2];
rz(0.56096983) q[3];
sx q[3];
rz(-1.4813981) q[3];
sx q[3];
rz(-0.62275732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65492594) q[0];
sx q[0];
rz(-2.4764562) q[0];
sx q[0];
rz(-1.2605865) q[0];
rz(-0.27443019) q[1];
sx q[1];
rz(-1.6703037) q[1];
sx q[1];
rz(0.71896499) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9372805) q[0];
sx q[0];
rz(-2.215909) q[0];
sx q[0];
rz(-0.44207032) q[0];
rz(-2.3173213) q[2];
sx q[2];
rz(-1.5780515) q[2];
sx q[2];
rz(1.5985009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5605041) q[1];
sx q[1];
rz(-1.4499272) q[1];
sx q[1];
rz(-0.64847364) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9419233) q[3];
sx q[3];
rz(-1.3667445) q[3];
sx q[3];
rz(0.95641092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4869953) q[2];
sx q[2];
rz(-1.8751112) q[2];
sx q[2];
rz(0.60301644) q[2];
rz(-1.4740137) q[3];
sx q[3];
rz(-1.0242198) q[3];
sx q[3];
rz(2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5176158) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(0.18727592) q[0];
rz(-2.1693443) q[1];
sx q[1];
rz(-0.15521237) q[1];
sx q[1];
rz(2.7211199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.792345) q[0];
sx q[0];
rz(-2.2195243) q[0];
sx q[0];
rz(-2.1029841) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0275317) q[2];
sx q[2];
rz(-1.8993371) q[2];
sx q[2];
rz(-1.3628886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9686465) q[1];
sx q[1];
rz(-1.1889646) q[1];
sx q[1];
rz(2.8416425) q[1];
rz(-pi) q[2];
rz(2.1450348) q[3];
sx q[3];
rz(-1.3769994) q[3];
sx q[3];
rz(2.7791948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.095852701) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(0.25071684) q[2];
rz(1.2480674) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(0.59984508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8641758) q[0];
sx q[0];
rz(-2.5499948) q[0];
sx q[0];
rz(-0.94171062) q[0];
rz(3.1031389) q[1];
sx q[1];
rz(-1.4310623) q[1];
sx q[1];
rz(3.0252735) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021892) q[0];
sx q[0];
rz(-0.84745126) q[0];
sx q[0];
rz(2.1795616) q[0];
rz(0.18949731) q[2];
sx q[2];
rz(-2.2431734) q[2];
sx q[2];
rz(1.3616691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1585576) q[1];
sx q[1];
rz(-1.3135034) q[1];
sx q[1];
rz(1.3844116) q[1];
rz(-pi) q[2];
rz(2.8781901) q[3];
sx q[3];
rz(-0.40168328) q[3];
sx q[3];
rz(-0.64611891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2945127) q[2];
sx q[2];
rz(-2.3307762) q[2];
sx q[2];
rz(-1.026356) q[2];
rz(-1.2096842) q[3];
sx q[3];
rz(-1.0299725) q[3];
sx q[3];
rz(-0.9616372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.472979) q[0];
sx q[0];
rz(-1.0740148) q[0];
sx q[0];
rz(0.37856722) q[0];
rz(2.0319132) q[1];
sx q[1];
rz(-0.30866426) q[1];
sx q[1];
rz(-0.027677061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4647163) q[0];
sx q[0];
rz(-2.4188571) q[0];
sx q[0];
rz(2.5739772) q[0];
rz(-pi) q[1];
rz(-0.1849298) q[2];
sx q[2];
rz(-1.7835296) q[2];
sx q[2];
rz(-0.30480584) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82904774) q[1];
sx q[1];
rz(-1.9282189) q[1];
sx q[1];
rz(2.9446359) q[1];
rz(-pi) q[2];
rz(-0.31881551) q[3];
sx q[3];
rz(-1.1622815) q[3];
sx q[3];
rz(-2.8466948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8370342) q[2];
sx q[2];
rz(-2.1874032) q[2];
sx q[2];
rz(2.7216116) q[2];
rz(2.3000681) q[3];
sx q[3];
rz(-2.1244815) q[3];
sx q[3];
rz(0.37874547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3751752) q[0];
sx q[0];
rz(-3.0247122) q[0];
sx q[0];
rz(0.28512678) q[0];
rz(0.20052234) q[1];
sx q[1];
rz(-1.2031809) q[1];
sx q[1];
rz(0.83555046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1932521) q[0];
sx q[0];
rz(-0.86328816) q[0];
sx q[0];
rz(2.18594) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1329536) q[2];
sx q[2];
rz(-0.64707478) q[2];
sx q[2];
rz(1.4467913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0802059) q[1];
sx q[1];
rz(-2.5701984) q[1];
sx q[1];
rz(-0.6462884) q[1];
rz(-1.6099168) q[3];
sx q[3];
rz(-1.8062544) q[3];
sx q[3];
rz(-0.92604107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7150813) q[2];
sx q[2];
rz(-1.0914404) q[2];
sx q[2];
rz(-2.4439028) q[2];
rz(2.6155124) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(-1.5066159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0195011) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(-1.9922235) q[1];
sx q[1];
rz(-2.414357) q[1];
sx q[1];
rz(2.3946708) q[1];
rz(0.89328881) q[2];
sx q[2];
rz(-1.2947686) q[2];
sx q[2];
rz(-0.99560621) q[2];
rz(1.2739194) q[3];
sx q[3];
rz(-2.0991201) q[3];
sx q[3];
rz(0.42602628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
