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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55104461) q[0];
sx q[0];
rz(-2.2529896) q[0];
sx q[0];
rz(2.5709573) q[0];
x q[1];
rz(-3.1093311) q[2];
sx q[2];
rz(-2.2164377) q[2];
sx q[2];
rz(2.1394859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1024485) q[1];
sx q[1];
rz(-1.7582201) q[1];
sx q[1];
rz(1.718912) q[1];
x q[2];
rz(1.5512439) q[3];
sx q[3];
rz(-2.3718155) q[3];
sx q[3];
rz(-2.584892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91780245) q[2];
sx q[2];
rz(-1.4897572) q[2];
sx q[2];
rz(-1.6415143) q[2];
rz(2.3598059) q[3];
sx q[3];
rz(-2.2512071) q[3];
sx q[3];
rz(-2.3894943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016945275) q[0];
sx q[0];
rz(-2.3638159) q[0];
sx q[0];
rz(0.97050226) q[0];
rz(-1.8151201) q[1];
sx q[1];
rz(-1.3423723) q[1];
sx q[1];
rz(2.2296947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058033179) q[0];
sx q[0];
rz(-0.61249295) q[0];
sx q[0];
rz(2.8508194) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7669589) q[2];
sx q[2];
rz(-1.4028143) q[2];
sx q[2];
rz(2.2998435) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7378547) q[1];
sx q[1];
rz(-2.4399839) q[1];
sx q[1];
rz(1.6340294) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5390057) q[3];
sx q[3];
rz(-0.45511757) q[3];
sx q[3];
rz(-1.2869175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9087002) q[2];
sx q[2];
rz(-2.0286655) q[2];
sx q[2];
rz(0.98928893) q[2];
rz(-2.7214637) q[3];
sx q[3];
rz(-2.1412854) q[3];
sx q[3];
rz(0.72107983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99783889) q[0];
sx q[0];
rz(-1.1622575) q[0];
sx q[0];
rz(-1.0379399) q[0];
rz(-0.85820091) q[1];
sx q[1];
rz(-0.52640262) q[1];
sx q[1];
rz(2.9720378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816681) q[0];
sx q[0];
rz(-1.5748275) q[0];
sx q[0];
rz(1.4272593) q[0];
rz(-pi) q[1];
rz(-2.8190024) q[2];
sx q[2];
rz(-1.1902404) q[2];
sx q[2];
rz(1.1104465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.922561) q[1];
sx q[1];
rz(-1.2777849) q[1];
sx q[1];
rz(2.0907563) q[1];
rz(-pi) q[2];
rz(-2.9870728) q[3];
sx q[3];
rz(-1.9294881) q[3];
sx q[3];
rz(-1.7254988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9618591) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(0.40599424) q[2];
rz(2.3640442) q[3];
sx q[3];
rz(-2.8113139) q[3];
sx q[3];
rz(-0.8146666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4985713) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(-0.92798573) q[0];
rz(0.058874933) q[1];
sx q[1];
rz(-1.4480271) q[1];
sx q[1];
rz(1.4307129) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74364036) q[0];
sx q[0];
rz(-2.7317606) q[0];
sx q[0];
rz(-1.7663971) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36305444) q[2];
sx q[2];
rz(-1.3842693) q[2];
sx q[2];
rz(-2.8721953) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8778363) q[1];
sx q[1];
rz(-2.0621513) q[1];
sx q[1];
rz(2.8884535) q[1];
rz(-2.1998422) q[3];
sx q[3];
rz(-1.3247196) q[3];
sx q[3];
rz(-1.6501282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4579939) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(1.9963473) q[2];
rz(-0.84732071) q[3];
sx q[3];
rz(-1.3342131) q[3];
sx q[3];
rz(-1.639521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(0.79967868) q[1];
sx q[1];
rz(-2.622602) q[1];
sx q[1];
rz(1.5934561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6131957) q[0];
sx q[0];
rz(-1.3269182) q[0];
sx q[0];
rz(0.24922483) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9351531) q[2];
sx q[2];
rz(-1.9597133) q[2];
sx q[2];
rz(-2.2909475) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1205129) q[1];
sx q[1];
rz(-1.9200293) q[1];
sx q[1];
rz(-1.9053188) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1372651) q[3];
sx q[3];
rz(-2.5133555) q[3];
sx q[3];
rz(-0.35997501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3758214) q[2];
sx q[2];
rz(-0.69362005) q[2];
sx q[2];
rz(3.0613464) q[2];
rz(-0.56096983) q[3];
sx q[3];
rz(-1.4813981) q[3];
sx q[3];
rz(-2.5188353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4866667) q[0];
sx q[0];
rz(-0.66513649) q[0];
sx q[0];
rz(1.2605865) q[0];
rz(-0.27443019) q[1];
sx q[1];
rz(-1.6703037) q[1];
sx q[1];
rz(-2.4226277) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2043122) q[0];
sx q[0];
rz(-2.215909) q[0];
sx q[0];
rz(-2.6995223) q[0];
rz(-pi) q[1];
rz(-1.5814797) q[2];
sx q[2];
rz(-2.3950393) q[2];
sx q[2];
rz(0.03554666) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.58108854) q[1];
sx q[1];
rz(-1.6916654) q[1];
sx q[1];
rz(0.64847364) q[1];
x q[2];
rz(1.0522862) q[3];
sx q[3];
rz(-0.42123367) q[3];
sx q[3];
rz(-3.0072834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4869953) q[2];
sx q[2];
rz(-1.2664814) q[2];
sx q[2];
rz(0.60301644) q[2];
rz(1.4740137) q[3];
sx q[3];
rz(-1.0242198) q[3];
sx q[3];
rz(-2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5176158) q[0];
sx q[0];
rz(-1.6083953) q[0];
sx q[0];
rz(-0.18727592) q[0];
rz(2.1693443) q[1];
sx q[1];
rz(-0.15521237) q[1];
sx q[1];
rz(0.42047277) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1218106) q[0];
sx q[0];
rz(-2.3276637) q[0];
sx q[0];
rz(-2.5518083) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37897972) q[2];
sx q[2];
rz(-1.0595269) q[2];
sx q[2];
rz(-0.015470964) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5198316) q[1];
sx q[1];
rz(-0.48096195) q[1];
sx q[1];
rz(2.2051808) q[1];
x q[2];
rz(-1.2240846) q[3];
sx q[3];
rz(-2.5390352) q[3];
sx q[3];
rz(-1.4977221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.04574) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(2.8908758) q[2];
rz(-1.2480674) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(-0.59984508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8641758) q[0];
sx q[0];
rz(-2.5499948) q[0];
sx q[0];
rz(-0.94171062) q[0];
rz(0.038453728) q[1];
sx q[1];
rz(-1.4310623) q[1];
sx q[1];
rz(-3.0252735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8909047) q[0];
sx q[0];
rz(-0.90849344) q[0];
sx q[0];
rz(0.57470365) q[0];
x q[1];
rz(-0.18949731) q[2];
sx q[2];
rz(-0.89841926) q[2];
sx q[2];
rz(-1.7799236) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6623936) q[1];
sx q[1];
rz(-2.8251007) q[1];
sx q[1];
rz(2.5280158) q[1];
rz(0.38920684) q[3];
sx q[3];
rz(-1.4688244) q[3];
sx q[3];
rz(-1.1679648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8470799) q[2];
sx q[2];
rz(-2.3307762) q[2];
sx q[2];
rz(-2.1152367) q[2];
rz(1.2096842) q[3];
sx q[3];
rz(-1.0299725) q[3];
sx q[3];
rz(-2.1799555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66861361) q[0];
sx q[0];
rz(-2.0675779) q[0];
sx q[0];
rz(-2.7630254) q[0];
rz(2.0319132) q[1];
sx q[1];
rz(-0.30866426) q[1];
sx q[1];
rz(3.1139156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6768764) q[0];
sx q[0];
rz(-0.72273556) q[0];
sx q[0];
rz(-2.5739772) q[0];
rz(-pi) q[1];
rz(0.86556025) q[2];
sx q[2];
rz(-2.8606374) q[2];
sx q[2];
rz(-2.11175) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3125449) q[1];
sx q[1];
rz(-1.2133737) q[1];
sx q[1];
rz(2.9446359) q[1];
x q[2];
rz(-2.1975193) q[3];
sx q[3];
rz(-0.51261307) q[3];
sx q[3];
rz(2.1533898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8370342) q[2];
sx q[2];
rz(-2.1874032) q[2];
sx q[2];
rz(0.41998106) q[2];
rz(-2.3000681) q[3];
sx q[3];
rz(-2.1244815) q[3];
sx q[3];
rz(-0.37874547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(-2.3751752) q[0];
sx q[0];
rz(-3.0247122) q[0];
sx q[0];
rz(-0.28512678) q[0];
rz(-2.9410703) q[1];
sx q[1];
rz(-1.9384117) q[1];
sx q[1];
rz(-0.83555046) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1211981) q[0];
sx q[0];
rz(-2.2403952) q[0];
sx q[0];
rz(0.59359896) q[0];
rz(-0.0086390583) q[2];
sx q[2];
rz(-2.4945179) q[2];
sx q[2];
rz(-1.6948014) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0802059) q[1];
sx q[1];
rz(-2.5701984) q[1];
sx q[1];
rz(-0.6462884) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5316758) q[3];
sx q[3];
rz(-1.3353383) q[3];
sx q[3];
rz(0.92604107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7150813) q[2];
sx q[2];
rz(-2.0501523) q[2];
sx q[2];
rz(2.4439028) q[2];
rz(0.52608025) q[3];
sx q[3];
rz(-2.3487921) q[3];
sx q[3];
rz(-1.5066159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0195011) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(-1.1493692) q[1];
sx q[1];
rz(-0.72723564) q[1];
sx q[1];
rz(-0.74692187) q[1];
rz(-2.7928945) q[2];
sx q[2];
rz(-0.92338466) q[2];
sx q[2];
rz(0.79104214) q[2];
rz(1.8676733) q[3];
sx q[3];
rz(-1.0424725) q[3];
sx q[3];
rz(-2.7155664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
