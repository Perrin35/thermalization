OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(0.36261121) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(-0.4904823) q[1];
sx q[1];
rz(2.7999556) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5255614) q[0];
sx q[0];
rz(-1.4481059) q[0];
sx q[0];
rz(-2.7659155) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5506637) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(0.17609827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4075549) q[1];
sx q[1];
rz(-1.8681548) q[1];
sx q[1];
rz(1.0337017) q[1];
x q[2];
rz(2.0658675) q[3];
sx q[3];
rz(-2.6945811) q[3];
sx q[3];
rz(1.3017553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59387702) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(2.1851052) q[2];
rz(-0.18125136) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(-1.6199002) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.99615) q[0];
rz(-1.068813) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(2.4157445) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23913357) q[0];
sx q[0];
rz(-2.4665138) q[0];
sx q[0];
rz(2.549987) q[0];
x q[1];
rz(-2.5088596) q[2];
sx q[2];
rz(-0.69721141) q[2];
sx q[2];
rz(2.7998507) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5686381) q[1];
sx q[1];
rz(-1.2619839) q[1];
sx q[1];
rz(-3.1152578) q[1];
rz(-2.5190982) q[3];
sx q[3];
rz(-1.7433634) q[3];
sx q[3];
rz(1.7164001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(0.99728161) q[2];
rz(-1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(-2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.41855758) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(-1.3035125) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(-2.2361141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83429503) q[0];
sx q[0];
rz(-0.30954888) q[0];
sx q[0];
rz(0.66972591) q[0];
rz(-pi) q[1];
rz(0.4586787) q[2];
sx q[2];
rz(-0.80692601) q[2];
sx q[2];
rz(-0.57070953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14721522) q[1];
sx q[1];
rz(-0.71502393) q[1];
sx q[1];
rz(1.6836402) q[1];
rz(-pi) q[2];
rz(0.043025322) q[3];
sx q[3];
rz(-2.283841) q[3];
sx q[3];
rz(-0.42792861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3811615) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(-0.9872438) q[2];
rz(0.045771249) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10953294) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(-3.0155244) q[0];
rz(-0.18445045) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.9741845) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4520893) q[0];
sx q[0];
rz(-0.79012442) q[0];
sx q[0];
rz(2.8566314) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7502968) q[2];
sx q[2];
rz(-2.0677535) q[2];
sx q[2];
rz(-2.3617982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0427525) q[1];
sx q[1];
rz(-0.71430695) q[1];
sx q[1];
rz(0.40978281) q[1];
rz(0.64658029) q[3];
sx q[3];
rz(-2.1255891) q[3];
sx q[3];
rz(0.085689714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2771153) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(-0.63151044) q[2];
rz(0.19550066) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(-0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1642078) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-2.6389129) q[0];
rz(-1.1133105) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(0.46708333) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53411667) q[0];
sx q[0];
rz(-2.7152938) q[0];
sx q[0];
rz(-2.0712453) q[0];
rz(-pi) q[1];
rz(2.5555243) q[2];
sx q[2];
rz(-1.7826826) q[2];
sx q[2];
rz(-0.52473247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0847229) q[1];
sx q[1];
rz(-1.1493756) q[1];
sx q[1];
rz(-0.23878581) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15848666) q[3];
sx q[3];
rz(-2.2500258) q[3];
sx q[3];
rz(-2.2331626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46665329) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(2.6494027) q[2];
rz(-0.28218937) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(-1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(1.4280691) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(1.496398) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8498358) q[0];
sx q[0];
rz(-1.0514326) q[0];
sx q[0];
rz(-2.915328) q[0];
rz(2.9005269) q[2];
sx q[2];
rz(-2.0965555) q[2];
sx q[2];
rz(0.17503967) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8138258) q[1];
sx q[1];
rz(-1.5944949) q[1];
sx q[1];
rz(1.2162672) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5580642) q[3];
sx q[3];
rz(-0.85789645) q[3];
sx q[3];
rz(1.7604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66203403) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(0.39624828) q[2];
rz(2.5475492) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783766) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(2.0492045) q[0];
rz(-1.7842402) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(3.1299652) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8788293) q[0];
sx q[0];
rz(-1.8352404) q[0];
sx q[0];
rz(-2.7468365) q[0];
x q[1];
rz(0.02545698) q[2];
sx q[2];
rz(-1.5989132) q[2];
sx q[2];
rz(-2.1643929) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7736518) q[1];
sx q[1];
rz(-2.5950187) q[1];
sx q[1];
rz(1.6163338) q[1];
rz(-1.2467975) q[3];
sx q[3];
rz(-0.9231336) q[3];
sx q[3];
rz(0.049829359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6860883) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(-0.28798506) q[2];
rz(1.9912432) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(-2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(-1.9294552) q[0];
rz(2.7638226) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15570116) q[0];
sx q[0];
rz(-0.37030664) q[0];
sx q[0];
rz(0.62472384) q[0];
rz(3.0817501) q[2];
sx q[2];
rz(-1.7027731) q[2];
sx q[2];
rz(-3.0795385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89477506) q[1];
sx q[1];
rz(-2.2404039) q[1];
sx q[1];
rz(-1.0424022) q[1];
rz(1.0687374) q[3];
sx q[3];
rz(-1.6747253) q[3];
sx q[3];
rz(-0.41307377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3528072) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(-1.8898233) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733474) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(2.9633203) q[0];
rz(3.0687304) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3423791) q[0];
sx q[0];
rz(-1.8633435) q[0];
sx q[0];
rz(0.53036687) q[0];
x q[1];
rz(-0.91782848) q[2];
sx q[2];
rz(-0.73790109) q[2];
sx q[2];
rz(-0.81673056) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8930248) q[1];
sx q[1];
rz(-2.5556459) q[1];
sx q[1];
rz(-0.47793169) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9411373) q[3];
sx q[3];
rz(-0.59560532) q[3];
sx q[3];
rz(-2.5411118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4420085) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(-1.0317831) q[2];
rz(-1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(-0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853302) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(1.306698) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(2.0526989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13578829) q[0];
sx q[0];
rz(-3.0514768) q[0];
sx q[0];
rz(2.3401005) q[0];
rz(-0.18386545) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(-2.0053787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2838193) q[1];
sx q[1];
rz(-2.3771264) q[1];
sx q[1];
rz(-2.43999) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18967929) q[3];
sx q[3];
rz(-1.2343725) q[3];
sx q[3];
rz(1.3996901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(-1.5169253) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(0.91296997) q[2];
sx q[2];
rz(-1.8987449) q[2];
sx q[2];
rz(1.434025) q[2];
rz(-1.2354479) q[3];
sx q[3];
rz(-1.609758) q[3];
sx q[3];
rz(2.5555425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
