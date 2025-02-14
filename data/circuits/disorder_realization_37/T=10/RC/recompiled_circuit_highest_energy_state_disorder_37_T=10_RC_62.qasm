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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(-2.6824644) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(-1.4067283) q[1];
sx q[1];
rz(-0.23174098) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2278175) q[0];
sx q[0];
rz(-1.639606) q[0];
sx q[0];
rz(2.7976996) q[0];
rz(-pi) q[1];
rz(-3.0475869) q[2];
sx q[2];
rz(-2.5931907) q[2];
sx q[2];
rz(1.3863871) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3331795) q[1];
sx q[1];
rz(-1.1929379) q[1];
sx q[1];
rz(2.4757132) q[1];
rz(-pi) q[2];
rz(-0.11756331) q[3];
sx q[3];
rz(-2.8010578) q[3];
sx q[3];
rz(-1.3887608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.090652466) q[2];
sx q[2];
rz(-0.53477627) q[2];
sx q[2];
rz(1.2789307) q[2];
rz(-0.88979641) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5559674) q[0];
sx q[0];
rz(-2.2323759) q[0];
sx q[0];
rz(0.92581785) q[0];
rz(-2.0766808) q[1];
sx q[1];
rz(-1.5771259) q[1];
sx q[1];
rz(3.056114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8877713) q[0];
sx q[0];
rz(-1.5740875) q[0];
sx q[0];
rz(2.7517767) q[0];
rz(0.34833724) q[2];
sx q[2];
rz(-2.1536963) q[2];
sx q[2];
rz(-2.9027241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8958455) q[1];
sx q[1];
rz(-1.3407882) q[1];
sx q[1];
rz(-1.0279015) q[1];
rz(-2.45473) q[3];
sx q[3];
rz(-1.9429038) q[3];
sx q[3];
rz(2.5290979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3257137) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(2.7549287) q[2];
rz(1.2169085) q[3];
sx q[3];
rz(-0.34438008) q[3];
sx q[3];
rz(-2.9215422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81095186) q[0];
sx q[0];
rz(-2.7702259) q[0];
sx q[0];
rz(2.6445342) q[0];
rz(2.0992384) q[1];
sx q[1];
rz(-2.1940239) q[1];
sx q[1];
rz(2.846948) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.153107) q[0];
sx q[0];
rz(-1.238958) q[0];
sx q[0];
rz(-1.4979304) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7435837) q[2];
sx q[2];
rz(-2.2026988) q[2];
sx q[2];
rz(1.2829219) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7859438) q[1];
sx q[1];
rz(-2.3501767) q[1];
sx q[1];
rz(-0.9407465) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14400836) q[3];
sx q[3];
rz(-0.97471182) q[3];
sx q[3];
rz(0.86897396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3543388) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(1.5617237) q[2];
rz(-2.0653557) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(-2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8266066) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(2.9122747) q[0];
rz(-1.6353105) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(3.07952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44767013) q[0];
sx q[0];
rz(-2.1248484) q[0];
sx q[0];
rz(2.5497132) q[0];
rz(-0.27917316) q[2];
sx q[2];
rz(-1.9491157) q[2];
sx q[2];
rz(-2.407848) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58846274) q[1];
sx q[1];
rz(-2.3847918) q[1];
sx q[1];
rz(-2.7314145) q[1];
x q[2];
rz(0.43020474) q[3];
sx q[3];
rz(-2.9211126) q[3];
sx q[3];
rz(-1.1663495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0486003) q[2];
sx q[2];
rz(-0.91602641) q[2];
sx q[2];
rz(1.830706) q[2];
rz(-0.038289573) q[3];
sx q[3];
rz(-1.7891276) q[3];
sx q[3];
rz(-2.766975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49486092) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(2.2485961) q[0];
rz(-1.0294186) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(0.095887862) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3678148) q[0];
sx q[0];
rz(-0.68499631) q[0];
sx q[0];
rz(-0.8789341) q[0];
x q[1];
rz(3.0711543) q[2];
sx q[2];
rz(-1.4820921) q[2];
sx q[2];
rz(-0.34102893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3600625) q[1];
sx q[1];
rz(-1.9340098) q[1];
sx q[1];
rz(-2.5468154) q[1];
rz(-2.5378102) q[3];
sx q[3];
rz(-0.60294916) q[3];
sx q[3];
rz(2.3100738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.15478495) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(0.42243877) q[3];
sx q[3];
rz(-1.1047624) q[3];
sx q[3];
rz(-0.5031684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7285889) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(0.11716209) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-2.2393176) q[1];
sx q[1];
rz(-1.0711627) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2881831) q[0];
sx q[0];
rz(-2.5700535) q[0];
sx q[0];
rz(-0.54067143) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6980459) q[2];
sx q[2];
rz(-2.5534592) q[2];
sx q[2];
rz(-0.085970446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92381645) q[1];
sx q[1];
rz(-0.72715302) q[1];
sx q[1];
rz(1.2283065) q[1];
rz(1.4070687) q[3];
sx q[3];
rz(-1.8843972) q[3];
sx q[3];
rz(-2.0255476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.57227197) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(0.039483698) q[2];
rz(-0.076400541) q[3];
sx q[3];
rz(-1.971784) q[3];
sx q[3];
rz(2.6881645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70306784) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(-1.9901265) q[1];
sx q[1];
rz(-1.6116424) q[1];
sx q[1];
rz(1.3075525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83043419) q[0];
sx q[0];
rz(-1.7896928) q[0];
sx q[0];
rz(-1.2375184) q[0];
x q[1];
rz(1.6436395) q[2];
sx q[2];
rz(-1.1955639) q[2];
sx q[2];
rz(-2.2323687) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7286018) q[1];
sx q[1];
rz(-2.6003693) q[1];
sx q[1];
rz(-1.9899998) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1767307) q[3];
sx q[3];
rz(-1.6814702) q[3];
sx q[3];
rz(2.0697274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82039708) q[2];
sx q[2];
rz(-0.68009192) q[2];
sx q[2];
rz(-0.55366984) q[2];
rz(0.6959483) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11947908) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(0.30690646) q[0];
rz(-1.8652929) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(-1.9006405) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0725482) q[0];
sx q[0];
rz(-2.6676151) q[0];
sx q[0];
rz(-1.9182253) q[0];
rz(-pi) q[1];
rz(-2.5503301) q[2];
sx q[2];
rz(-0.68702543) q[2];
sx q[2];
rz(-2.9369773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13682374) q[1];
sx q[1];
rz(-1.7589749) q[1];
sx q[1];
rz(-0.31463639) q[1];
x q[2];
rz(-1.1697024) q[3];
sx q[3];
rz(-2.1425284) q[3];
sx q[3];
rz(-2.9740262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3500195) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(-1.8512858) q[2];
rz(0.6811412) q[3];
sx q[3];
rz(-0.9001503) q[3];
sx q[3];
rz(-1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27555585) q[0];
sx q[0];
rz(-1.0365726) q[0];
sx q[0];
rz(1.9679029) q[0];
rz(-2.5545919) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(-3.0618844) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3143334) q[0];
sx q[0];
rz(-1.9406576) q[0];
sx q[0];
rz(-0.24263675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22144145) q[2];
sx q[2];
rz(-1.6802854) q[2];
sx q[2];
rz(-1.8851282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.6882872) q[1];
sx q[1];
rz(-2.2548179) q[1];
sx q[1];
rz(-0.75835336) q[1];
x q[2];
rz(-2.856273) q[3];
sx q[3];
rz(-0.82672182) q[3];
sx q[3];
rz(-0.43275012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54032636) q[2];
sx q[2];
rz(-1.2286295) q[2];
sx q[2];
rz(1.3339174) q[2];
rz(1.5506844) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(-0.18868119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.659336) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(1.2905066) q[0];
rz(0.036529649) q[1];
sx q[1];
rz(-1.1643658) q[1];
sx q[1];
rz(1.0443002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1636476) q[0];
sx q[0];
rz(-0.53102101) q[0];
sx q[0];
rz(0.95300302) q[0];
x q[1];
rz(-0.26728435) q[2];
sx q[2];
rz(-1.2517831) q[2];
sx q[2];
rz(2.6397982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86856213) q[1];
sx q[1];
rz(-0.84914637) q[1];
sx q[1];
rz(2.9356586) q[1];
rz(-2.812522) q[3];
sx q[3];
rz(-1.4621203) q[3];
sx q[3];
rz(1.4377126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2269939) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(-1.8168137) q[2];
rz(-0.75657183) q[3];
sx q[3];
rz(-2.2533267) q[3];
sx q[3];
rz(-2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.666438) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(1.4229763) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(-2.4074211) q[2];
sx q[2];
rz(-1.894507) q[2];
sx q[2];
rz(1.4121216) q[2];
rz(2.041009) q[3];
sx q[3];
rz(-2.1775424) q[3];
sx q[3];
rz(1.2829124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
