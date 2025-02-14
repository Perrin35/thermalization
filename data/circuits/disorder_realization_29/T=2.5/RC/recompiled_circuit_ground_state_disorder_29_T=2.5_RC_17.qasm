OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.27958265) q[0];
sx q[0];
rz(-2.6065338) q[0];
sx q[0];
rz(-0.23393272) q[0];
rz(-2.2995931) q[1];
sx q[1];
rz(-1.41398) q[1];
sx q[1];
rz(1.5185305) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5748225) q[0];
sx q[0];
rz(-1.8810417) q[0];
sx q[0];
rz(-0.53939087) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0064924) q[2];
sx q[2];
rz(-0.47347906) q[2];
sx q[2];
rz(-1.4896637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.17659345) q[1];
sx q[1];
rz(-1.0142583) q[1];
sx q[1];
rz(0.92265572) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2330448) q[3];
sx q[3];
rz(-1.1881785) q[3];
sx q[3];
rz(0.97272129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6003549) q[2];
sx q[2];
rz(-1.7011832) q[2];
sx q[2];
rz(2.3360628) q[2];
rz(2.7496998) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(-2.384757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517752) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(-2.7031194) q[0];
rz(-1.27502) q[1];
sx q[1];
rz(-1.6908815) q[1];
sx q[1];
rz(3.0335887) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1718801) q[0];
sx q[0];
rz(-1.4599271) q[0];
sx q[0];
rz(3.1165439) q[0];
rz(-pi) q[1];
rz(-0.63814068) q[2];
sx q[2];
rz(-2.2109988) q[2];
sx q[2];
rz(-1.8162883) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9671772) q[1];
sx q[1];
rz(-1.6545873) q[1];
sx q[1];
rz(-1.8248952) q[1];
rz(-pi) q[2];
rz(-1.9181973) q[3];
sx q[3];
rz(-0.51651556) q[3];
sx q[3];
rz(1.7361189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3652304) q[2];
sx q[2];
rz(-0.47351101) q[2];
sx q[2];
rz(-3.0253809) q[2];
rz(2.727437) q[3];
sx q[3];
rz(-2.0678346) q[3];
sx q[3];
rz(-2.3381086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6666343) q[0];
sx q[0];
rz(-1.8284429) q[0];
sx q[0];
rz(0.83589244) q[0];
rz(1.6235141) q[1];
sx q[1];
rz(-1.6267136) q[1];
sx q[1];
rz(-1.2380884) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94592124) q[0];
sx q[0];
rz(-2.4950881) q[0];
sx q[0];
rz(-2.2158428) q[0];
rz(-pi) q[1];
rz(2.798271) q[2];
sx q[2];
rz(-1.3385217) q[2];
sx q[2];
rz(-2.6417123) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.064621335) q[1];
sx q[1];
rz(-3.0107452) q[1];
sx q[1];
rz(-2.0276643) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79583056) q[3];
sx q[3];
rz(-2.748877) q[3];
sx q[3];
rz(-2.8098729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28430024) q[2];
sx q[2];
rz(-0.14480536) q[2];
sx q[2];
rz(2.4227552) q[2];
rz(-2.5607064) q[3];
sx q[3];
rz(-1.418117) q[3];
sx q[3];
rz(0.7287997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.523664) q[0];
sx q[0];
rz(-1.5939465) q[0];
sx q[0];
rz(-2.0148328) q[0];
rz(1.2976546) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(-2.0984971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8617423) q[0];
sx q[0];
rz(-1.6011097) q[0];
sx q[0];
rz(1.3002731) q[0];
x q[1];
rz(-0.2506855) q[2];
sx q[2];
rz(-1.0480685) q[2];
sx q[2];
rz(-1.2831068) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7876387) q[1];
sx q[1];
rz(-1.8754205) q[1];
sx q[1];
rz(2.980176) q[1];
rz(-pi) q[2];
rz(0.11123379) q[3];
sx q[3];
rz(-2.4832442) q[3];
sx q[3];
rz(-1.2963795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2780693) q[2];
sx q[2];
rz(-1.7570644) q[2];
sx q[2];
rz(-1.0370022) q[2];
rz(0.95748025) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(0.83468848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0015513) q[0];
sx q[0];
rz(-2.934444) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(2.7032848) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(-1.3137438) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2976332) q[0];
sx q[0];
rz(-2.5974413) q[0];
sx q[0];
rz(2.0354969) q[0];
x q[1];
rz(-1.6694267) q[2];
sx q[2];
rz(-1.5336138) q[2];
sx q[2];
rz(1.94869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.61699762) q[1];
sx q[1];
rz(-1.4819078) q[1];
sx q[1];
rz(1.7434381) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3145218) q[3];
sx q[3];
rz(-1.7328784) q[3];
sx q[3];
rz(0.42419285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4801415) q[2];
sx q[2];
rz(-1.9037312) q[2];
sx q[2];
rz(-0.56337774) q[2];
rz(-1.2207458) q[3];
sx q[3];
rz(-1.3910339) q[3];
sx q[3];
rz(-0.63108546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1221984) q[0];
sx q[0];
rz(-2.4833184) q[0];
sx q[0];
rz(-3.1296545) q[0];
rz(2.7519233) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(-0.24519244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48680533) q[0];
sx q[0];
rz(-2.0005032) q[0];
sx q[0];
rz(-0.91581099) q[0];
rz(-pi) q[1];
rz(0.25845627) q[2];
sx q[2];
rz(-1.788013) q[2];
sx q[2];
rz(2.1803792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4575544) q[1];
sx q[1];
rz(-1.1361309) q[1];
sx q[1];
rz(-0.60575374) q[1];
rz(-pi) q[2];
rz(0.61556365) q[3];
sx q[3];
rz(-2.6608753) q[3];
sx q[3];
rz(1.2445039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37990722) q[2];
sx q[2];
rz(-1.8319538) q[2];
sx q[2];
rz(1.7725819) q[2];
rz(-2.6109429) q[3];
sx q[3];
rz(-2.4574418) q[3];
sx q[3];
rz(-2.0556889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.948792) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(-0.2970933) q[0];
rz(-1.6311215) q[1];
sx q[1];
rz(-2.511697) q[1];
sx q[1];
rz(0.43103257) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1566136) q[0];
sx q[0];
rz(-0.3604069) q[0];
sx q[0];
rz(0.63215881) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9210429) q[2];
sx q[2];
rz(-0.23633453) q[2];
sx q[2];
rz(-2.5823809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4342793) q[1];
sx q[1];
rz(-2.523604) q[1];
sx q[1];
rz(0.98843482) q[1];
x q[2];
rz(-2.0991574) q[3];
sx q[3];
rz(-1.3155121) q[3];
sx q[3];
rz(2.0057037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83583528) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(-0.26710278) q[2];
rz(0.032020656) q[3];
sx q[3];
rz(-1.9710385) q[3];
sx q[3];
rz(-3.035868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.06374643) q[0];
sx q[0];
rz(-1.8505322) q[0];
sx q[0];
rz(2.5015976) q[0];
rz(0.85482875) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(1.7291501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4025637) q[0];
sx q[0];
rz(-1.978658) q[0];
sx q[0];
rz(0.5824851) q[0];
rz(1.0785901) q[2];
sx q[2];
rz(-1.6160496) q[2];
sx q[2];
rz(2.947888) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.392726) q[1];
sx q[1];
rz(-2.7905474) q[1];
sx q[1];
rz(0.24715329) q[1];
rz(1.1453923) q[3];
sx q[3];
rz(-1.914905) q[3];
sx q[3];
rz(0.95208012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5414446) q[2];
sx q[2];
rz(-1.7057799) q[2];
sx q[2];
rz(2.7195462) q[2];
rz(-0.47354928) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(-2.9464909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37050978) q[0];
sx q[0];
rz(-0.93733731) q[0];
sx q[0];
rz(-1.7899845) q[0];
rz(2.8260258) q[1];
sx q[1];
rz(-1.6866997) q[1];
sx q[1];
rz(-1.1531166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6692659) q[0];
sx q[0];
rz(-1.5184222) q[0];
sx q[0];
rz(0.050950296) q[0];
rz(-pi) q[1];
rz(-0.53535992) q[2];
sx q[2];
rz(-1.6125814) q[2];
sx q[2];
rz(0.9524065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0259195) q[1];
sx q[1];
rz(-2.5458836) q[1];
sx q[1];
rz(2.9133948) q[1];
x q[2];
rz(-0.38001506) q[3];
sx q[3];
rz(-1.9560062) q[3];
sx q[3];
rz(-2.4917701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45712581) q[2];
sx q[2];
rz(-1.212333) q[2];
sx q[2];
rz(-0.37377629) q[2];
rz(1.2260381) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(-1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1821197) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(1.3207588) q[0];
rz(-2.2044115) q[1];
sx q[1];
rz(-1.8056185) q[1];
sx q[1];
rz(-2.6722867) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.824911) q[0];
sx q[0];
rz(-1.6880205) q[0];
sx q[0];
rz(0.076923142) q[0];
rz(-1.2861757) q[2];
sx q[2];
rz(-0.60065833) q[2];
sx q[2];
rz(-2.9702368) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9633858) q[1];
sx q[1];
rz(-1.4904516) q[1];
sx q[1];
rz(-0.15060138) q[1];
rz(-1.3239884) q[3];
sx q[3];
rz(-1.863552) q[3];
sx q[3];
rz(-2.5250662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70539537) q[2];
sx q[2];
rz(-2.2787978) q[2];
sx q[2];
rz(1.2104642) q[2];
rz(2.5734899) q[3];
sx q[3];
rz(-1.7567239) q[3];
sx q[3];
rz(-0.97584045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9217011) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(2.2484491) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(-2.91165) q[2];
sx q[2];
rz(-1.8529006) q[2];
sx q[2];
rz(0.012250031) q[2];
rz(1.0004956) q[3];
sx q[3];
rz(-1.9420997) q[3];
sx q[3];
rz(1.6122769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
