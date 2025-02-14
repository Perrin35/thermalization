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
rz(-2.6920707) q[0];
sx q[0];
rz(-1.7188526) q[0];
sx q[0];
rz(2.0916405) q[0];
rz(2.6226251) q[1];
sx q[1];
rz(1.5331886) q[1];
sx q[1];
rz(9.4756995) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22171709) q[0];
sx q[0];
rz(-1.1394986) q[0];
sx q[0];
rz(-1.2141411) q[0];
rz(-2.9257183) q[2];
sx q[2];
rz(-1.2435438) q[2];
sx q[2];
rz(2.808771) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.66827735) q[1];
sx q[1];
rz(-0.44679579) q[1];
sx q[1];
rz(0.75832851) q[1];
rz(-2.256278) q[3];
sx q[3];
rz(-1.3847794) q[3];
sx q[3];
rz(-1.3929588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.27429399) q[2];
sx q[2];
rz(-0.87624246) q[2];
sx q[2];
rz(1.743861) q[2];
rz(-0.45462576) q[3];
sx q[3];
rz(-0.8780829) q[3];
sx q[3];
rz(1.4286058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6468663) q[0];
sx q[0];
rz(-0.88748256) q[0];
sx q[0];
rz(-2.5383762) q[0];
rz(-1.8869205) q[1];
sx q[1];
rz(-0.61379495) q[1];
sx q[1];
rz(2.5616554) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2771526) q[0];
sx q[0];
rz(-1.7371558) q[0];
sx q[0];
rz(-0.54085462) q[0];
rz(-pi) q[1];
rz(-2.4926212) q[2];
sx q[2];
rz(-0.67788092) q[2];
sx q[2];
rz(-3.0650795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.825527) q[1];
sx q[1];
rz(-0.60548079) q[1];
sx q[1];
rz(1.6405296) q[1];
x q[2];
rz(-2.5865057) q[3];
sx q[3];
rz(-0.64315851) q[3];
sx q[3];
rz(-0.023338524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6569528) q[2];
sx q[2];
rz(-2.3886949) q[2];
sx q[2];
rz(2.7638655) q[2];
rz(1.7083302) q[3];
sx q[3];
rz(-1.6323615) q[3];
sx q[3];
rz(1.9506955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3483873) q[0];
sx q[0];
rz(-0.054940104) q[0];
sx q[0];
rz(1.4980263) q[0];
rz(1.161423) q[1];
sx q[1];
rz(-1.252251) q[1];
sx q[1];
rz(-0.84322554) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7347618) q[0];
sx q[0];
rz(-0.94862445) q[0];
sx q[0];
rz(-0.0733331) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9887504) q[2];
sx q[2];
rz(-1.770263) q[2];
sx q[2];
rz(-2.0355088) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49065152) q[1];
sx q[1];
rz(-2.2434714) q[1];
sx q[1];
rz(2.9789799) q[1];
x q[2];
rz(-2.7494861) q[3];
sx q[3];
rz(-2.4014353) q[3];
sx q[3];
rz(1.3010493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8680385) q[2];
sx q[2];
rz(-0.66323438) q[2];
sx q[2];
rz(2.3466477) q[2];
rz(2.3718209) q[3];
sx q[3];
rz(-0.92307463) q[3];
sx q[3];
rz(-2.9845089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.58761) q[0];
sx q[0];
rz(-1.8164604) q[0];
sx q[0];
rz(-0.28833589) q[0];
rz(0.68710697) q[1];
sx q[1];
rz(-1.4930875) q[1];
sx q[1];
rz(-1.7126602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9754959) q[0];
sx q[0];
rz(-1.5773443) q[0];
sx q[0];
rz(-3.1295144) q[0];
x q[1];
rz(-1.4304787) q[2];
sx q[2];
rz(-1.9694917) q[2];
sx q[2];
rz(1.2600419) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1667119) q[1];
sx q[1];
rz(-2.877241) q[1];
sx q[1];
rz(-0.2522142) q[1];
x q[2];
rz(0.89844729) q[3];
sx q[3];
rz(-1.0440376) q[3];
sx q[3];
rz(-1.4970327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.571542) q[2];
sx q[2];
rz(-2.1752581) q[2];
sx q[2];
rz(-2.9759882) q[2];
rz(0.064519493) q[3];
sx q[3];
rz(-3.0118628) q[3];
sx q[3];
rz(1.8701514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22873838) q[0];
sx q[0];
rz(-1.1930635) q[0];
sx q[0];
rz(-2.2304529) q[0];
rz(-0.97081026) q[1];
sx q[1];
rz(-1.3263005) q[1];
sx q[1];
rz(-0.86311805) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2902888) q[0];
sx q[0];
rz(-0.57803854) q[0];
sx q[0];
rz(2.1208288) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4865033) q[2];
sx q[2];
rz(-2.8721923) q[2];
sx q[2];
rz(-1.0625372) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1974302) q[1];
sx q[1];
rz(-1.3384377) q[1];
sx q[1];
rz(-2.5974823) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.04539312) q[3];
sx q[3];
rz(-1.2587341) q[3];
sx q[3];
rz(1.6873716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0866278) q[2];
sx q[2];
rz(-1.2530155) q[2];
sx q[2];
rz(-0.48913726) q[2];
rz(2.1431811) q[3];
sx q[3];
rz(-2.9676134) q[3];
sx q[3];
rz(-3.0424931) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973706) q[0];
sx q[0];
rz(-0.64592823) q[0];
sx q[0];
rz(-0.55322629) q[0];
rz(-0.79752254) q[1];
sx q[1];
rz(-2.1452466) q[1];
sx q[1];
rz(-1.9901989) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40276179) q[0];
sx q[0];
rz(-0.58332764) q[0];
sx q[0];
rz(2.8999694) q[0];
rz(-0.2739353) q[2];
sx q[2];
rz(-2.5735852) q[2];
sx q[2];
rz(0.97962475) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0480736) q[1];
sx q[1];
rz(-2.8052605) q[1];
sx q[1];
rz(-0.63365714) q[1];
x q[2];
rz(-0.38383936) q[3];
sx q[3];
rz(-2.1999) q[3];
sx q[3];
rz(1.5590547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3449751) q[2];
sx q[2];
rz(-1.4130219) q[2];
sx q[2];
rz(0.83149347) q[2];
rz(1.8031395) q[3];
sx q[3];
rz(-1.0748539) q[3];
sx q[3];
rz(2.2510546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.3727342) q[0];
sx q[0];
rz(-1.9356198) q[0];
sx q[0];
rz(-0.15705577) q[0];
rz(1.318469) q[1];
sx q[1];
rz(-2.1993957) q[1];
sx q[1];
rz(-2.7838321) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6899932) q[0];
sx q[0];
rz(-0.66339657) q[0];
sx q[0];
rz(0.76865743) q[0];
rz(-2.2382508) q[2];
sx q[2];
rz(-1.2111665) q[2];
sx q[2];
rz(-2.0697167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5562812) q[1];
sx q[1];
rz(-1.0204633) q[1];
sx q[1];
rz(-2.354391) q[1];
rz(-pi) q[2];
rz(-1.412343) q[3];
sx q[3];
rz(-1.0020578) q[3];
sx q[3];
rz(-0.67711745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1498108) q[2];
sx q[2];
rz(-1.8842183) q[2];
sx q[2];
rz(0.020616654) q[2];
rz(-1.2425544) q[3];
sx q[3];
rz(-0.81762448) q[3];
sx q[3];
rz(0.29153618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6108516) q[0];
sx q[0];
rz(-1.1854956) q[0];
sx q[0];
rz(0.32671842) q[0];
rz(0.84699455) q[1];
sx q[1];
rz(-2.0920483) q[1];
sx q[1];
rz(-1.0583896) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9237019) q[0];
sx q[0];
rz(-1.0578706) q[0];
sx q[0];
rz(-3.0615527) q[0];
rz(-pi) q[1];
rz(1.3943259) q[2];
sx q[2];
rz(-2.2024849) q[2];
sx q[2];
rz(1.2588009) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8673593) q[1];
sx q[1];
rz(-1.6489442) q[1];
sx q[1];
rz(-1.9539425) q[1];
rz(-pi) q[2];
rz(2.045481) q[3];
sx q[3];
rz(-0.62020909) q[3];
sx q[3];
rz(-1.1720464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1191001) q[2];
sx q[2];
rz(-2.7969226) q[2];
sx q[2];
rz(1.1723088) q[2];
rz(-3.1398224) q[3];
sx q[3];
rz(-2.6383196) q[3];
sx q[3];
rz(-0.9790023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718647) q[0];
sx q[0];
rz(-0.80369049) q[0];
sx q[0];
rz(1.3386238) q[0];
rz(1.9610693) q[1];
sx q[1];
rz(-1.0931284) q[1];
sx q[1];
rz(0.48042935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051078253) q[0];
sx q[0];
rz(-1.773293) q[0];
sx q[0];
rz(-1.1064873) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6132203) q[2];
sx q[2];
rz(-1.7406751) q[2];
sx q[2];
rz(0.9695878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5753382) q[1];
sx q[1];
rz(-1.6702139) q[1];
sx q[1];
rz(-1.0097617) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3442743) q[3];
sx q[3];
rz(-1.9206646) q[3];
sx q[3];
rz(-2.1893152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7052445) q[2];
sx q[2];
rz(-0.99506012) q[2];
sx q[2];
rz(-2.056541) q[2];
rz(2.0294225) q[3];
sx q[3];
rz(-0.89434904) q[3];
sx q[3];
rz(2.3280242) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28806624) q[0];
sx q[0];
rz(-0.94678322) q[0];
sx q[0];
rz(-2.1296401) q[0];
rz(-0.90156737) q[1];
sx q[1];
rz(-2.8755867) q[1];
sx q[1];
rz(2.576135) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31663592) q[0];
sx q[0];
rz(-0.68616435) q[0];
sx q[0];
rz(-1.2385204) q[0];
x q[1];
rz(-1.0305963) q[2];
sx q[2];
rz(-2.6517617) q[2];
sx q[2];
rz(2.4649807) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31670891) q[1];
sx q[1];
rz(-1.4780586) q[1];
sx q[1];
rz(2.2888661) q[1];
x q[2];
rz(-0.10528414) q[3];
sx q[3];
rz(-1.255688) q[3];
sx q[3];
rz(1.1694825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3507583) q[2];
sx q[2];
rz(-1.5956722) q[2];
sx q[2];
rz(-1.2606384) q[2];
rz(-0.44220051) q[3];
sx q[3];
rz(-2.111777) q[3];
sx q[3];
rz(1.7650167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5762536) q[0];
sx q[0];
rz(-2.2591142) q[0];
sx q[0];
rz(2.2478065) q[0];
rz(-1.5743938) q[1];
sx q[1];
rz(-1.6849453) q[1];
sx q[1];
rz(-1.4243855) q[1];
rz(1.2906176) q[2];
sx q[2];
rz(-1.9815255) q[2];
sx q[2];
rz(-1.6407914) q[2];
rz(1.1461729) q[3];
sx q[3];
rz(-2.7594447) q[3];
sx q[3];
rz(-1.2654163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
