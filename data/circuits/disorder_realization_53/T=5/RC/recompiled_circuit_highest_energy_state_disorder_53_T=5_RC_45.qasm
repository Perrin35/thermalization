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
rz(-0.1102912) q[0];
sx q[0];
rz(4.9977367) q[0];
sx q[0];
rz(9.1060299) q[0];
rz(1.1822074) q[1];
sx q[1];
rz(-1.6637586) q[1];
sx q[1];
rz(-1.1629265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1308474) q[0];
sx q[0];
rz(-1.3358572) q[0];
sx q[0];
rz(1.9128591) q[0];
rz(3.0672795) q[2];
sx q[2];
rz(-0.19467672) q[2];
sx q[2];
rz(-3.1271324) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8713177) q[1];
sx q[1];
rz(-1.2658889) q[1];
sx q[1];
rz(-0.31590806) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5507023) q[3];
sx q[3];
rz(-2.1455975) q[3];
sx q[3];
rz(-1.2736831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9691539) q[2];
sx q[2];
rz(-2.6814851) q[2];
sx q[2];
rz(-0.13574204) q[2];
rz(2.8067449) q[3];
sx q[3];
rz(-1.2232774) q[3];
sx q[3];
rz(0.39723435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44884509) q[0];
sx q[0];
rz(-2.1362342) q[0];
sx q[0];
rz(-2.1951065) q[0];
rz(-1.954156) q[1];
sx q[1];
rz(-0.58381909) q[1];
sx q[1];
rz(0.28396398) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0683354) q[0];
sx q[0];
rz(-2.8767902) q[0];
sx q[0];
rz(-0.7629443) q[0];
rz(-pi) q[1];
rz(-1.8503485) q[2];
sx q[2];
rz(-2.4233305) q[2];
sx q[2];
rz(-0.58092434) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75278209) q[1];
sx q[1];
rz(-1.4335263) q[1];
sx q[1];
rz(-2.7146882) q[1];
rz(1.918774) q[3];
sx q[3];
rz(-1.9377515) q[3];
sx q[3];
rz(1.8311794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2076608) q[2];
sx q[2];
rz(-1.8457103) q[2];
sx q[2];
rz(-2.3409823) q[2];
rz(-1.8396395) q[3];
sx q[3];
rz(-1.1336361) q[3];
sx q[3];
rz(3.083526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0170853) q[0];
sx q[0];
rz(-2.7094816) q[0];
sx q[0];
rz(-2.3486163) q[0];
rz(-2.4555581) q[1];
sx q[1];
rz(-1.425309) q[1];
sx q[1];
rz(-2.3763903) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3065465) q[0];
sx q[0];
rz(-1.5126499) q[0];
sx q[0];
rz(1.3594199) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5386109) q[2];
sx q[2];
rz(-1.5192521) q[2];
sx q[2];
rz(-2.3508765) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6531892) q[1];
sx q[1];
rz(-2.326366) q[1];
sx q[1];
rz(-1.1865739) q[1];
x q[2];
rz(2.3451557) q[3];
sx q[3];
rz(-2.0579866) q[3];
sx q[3];
rz(1.1034484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5746295) q[2];
sx q[2];
rz(-2.4506863) q[2];
sx q[2];
rz(1.8951269) q[2];
rz(1.1860819) q[3];
sx q[3];
rz(-1.3776774) q[3];
sx q[3];
rz(-2.932909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801159) q[0];
sx q[0];
rz(-1.8210541) q[0];
sx q[0];
rz(-3.0082974) q[0];
rz(-2.8721299) q[1];
sx q[1];
rz(-2.2270484) q[1];
sx q[1];
rz(-2.6752313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15233697) q[0];
sx q[0];
rz(-1.6902802) q[0];
sx q[0];
rz(2.5463922) q[0];
x q[1];
rz(0.94302098) q[2];
sx q[2];
rz(-1.7258712) q[2];
sx q[2];
rz(0.53395203) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.129648) q[1];
sx q[1];
rz(-2.5135871) q[1];
sx q[1];
rz(2.8824214) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5047938) q[3];
sx q[3];
rz(-1.0579234) q[3];
sx q[3];
rz(-2.6309225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2103297) q[2];
sx q[2];
rz(-0.42079058) q[2];
sx q[2];
rz(2.8974864) q[2];
rz(-1.7804451) q[3];
sx q[3];
rz(-2.1972392) q[3];
sx q[3];
rz(1.9349499) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72531438) q[0];
sx q[0];
rz(-0.84753528) q[0];
sx q[0];
rz(-0.086294802) q[0];
rz(1.3823973) q[1];
sx q[1];
rz(-2.081213) q[1];
sx q[1];
rz(1.28654) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2479682) q[0];
sx q[0];
rz(-2.7148406) q[0];
sx q[0];
rz(-1.990699) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3628144) q[2];
sx q[2];
rz(-1.7302824) q[2];
sx q[2];
rz(0.82568141) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5865422) q[1];
sx q[1];
rz(-0.70603958) q[1];
sx q[1];
rz(-1.6976188) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52545737) q[3];
sx q[3];
rz(-2.1920929) q[3];
sx q[3];
rz(0.8233101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7354108) q[2];
sx q[2];
rz(-1.4124796) q[2];
sx q[2];
rz(2.0162876) q[2];
rz(0.79648894) q[3];
sx q[3];
rz(-0.73553604) q[3];
sx q[3];
rz(-2.9430732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7859802) q[0];
sx q[0];
rz(-0.2934083) q[0];
sx q[0];
rz(0.32817131) q[0];
rz(-0.38901058) q[1];
sx q[1];
rz(-1.5704472) q[1];
sx q[1];
rz(1.6860115) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9756115) q[0];
sx q[0];
rz(-1.360438) q[0];
sx q[0];
rz(1.4882404) q[0];
rz(-pi) q[1];
rz(-3.0148618) q[2];
sx q[2];
rz(-1.4369082) q[2];
sx q[2];
rz(0.19370475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4715188) q[1];
sx q[1];
rz(-0.68738261) q[1];
sx q[1];
rz(1.2015461) q[1];
x q[2];
rz(-1.2120655) q[3];
sx q[3];
rz(-1.609612) q[3];
sx q[3];
rz(1.9129378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97633156) q[2];
sx q[2];
rz(-0.92504048) q[2];
sx q[2];
rz(-0.46249214) q[2];
rz(-1.0235323) q[3];
sx q[3];
rz(-2.0868802) q[3];
sx q[3];
rz(2.3033477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3384712) q[0];
sx q[0];
rz(-1.7110889) q[0];
sx q[0];
rz(2.9679003) q[0];
rz(2.9202785) q[1];
sx q[1];
rz(-0.68562713) q[1];
sx q[1];
rz(1.7783222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46745978) q[0];
sx q[0];
rz(-1.9741575) q[0];
sx q[0];
rz(2.9441071) q[0];
rz(-pi) q[1];
rz(-3.0603833) q[2];
sx q[2];
rz(-1.6575362) q[2];
sx q[2];
rz(-2.7600811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9611836) q[1];
sx q[1];
rz(-0.81432691) q[1];
sx q[1];
rz(2.4717038) q[1];
rz(-pi) q[2];
rz(1.8425499) q[3];
sx q[3];
rz(-1.8044961) q[3];
sx q[3];
rz(2.8444447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8814016) q[2];
sx q[2];
rz(-1.9804201) q[2];
sx q[2];
rz(-1.4455522) q[2];
rz(-0.77406231) q[3];
sx q[3];
rz(-0.71662199) q[3];
sx q[3];
rz(1.6707481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5792907) q[0];
sx q[0];
rz(-2.2023872) q[0];
sx q[0];
rz(-0.24765177) q[0];
rz(0.66894764) q[1];
sx q[1];
rz(-1.1890143) q[1];
sx q[1];
rz(-0.98562366) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068531903) q[0];
sx q[0];
rz(-1.4559696) q[0];
sx q[0];
rz(2.5825744) q[0];
rz(-1.2560166) q[2];
sx q[2];
rz(-2.1119962) q[2];
sx q[2];
rz(-1.4044631) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0874153) q[1];
sx q[1];
rz(-1.8761411) q[1];
sx q[1];
rz(-1.1609689) q[1];
x q[2];
rz(2.5466315) q[3];
sx q[3];
rz(-1.5719126) q[3];
sx q[3];
rz(3.0240455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0011657) q[2];
sx q[2];
rz(-0.75158921) q[2];
sx q[2];
rz(-1.8482194) q[2];
rz(-1.9017396) q[3];
sx q[3];
rz(-0.69609061) q[3];
sx q[3];
rz(-2.4333439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.86466113) q[0];
sx q[0];
rz(-1.4286574) q[0];
sx q[0];
rz(2.1760333) q[0];
rz(-2.7652265) q[1];
sx q[1];
rz(-1.3220359) q[1];
sx q[1];
rz(1.9115062) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7853785) q[0];
sx q[0];
rz(-1.5206771) q[0];
sx q[0];
rz(-0.26595195) q[0];
rz(1.240991) q[2];
sx q[2];
rz(-1.1993186) q[2];
sx q[2];
rz(2.7895751) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9113637) q[1];
sx q[1];
rz(-1.5004284) q[1];
sx q[1];
rz(2.5350948) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3283703) q[3];
sx q[3];
rz(-1.5322161) q[3];
sx q[3];
rz(0.67594066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2824715) q[2];
sx q[2];
rz(-2.4716447) q[2];
sx q[2];
rz(0.13151375) q[2];
rz(-0.48111835) q[3];
sx q[3];
rz(-0.93004623) q[3];
sx q[3];
rz(0.49829811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1199101) q[0];
sx q[0];
rz(-2.5639738) q[0];
sx q[0];
rz(-0.48900327) q[0];
rz(1.3267714) q[1];
sx q[1];
rz(-1.2811456) q[1];
sx q[1];
rz(0.16924032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45825935) q[0];
sx q[0];
rz(-1.1348327) q[0];
sx q[0];
rz(0.13157121) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18074482) q[2];
sx q[2];
rz(-1.6588447) q[2];
sx q[2];
rz(-3.0918025) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9368912) q[1];
sx q[1];
rz(-1.1549779) q[1];
sx q[1];
rz(-0.73483006) q[1];
rz(-pi) q[2];
rz(-0.80767969) q[3];
sx q[3];
rz(-0.39356222) q[3];
sx q[3];
rz(-0.53787947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5532316) q[2];
sx q[2];
rz(-1.0855805) q[2];
sx q[2];
rz(-0.21151839) q[2];
rz(1.616098) q[3];
sx q[3];
rz(-2.1414521) q[3];
sx q[3];
rz(0.26743993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.78032988) q[0];
sx q[0];
rz(-1.7050671) q[0];
sx q[0];
rz(0.29062511) q[0];
rz(0.81193874) q[1];
sx q[1];
rz(-0.62807905) q[1];
sx q[1];
rz(-1.5608578) q[1];
rz(2.7079034) q[2];
sx q[2];
rz(-1.5530752) q[2];
sx q[2];
rz(-1.8032522) q[2];
rz(2.3216861) q[3];
sx q[3];
rz(-1.0884566) q[3];
sx q[3];
rz(1.4940445) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
