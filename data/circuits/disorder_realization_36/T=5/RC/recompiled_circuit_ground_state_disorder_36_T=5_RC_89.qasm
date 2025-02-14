OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9517684) q[0];
sx q[0];
rz(3.8444001) q[0];
sx q[0];
rz(8.202717) q[0];
rz(-3.1932073) q[1];
sx q[1];
rz(1.5049223) q[1];
sx q[1];
rz(7.7636889) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015704) q[0];
sx q[0];
rz(-0.68212494) q[0];
sx q[0];
rz(1.8322102) q[0];
rz(-2.1406581) q[2];
sx q[2];
rz(-0.20059948) q[2];
sx q[2];
rz(-2.4127485) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.192321) q[1];
sx q[1];
rz(-1.1484129) q[1];
sx q[1];
rz(1.3793089) q[1];
x q[2];
rz(1.0869671) q[3];
sx q[3];
rz(-1.1796406) q[3];
sx q[3];
rz(0.018875518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9606955) q[2];
sx q[2];
rz(-0.81520671) q[2];
sx q[2];
rz(-2.1107471) q[2];
rz(2.893462) q[3];
sx q[3];
rz(-1.3458359) q[3];
sx q[3];
rz(2.8743751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2613075) q[0];
sx q[0];
rz(-1.945865) q[0];
sx q[0];
rz(-1.1241359) q[0];
rz(1.0719489) q[1];
sx q[1];
rz(-1.5644466) q[1];
sx q[1];
rz(2.7148278) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65157344) q[0];
sx q[0];
rz(-2.314004) q[0];
sx q[0];
rz(-0.0036959582) q[0];
x q[1];
rz(-2.7960504) q[2];
sx q[2];
rz(-1.3239408) q[2];
sx q[2];
rz(3.0902362) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71449762) q[1];
sx q[1];
rz(-0.87676226) q[1];
sx q[1];
rz(-1.7114033) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9396966) q[3];
sx q[3];
rz(-1.1888973) q[3];
sx q[3];
rz(-0.61743067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65853226) q[2];
sx q[2];
rz(-2.6250562) q[2];
sx q[2];
rz(3.0413682) q[2];
rz(-2.0641067) q[3];
sx q[3];
rz(-1.5651549) q[3];
sx q[3];
rz(-1.9061576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.539262) q[0];
sx q[0];
rz(-2.1644008) q[0];
sx q[0];
rz(1.4343028) q[0];
rz(-2.3151248) q[1];
sx q[1];
rz(-1.4974599) q[1];
sx q[1];
rz(-2.7109587) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185558) q[0];
sx q[0];
rz(-1.7346177) q[0];
sx q[0];
rz(-0.27351606) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13027066) q[2];
sx q[2];
rz(-1.8415057) q[2];
sx q[2];
rz(-0.66616466) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7942126) q[1];
sx q[1];
rz(-0.93232226) q[1];
sx q[1];
rz(2.641851) q[1];
rz(-2.1149733) q[3];
sx q[3];
rz(-1.2261571) q[3];
sx q[3];
rz(2.1066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87381252) q[2];
sx q[2];
rz(-1.1649106) q[2];
sx q[2];
rz(2.5679892) q[2];
rz(0.38392797) q[3];
sx q[3];
rz(-1.826518) q[3];
sx q[3];
rz(0.65281502) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4189932) q[0];
sx q[0];
rz(-2.3498131) q[0];
sx q[0];
rz(2.0844039) q[0];
rz(-0.81622299) q[1];
sx q[1];
rz(-2.1482601) q[1];
sx q[1];
rz(-2.9073471) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.809991) q[0];
sx q[0];
rz(-1.9854641) q[0];
sx q[0];
rz(3.0219487) q[0];
rz(-pi) q[1];
rz(-2.2322237) q[2];
sx q[2];
rz(-2.5071835) q[2];
sx q[2];
rz(0.4935136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1409451) q[1];
sx q[1];
rz(-2.0785626) q[1];
sx q[1];
rz(-2.341048) q[1];
rz(-pi) q[2];
rz(-2.6654763) q[3];
sx q[3];
rz(-2.0272539) q[3];
sx q[3];
rz(1.2904087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.28698841) q[2];
sx q[2];
rz(-1.333326) q[2];
sx q[2];
rz(0.90281478) q[2];
rz(-1.3826987) q[3];
sx q[3];
rz(-2.8185676) q[3];
sx q[3];
rz(0.84669101) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46120241) q[0];
sx q[0];
rz(-1.8620551) q[0];
sx q[0];
rz(2.4317106) q[0];
rz(-2.6704125) q[1];
sx q[1];
rz(-1.9147562) q[1];
sx q[1];
rz(-0.064362854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.956139) q[0];
sx q[0];
rz(-0.50489932) q[0];
sx q[0];
rz(-1.8527777) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9079403) q[2];
sx q[2];
rz(-2.0548247) q[2];
sx q[2];
rz(2.7357227) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2272294) q[1];
sx q[1];
rz(-2.0733979) q[1];
sx q[1];
rz(2.05343) q[1];
x q[2];
rz(1.5803171) q[3];
sx q[3];
rz(-1.4839982) q[3];
sx q[3];
rz(0.95918005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71089661) q[2];
sx q[2];
rz(-0.5527834) q[2];
sx q[2];
rz(-0.95332471) q[2];
rz(2.583875) q[3];
sx q[3];
rz(-1.4671114) q[3];
sx q[3];
rz(2.6127156) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3144658) q[0];
sx q[0];
rz(-2.0966661) q[0];
sx q[0];
rz(0.99660981) q[0];
rz(-0.059545513) q[1];
sx q[1];
rz(-0.98809067) q[1];
sx q[1];
rz(1.7893808) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9080862) q[0];
sx q[0];
rz(-0.75859374) q[0];
sx q[0];
rz(0.80484606) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33783317) q[2];
sx q[2];
rz(-2.3408836) q[2];
sx q[2];
rz(1.5429301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88207605) q[1];
sx q[1];
rz(-1.8456371) q[1];
sx q[1];
rz(2.8311353) q[1];
rz(0.91259802) q[3];
sx q[3];
rz(-1.8479947) q[3];
sx q[3];
rz(0.67597055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38702866) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(-0.55956364) q[2];
rz(-2.3061421) q[3];
sx q[3];
rz(-2.7159034) q[3];
sx q[3];
rz(-1.3028418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.68592042) q[0];
sx q[0];
rz(-0.73021013) q[0];
sx q[0];
rz(-0.36668229) q[0];
rz(-0.8815676) q[1];
sx q[1];
rz(-0.63789788) q[1];
sx q[1];
rz(1.4311904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.139643) q[0];
sx q[0];
rz(-2.3649923) q[0];
sx q[0];
rz(-2.5703672) q[0];
rz(-pi) q[1];
rz(-1.9164247) q[2];
sx q[2];
rz(-0.91424886) q[2];
sx q[2];
rz(-0.57166568) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28364946) q[1];
sx q[1];
rz(-1.4542631) q[1];
sx q[1];
rz(-1.5416508) q[1];
rz(-1.8058947) q[3];
sx q[3];
rz(-2.2155016) q[3];
sx q[3];
rz(0.39261445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7876579) q[2];
sx q[2];
rz(-2.6255609) q[2];
sx q[2];
rz(-0.77989522) q[2];
rz(-2.6025313) q[3];
sx q[3];
rz(-1.5122248) q[3];
sx q[3];
rz(2.5299634) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270545) q[0];
sx q[0];
rz(-0.67033613) q[0];
sx q[0];
rz(2.3624453) q[0];
rz(-0.51891333) q[1];
sx q[1];
rz(-1.8592368) q[1];
sx q[1];
rz(-0.40294495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76501827) q[0];
sx q[0];
rz(-2.2857762) q[0];
sx q[0];
rz(0.9299703) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31488089) q[2];
sx q[2];
rz(-1.1078896) q[2];
sx q[2];
rz(-0.34466448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.886595) q[1];
sx q[1];
rz(-0.95980058) q[1];
sx q[1];
rz(-2.0605036) q[1];
x q[2];
rz(-2.8759148) q[3];
sx q[3];
rz(-1.0138265) q[3];
sx q[3];
rz(1.1204615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11401033) q[2];
sx q[2];
rz(-1.9738013) q[2];
sx q[2];
rz(-1.1172509) q[2];
rz(-2.4452325) q[3];
sx q[3];
rz(-1.1098692) q[3];
sx q[3];
rz(-1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6663412) q[0];
sx q[0];
rz(-0.23463686) q[0];
sx q[0];
rz(-0.41467211) q[0];
rz(-1.0617725) q[1];
sx q[1];
rz(-2.2551408) q[1];
sx q[1];
rz(0.83962238) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0939621) q[0];
sx q[0];
rz(-2.5174401) q[0];
sx q[0];
rz(-2.312129) q[0];
rz(-pi) q[1];
rz(-2.1346852) q[2];
sx q[2];
rz(-0.93265763) q[2];
sx q[2];
rz(1.9319122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66218195) q[1];
sx q[1];
rz(-0.67711293) q[1];
sx q[1];
rz(0.22775316) q[1];
x q[2];
rz(-2.5483545) q[3];
sx q[3];
rz(-2.7194033) q[3];
sx q[3];
rz(-1.4271023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35905579) q[2];
sx q[2];
rz(-1.3691207) q[2];
sx q[2];
rz(0.071652023) q[2];
rz(3.0960633) q[3];
sx q[3];
rz(-1.1979016) q[3];
sx q[3];
rz(0.45904747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5422106) q[0];
sx q[0];
rz(-0.57294232) q[0];
sx q[0];
rz(-1.0546767) q[0];
rz(-3.1090464) q[1];
sx q[1];
rz(-0.97144214) q[1];
sx q[1];
rz(0.5221101) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73682937) q[0];
sx q[0];
rz(-1.0008345) q[0];
sx q[0];
rz(1.9362157) q[0];
x q[1];
rz(-2.5100461) q[2];
sx q[2];
rz(-0.8435404) q[2];
sx q[2];
rz(0.90383672) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5473816) q[1];
sx q[1];
rz(-1.1923596) q[1];
sx q[1];
rz(-0.40706472) q[1];
rz(1.9315835) q[3];
sx q[3];
rz(-2.5639002) q[3];
sx q[3];
rz(-0.1298187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1471499) q[2];
sx q[2];
rz(-2.9394737) q[2];
sx q[2];
rz(1.7333376) q[2];
rz(1.5857006) q[3];
sx q[3];
rz(-0.2318016) q[3];
sx q[3];
rz(-0.91109341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0178575) q[0];
sx q[0];
rz(-1.3029079) q[0];
sx q[0];
rz(-2.2444176) q[0];
rz(-2.1625715) q[1];
sx q[1];
rz(-1.1430102) q[1];
sx q[1];
rz(2.7390726) q[1];
rz(-1.903309) q[2];
sx q[2];
rz(-1.3855743) q[2];
sx q[2];
rz(-0.34045548) q[2];
rz(-2.9036941) q[3];
sx q[3];
rz(-2.4645505) q[3];
sx q[3];
rz(-0.6445618) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
