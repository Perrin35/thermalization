OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1898243) q[0];
sx q[0];
rz(-0.70280743) q[0];
sx q[0];
rz(-1.9195317) q[0];
rz(3.089978) q[1];
sx q[1];
rz(-1.6366704) q[1];
sx q[1];
rz(1.6610891) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015704) q[0];
sx q[0];
rz(-0.68212494) q[0];
sx q[0];
rz(-1.3093824) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7403549) q[2];
sx q[2];
rz(-1.6785067) q[2];
sx q[2];
rz(-2.8603399) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29917078) q[1];
sx q[1];
rz(-1.7452733) q[1];
sx q[1];
rz(-2.7122696) q[1];
rz(1.0869671) q[3];
sx q[3];
rz(-1.1796406) q[3];
sx q[3];
rz(-3.1227171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1808971) q[2];
sx q[2];
rz(-2.3263859) q[2];
sx q[2];
rz(-1.0308456) q[2];
rz(2.893462) q[3];
sx q[3];
rz(-1.7957567) q[3];
sx q[3];
rz(0.26721755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88028518) q[0];
sx q[0];
rz(-1.1957276) q[0];
sx q[0];
rz(-2.0174568) q[0];
rz(1.0719489) q[1];
sx q[1];
rz(-1.5644466) q[1];
sx q[1];
rz(2.7148278) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.916722) q[0];
sx q[0];
rz(-1.568075) q[0];
sx q[0];
rz(-0.82758521) q[0];
x q[1];
rz(-0.34554225) q[2];
sx q[2];
rz(-1.3239408) q[2];
sx q[2];
rz(-3.0902362) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71449762) q[1];
sx q[1];
rz(-0.87676226) q[1];
sx q[1];
rz(1.7114033) q[1];
rz(-2.9396966) q[3];
sx q[3];
rz(-1.9526953) q[3];
sx q[3];
rz(2.524162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4830604) q[2];
sx q[2];
rz(-2.6250562) q[2];
sx q[2];
rz(3.0413682) q[2];
rz(1.077486) q[3];
sx q[3];
rz(-1.5651549) q[3];
sx q[3];
rz(1.235435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.539262) q[0];
sx q[0];
rz(-0.97719181) q[0];
sx q[0];
rz(1.4343028) q[0];
rz(0.82646787) q[1];
sx q[1];
rz(-1.4974599) q[1];
sx q[1];
rz(0.43063393) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62118351) q[0];
sx q[0];
rz(-0.31776515) q[0];
sx q[0];
rz(-2.5924224) q[0];
rz(-pi) q[1];
rz(-1.8437064) q[2];
sx q[2];
rz(-1.6962972) q[2];
sx q[2];
rz(-2.2019406) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0525749) q[1];
sx q[1];
rz(-0.7886501) q[1];
sx q[1];
rz(0.99747212) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39726302) q[3];
sx q[3];
rz(-2.0797585) q[3];
sx q[3];
rz(-2.404117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87381252) q[2];
sx q[2];
rz(-1.9766821) q[2];
sx q[2];
rz(2.5679892) q[2];
rz(-2.7576647) q[3];
sx q[3];
rz(-1.826518) q[3];
sx q[3];
rz(-2.4887776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7225994) q[0];
sx q[0];
rz(-0.79177952) q[0];
sx q[0];
rz(1.0571887) q[0];
rz(-2.3253697) q[1];
sx q[1];
rz(-2.1482601) q[1];
sx q[1];
rz(-0.23424558) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5200096) q[0];
sx q[0];
rz(-2.710973) q[0];
sx q[0];
rz(-1.3060116) q[0];
rz(-pi) q[1];
rz(-0.90936898) q[2];
sx q[2];
rz(-2.5071835) q[2];
sx q[2];
rz(2.6480791) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.034742753) q[1];
sx q[1];
rz(-2.2487246) q[1];
sx q[1];
rz(0.89660109) q[1];
rz(-pi) q[2];
rz(0.81984249) q[3];
sx q[3];
rz(-0.64717918) q[3];
sx q[3];
rz(-2.1539713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8546042) q[2];
sx q[2];
rz(-1.333326) q[2];
sx q[2];
rz(0.90281478) q[2];
rz(-1.758894) q[3];
sx q[3];
rz(-2.8185676) q[3];
sx q[3];
rz(2.2949016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.46120241) q[0];
sx q[0];
rz(-1.8620551) q[0];
sx q[0];
rz(-2.4317106) q[0];
rz(0.4711802) q[1];
sx q[1];
rz(-1.9147562) q[1];
sx q[1];
rz(3.0772298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.956139) q[0];
sx q[0];
rz(-0.50489932) q[0];
sx q[0];
rz(1.8527777) q[0];
rz(-pi) q[1];
rz(-1.9856312) q[2];
sx q[2];
rz(-2.6081787) q[2];
sx q[2];
rz(-2.2629625) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91436323) q[1];
sx q[1];
rz(-1.0681947) q[1];
sx q[1];
rz(1.0881626) q[1];
rz(-pi) q[2];
rz(-0.086802089) q[3];
sx q[3];
rz(-1.5802812) q[3];
sx q[3];
rz(-2.529151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71089661) q[2];
sx q[2];
rz(-0.5527834) q[2];
sx q[2];
rz(0.95332471) q[2];
rz(-2.583875) q[3];
sx q[3];
rz(-1.6744813) q[3];
sx q[3];
rz(-0.52887708) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8271269) q[0];
sx q[0];
rz(-1.0449266) q[0];
sx q[0];
rz(-0.99660981) q[0];
rz(-3.0820471) q[1];
sx q[1];
rz(-2.153502) q[1];
sx q[1];
rz(1.7893808) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9080862) q[0];
sx q[0];
rz(-2.3829989) q[0];
sx q[0];
rz(-2.3367466) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77162051) q[2];
sx q[2];
rz(-1.3305656) q[2];
sx q[2];
rz(-0.26773237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.775573) q[1];
sx q[1];
rz(-1.2723574) q[1];
sx q[1];
rz(1.2828903) q[1];
rz(-pi) q[2];
rz(1.1354225) q[3];
sx q[3];
rz(-0.7061031) q[3];
sx q[3];
rz(-2.5869334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38702866) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(-2.582029) q[2];
rz(-0.83545056) q[3];
sx q[3];
rz(-0.42568922) q[3];
sx q[3];
rz(-1.3028418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68592042) q[0];
sx q[0];
rz(-0.73021013) q[0];
sx q[0];
rz(-2.7749104) q[0];
rz(0.8815676) q[1];
sx q[1];
rz(-0.63789788) q[1];
sx q[1];
rz(1.7104023) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.139643) q[0];
sx q[0];
rz(-2.3649923) q[0];
sx q[0];
rz(0.5712255) q[0];
rz(-pi) q[1];
rz(-1.9164247) q[2];
sx q[2];
rz(-0.91424886) q[2];
sx q[2];
rz(2.569927) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6122653) q[1];
sx q[1];
rz(-3.021486) q[1];
sx q[1];
rz(-0.24397759) q[1];
rz(-pi) q[2];
rz(-1.335698) q[3];
sx q[3];
rz(-0.92609105) q[3];
sx q[3];
rz(0.39261445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7876579) q[2];
sx q[2];
rz(-2.6255609) q[2];
sx q[2];
rz(2.3616974) q[2];
rz(-2.6025313) q[3];
sx q[3];
rz(-1.6293679) q[3];
sx q[3];
rz(-2.5299634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014538177) q[0];
sx q[0];
rz(-2.4712565) q[0];
sx q[0];
rz(-2.3624453) q[0];
rz(-2.6226793) q[1];
sx q[1];
rz(-1.8592368) q[1];
sx q[1];
rz(0.40294495) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8810324) q[0];
sx q[0];
rz(-1.1024108) q[0];
sx q[0];
rz(2.3163176) q[0];
x q[1];
rz(-2.8267118) q[2];
sx q[2];
rz(-2.0337031) q[2];
sx q[2];
rz(-0.34466448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1379304) q[1];
sx q[1];
rz(-0.76293412) q[1];
sx q[1];
rz(-0.59138815) q[1];
rz(-0.99767606) q[3];
sx q[3];
rz(-1.7955639) q[3];
sx q[3];
rz(2.5483957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0275823) q[2];
sx q[2];
rz(-1.1677914) q[2];
sx q[2];
rz(-1.1172509) q[2];
rz(0.69636017) q[3];
sx q[3];
rz(-1.1098692) q[3];
sx q[3];
rz(-1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6663412) q[0];
sx q[0];
rz(-0.23463686) q[0];
sx q[0];
rz(0.41467211) q[0];
rz(-1.0617725) q[1];
sx q[1];
rz(-2.2551408) q[1];
sx q[1];
rz(-2.3019703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79796004) q[0];
sx q[0];
rz(-2.0164444) q[0];
sx q[0];
rz(-0.4526505) q[0];
rz(-0.62445415) q[2];
sx q[2];
rz(-0.82459282) q[2];
sx q[2];
rz(-2.7471682) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.95123467) q[1];
sx q[1];
rz(-2.2273185) q[1];
sx q[1];
rz(1.3912398) q[1];
rz(-pi) q[2];
rz(-2.5483545) q[3];
sx q[3];
rz(-2.7194033) q[3];
sx q[3];
rz(1.7144904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35905579) q[2];
sx q[2];
rz(-1.3691207) q[2];
sx q[2];
rz(0.071652023) q[2];
rz(3.0960633) q[3];
sx q[3];
rz(-1.9436911) q[3];
sx q[3];
rz(2.6825452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5422106) q[0];
sx q[0];
rz(-0.57294232) q[0];
sx q[0];
rz(1.0546767) q[0];
rz(3.1090464) q[1];
sx q[1];
rz(-0.97144214) q[1];
sx q[1];
rz(-0.5221101) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73682937) q[0];
sx q[0];
rz(-2.1407581) q[0];
sx q[0];
rz(1.205377) q[0];
rz(-pi) q[1];
rz(-0.98507757) q[2];
sx q[2];
rz(-0.92364468) q[2];
sx q[2];
rz(-0.070731846) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0070399) q[1];
sx q[1];
rz(-1.1940446) q[1];
sx q[1];
rz(1.1621848) q[1];
rz(-pi) q[2];
rz(1.9315835) q[3];
sx q[3];
rz(-0.57769247) q[3];
sx q[3];
rz(-3.011774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9944428) q[2];
sx q[2];
rz(-0.20211896) q[2];
sx q[2];
rz(-1.7333376) q[2];
rz(1.555892) q[3];
sx q[3];
rz(-2.9097911) q[3];
sx q[3];
rz(-0.91109341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0178575) q[0];
sx q[0];
rz(-1.3029079) q[0];
sx q[0];
rz(-2.2444176) q[0];
rz(-0.97902117) q[1];
sx q[1];
rz(-1.9985825) q[1];
sx q[1];
rz(-0.40252007) q[1];
rz(-0.19569068) q[2];
sx q[2];
rz(-1.8974081) q[2];
sx q[2];
rz(-1.8477389) q[2];
rz(-1.3835945) q[3];
sx q[3];
rz(-0.91619195) q[3];
sx q[3];
rz(-0.34294101) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
