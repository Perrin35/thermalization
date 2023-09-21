OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(-1.3949431) q[0];
sx q[0];
rz(-3.1403132) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(4.2445634) q[1];
sx q[1];
rz(7.0581262) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476386) q[0];
sx q[0];
rz(-1.5210266) q[0];
sx q[0];
rz(2.9988891) q[0];
rz(2.7156419) q[2];
sx q[2];
rz(-0.89086878) q[2];
sx q[2];
rz(-1.3198864) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.98382681) q[1];
sx q[1];
rz(-2.9938934) q[1];
sx q[1];
rz(1.8450792) q[1];
rz(-pi) q[2];
rz(-1.7232056) q[3];
sx q[3];
rz(-1.2890352) q[3];
sx q[3];
rz(1.2581274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(1.1958896) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(1.0999854) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(-1.7659448) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7398867) q[0];
sx q[0];
rz(-1.6580083) q[0];
sx q[0];
rz(0.4102104) q[0];
rz(1.3324845) q[2];
sx q[2];
rz(-1.6147699) q[2];
sx q[2];
rz(-2.6251052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34054204) q[1];
sx q[1];
rz(-1.2060296) q[1];
sx q[1];
rz(1.7629452) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0248915) q[3];
sx q[3];
rz(-2.5149269) q[3];
sx q[3];
rz(2.6729667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21330825) q[2];
sx q[2];
rz(-1.7958612) q[2];
sx q[2];
rz(2.6518872) q[2];
rz(2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-1.0666696) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5415444) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(-1.2352357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89861682) q[0];
sx q[0];
rz(-1.9447864) q[0];
sx q[0];
rz(0.6181194) q[0];
rz(-1.3182993) q[2];
sx q[2];
rz(-1.9324979) q[2];
sx q[2];
rz(-0.30644882) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8894549) q[1];
sx q[1];
rz(-0.368202) q[1];
sx q[1];
rz(-0.77106573) q[1];
rz(-2.5127605) q[3];
sx q[3];
rz(-1.9753846) q[3];
sx q[3];
rz(1.7189327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(0.60825545) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075994611) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(2.3377989) q[0];
rz(-2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-3.0217357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.170341) q[0];
sx q[0];
rz(-2.5345638) q[0];
sx q[0];
rz(-1.5081362) q[0];
x q[1];
rz(0.36977936) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(-1.0880926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2377492) q[1];
sx q[1];
rz(-1.8029873) q[1];
sx q[1];
rz(-3.1234427) q[1];
rz(-pi) q[2];
rz(1.0355789) q[3];
sx q[3];
rz(-2.1018627) q[3];
sx q[3];
rz(0.5862743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5084761) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478093) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(0.48686349) q[0];
rz(-0.72987366) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(1.9015076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86868984) q[0];
sx q[0];
rz(-1.5575952) q[0];
sx q[0];
rz(0.23916434) q[0];
x q[1];
rz(-1.3081595) q[2];
sx q[2];
rz(-1.3381759) q[2];
sx q[2];
rz(-2.5831985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33927321) q[1];
sx q[1];
rz(-0.60802751) q[1];
sx q[1];
rz(0.1165216) q[1];
x q[2];
rz(-2.2206743) q[3];
sx q[3];
rz(-2.4267567) q[3];
sx q[3];
rz(2.9361847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(-0.05274996) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.4809158) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67985095) q[0];
sx q[0];
rz(-1.6438419) q[0];
sx q[0];
rz(-1.9528271) q[0];
rz(-pi) q[1];
rz(1.5691721) q[2];
sx q[2];
rz(-2.6964028) q[2];
sx q[2];
rz(-1.0908529) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.90658014) q[1];
sx q[1];
rz(-1.7791516) q[1];
sx q[1];
rz(2.1085897) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5105641) q[3];
sx q[3];
rz(-2.4619953) q[3];
sx q[3];
rz(-0.70590245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(2.0810614) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(-2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19787191) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(-1.6725756) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(1.3247103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0990471) q[0];
sx q[0];
rz(-1.5481661) q[0];
sx q[0];
rz(3.0142473) q[0];
x q[1];
rz(2.1340738) q[2];
sx q[2];
rz(-2.0118606) q[2];
sx q[2];
rz(1.2457459) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8285671) q[1];
sx q[1];
rz(-1.9777858) q[1];
sx q[1];
rz(-1.160497) q[1];
x q[2];
rz(2.0539829) q[3];
sx q[3];
rz(-0.71648589) q[3];
sx q[3];
rz(-1.0928175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(2.1929072) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(-2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396486) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(-2.6810714) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(-1.2896279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136918) q[0];
sx q[0];
rz(-2.0223589) q[0];
sx q[0];
rz(0.62664647) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3233534) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(0.73275369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21266567) q[1];
sx q[1];
rz(-1.3839098) q[1];
sx q[1];
rz(-3.0287663) q[1];
x q[2];
rz(-1.4054221) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(1.0429494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37725317) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(-2.8073231) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(2.8589378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7455505) q[0];
sx q[0];
rz(-2.1746785) q[0];
sx q[0];
rz(1.0066443) q[0];
x q[1];
rz(-1.7868144) q[2];
sx q[2];
rz(-1.8945433) q[2];
sx q[2];
rz(0.87460364) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.64687318) q[1];
sx q[1];
rz(-2.6962526) q[1];
sx q[1];
rz(0.41782197) q[1];
rz(-pi) q[2];
rz(1.71404) q[3];
sx q[3];
rz(-0.69056615) q[3];
sx q[3];
rz(2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9294372) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(2.5575976) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.042645) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(3.0264061) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(0.68181109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0081351) q[0];
sx q[0];
rz(-1.4032161) q[0];
sx q[0];
rz(2.125678) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.290906) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(2.1249287) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1055923) q[1];
sx q[1];
rz(-0.75165527) q[1];
sx q[1];
rz(2.4113301) q[1];
rz(2.6184611) q[3];
sx q[3];
rz(-2.0683859) q[3];
sx q[3];
rz(-2.350972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8455785) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(-2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(0.63275679) q[2];
sx q[2];
rz(-0.7706332) q[2];
sx q[2];
rz(2.3494233) q[2];
rz(2.6736005) q[3];
sx q[3];
rz(-2.7157126) q[3];
sx q[3];
rz(-1.4117905) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];