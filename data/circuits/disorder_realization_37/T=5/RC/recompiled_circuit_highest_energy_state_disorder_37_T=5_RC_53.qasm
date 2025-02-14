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
rz(-0.77684075) q[0];
sx q[0];
rz(-0.87640327) q[0];
sx q[0];
rz(-0.25282282) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(4.0568772) q[1];
sx q[1];
rz(8.6203909) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0067449) q[0];
sx q[0];
rz(-0.50323707) q[0];
sx q[0];
rz(-0.87849599) q[0];
x q[1];
rz(2.6635799) q[2];
sx q[2];
rz(-0.4768663) q[2];
sx q[2];
rz(3.06682) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8411257) q[1];
sx q[1];
rz(-1.2004832) q[1];
sx q[1];
rz(-1.200586) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7976296) q[3];
sx q[3];
rz(-1.6795782) q[3];
sx q[3];
rz(-3.1358842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.09314166) q[2];
sx q[2];
rz(-0.96258771) q[2];
sx q[2];
rz(-2.2691881) q[2];
rz(-0.33556542) q[3];
sx q[3];
rz(-1.7458956) q[3];
sx q[3];
rz(3.057835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5035079) q[0];
sx q[0];
rz(-0.26569772) q[0];
sx q[0];
rz(0.22915325) q[0];
rz(-2.9511662) q[1];
sx q[1];
rz(-1.3399905) q[1];
sx q[1];
rz(-2.4280039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.021727) q[0];
sx q[0];
rz(-1.83868) q[0];
sx q[0];
rz(-0.88587205) q[0];
rz(-pi) q[1];
rz(-2.607478) q[2];
sx q[2];
rz(-0.91098173) q[2];
sx q[2];
rz(2.5995863) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.16878) q[1];
sx q[1];
rz(-2.3832364) q[1];
sx q[1];
rz(-0.43557628) q[1];
rz(-1.7341033) q[3];
sx q[3];
rz(-2.3288832) q[3];
sx q[3];
rz(-0.14369609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4552292) q[2];
sx q[2];
rz(-0.31193048) q[2];
sx q[2];
rz(-2.6658106) q[2];
rz(-2.0717715) q[3];
sx q[3];
rz(-2.1333623) q[3];
sx q[3];
rz(-2.3750335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0692724) q[0];
sx q[0];
rz(-1.9444436) q[0];
sx q[0];
rz(-1.9092165) q[0];
rz(0.45066372) q[1];
sx q[1];
rz(-1.181299) q[1];
sx q[1];
rz(0.88952363) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8101533) q[0];
sx q[0];
rz(-2.4001849) q[0];
sx q[0];
rz(-2.5196694) q[0];
rz(-pi) q[1];
rz(2.6085147) q[2];
sx q[2];
rz(-2.0642082) q[2];
sx q[2];
rz(0.62668884) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7550068) q[1];
sx q[1];
rz(-2.4535547) q[1];
sx q[1];
rz(-1.2351456) q[1];
rz(-0.74340762) q[3];
sx q[3];
rz(-0.49170676) q[3];
sx q[3];
rz(1.7984901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4543317) q[2];
sx q[2];
rz(-2.7165518) q[2];
sx q[2];
rz(1.8827776) q[2];
rz(-1.2339633) q[3];
sx q[3];
rz(-1.1326658) q[3];
sx q[3];
rz(-3.1335355) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71037978) q[0];
sx q[0];
rz(-2.2358535) q[0];
sx q[0];
rz(1.6718965) q[0];
rz(1.1471033) q[1];
sx q[1];
rz(-1.0018145) q[1];
sx q[1];
rz(2.240644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1484447) q[0];
sx q[0];
rz(-1.6560418) q[0];
sx q[0];
rz(0.31384018) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1828762) q[2];
sx q[2];
rz(-2.0953853) q[2];
sx q[2];
rz(-2.1647349) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1114677) q[1];
sx q[1];
rz(-0.84967597) q[1];
sx q[1];
rz(-1.6961873) q[1];
rz(-2.9649686) q[3];
sx q[3];
rz(-2.400853) q[3];
sx q[3];
rz(2.3185042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6441696) q[2];
sx q[2];
rz(-1.3966509) q[2];
sx q[2];
rz(0.72745848) q[2];
rz(2.4265031) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(-0.49811825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46399507) q[0];
sx q[0];
rz(-1.7514739) q[0];
sx q[0];
rz(-2.4053307) q[0];
rz(-0.62034208) q[1];
sx q[1];
rz(-0.54519975) q[1];
sx q[1];
rz(1.6166519) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0686091) q[0];
sx q[0];
rz(-1.6260379) q[0];
sx q[0];
rz(-1.5609571) q[0];
rz(-pi) q[1];
rz(0.23133607) q[2];
sx q[2];
rz(-1.243719) q[2];
sx q[2];
rz(-1.2755204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.273539) q[1];
sx q[1];
rz(-2.2333849) q[1];
sx q[1];
rz(2.6602547) q[1];
x q[2];
rz(-2.3297263) q[3];
sx q[3];
rz(-0.97211876) q[3];
sx q[3];
rz(0.52385274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5357431) q[2];
sx q[2];
rz(-0.77631408) q[2];
sx q[2];
rz(-2.1898451) q[2];
rz(-0.4176628) q[3];
sx q[3];
rz(-2.8916736) q[3];
sx q[3];
rz(-1.1550268) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34053892) q[0];
sx q[0];
rz(-1.8908353) q[0];
sx q[0];
rz(2.387555) q[0];
rz(-2.2484089) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(2.975614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8360031) q[0];
sx q[0];
rz(-0.89035946) q[0];
sx q[0];
rz(-1.7050118) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31678172) q[2];
sx q[2];
rz(-0.52826476) q[2];
sx q[2];
rz(1.007391) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41478911) q[1];
sx q[1];
rz(-0.94894743) q[1];
sx q[1];
rz(-1.0110823) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13940553) q[3];
sx q[3];
rz(-1.8434815) q[3];
sx q[3];
rz(2.8806339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7665427) q[2];
sx q[2];
rz(-2.0039717) q[2];
sx q[2];
rz(0.005391187) q[2];
rz(3.052875) q[3];
sx q[3];
rz(-1.6088477) q[3];
sx q[3];
rz(-2.3458792) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8832815) q[0];
sx q[0];
rz(-1.9323876) q[0];
sx q[0];
rz(-0.2051556) q[0];
rz(-0.908665) q[1];
sx q[1];
rz(-0.87037194) q[1];
sx q[1];
rz(0.044513449) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9712898) q[0];
sx q[0];
rz(-1.7006386) q[0];
sx q[0];
rz(1.4128039) q[0];
x q[1];
rz(-1.7111383) q[2];
sx q[2];
rz(-1.5333041) q[2];
sx q[2];
rz(0.692779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2878814) q[1];
sx q[1];
rz(-1.3559623) q[1];
sx q[1];
rz(0.98835215) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3839339) q[3];
sx q[3];
rz(-3.0198041) q[3];
sx q[3];
rz(0.78021061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.32555106) q[2];
sx q[2];
rz(-2.6191235) q[2];
sx q[2];
rz(-0.52118707) q[2];
rz(-0.19736396) q[3];
sx q[3];
rz(-1.7557996) q[3];
sx q[3];
rz(0.50281966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7407783) q[0];
sx q[0];
rz(-2.9476808) q[0];
sx q[0];
rz(-0.25522301) q[0];
rz(0.1420282) q[1];
sx q[1];
rz(-1.4165001) q[1];
sx q[1];
rz(2.7878888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86310327) q[0];
sx q[0];
rz(-1.2444087) q[0];
sx q[0];
rz(-1.9410237) q[0];
x q[1];
rz(2.1944517) q[2];
sx q[2];
rz(-1.9883687) q[2];
sx q[2];
rz(3.0018978) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4507915) q[1];
sx q[1];
rz(-2.004068) q[1];
sx q[1];
rz(2.7669897) q[1];
x q[2];
rz(-2.1347624) q[3];
sx q[3];
rz(-2.2368397) q[3];
sx q[3];
rz(0.49426916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3069309) q[2];
sx q[2];
rz(-2.8262409) q[2];
sx q[2];
rz(1.5241415) q[2];
rz(0.8664242) q[3];
sx q[3];
rz(-1.4640936) q[3];
sx q[3];
rz(0.58969897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6397112) q[0];
sx q[0];
rz(-2.9852133) q[0];
sx q[0];
rz(-1.7891275) q[0];
rz(1.9068708) q[1];
sx q[1];
rz(-2.9094978) q[1];
sx q[1];
rz(-2.629705) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8361266) q[0];
sx q[0];
rz(-2.9752467) q[0];
sx q[0];
rz(-1.8977099) q[0];
rz(-pi) q[1];
rz(2.9779766) q[2];
sx q[2];
rz(-2.3694042) q[2];
sx q[2];
rz(0.26547394) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33495397) q[1];
sx q[1];
rz(-0.52577621) q[1];
sx q[1];
rz(-0.61383617) q[1];
rz(-pi) q[2];
rz(2.0824349) q[3];
sx q[3];
rz(-0.64877227) q[3];
sx q[3];
rz(1.8442475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1341683) q[2];
sx q[2];
rz(-1.8210501) q[2];
sx q[2];
rz(2.7443938) q[2];
rz(2.126179) q[3];
sx q[3];
rz(-1.0011324) q[3];
sx q[3];
rz(2.5889682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7427202) q[0];
sx q[0];
rz(-1.2074559) q[0];
sx q[0];
rz(-2.7040828) q[0];
rz(0.73879009) q[1];
sx q[1];
rz(-0.80016017) q[1];
sx q[1];
rz(-1.4791666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9112944) q[0];
sx q[0];
rz(-2.0574967) q[0];
sx q[0];
rz(-1.6337956) q[0];
rz(-pi) q[1];
rz(-0.11585856) q[2];
sx q[2];
rz(-1.2136974) q[2];
sx q[2];
rz(-3.0810205) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91875118) q[1];
sx q[1];
rz(-1.4645618) q[1];
sx q[1];
rz(-1.5083471) q[1];
rz(-pi) q[2];
rz(0.78082943) q[3];
sx q[3];
rz(-2.1159647) q[3];
sx q[3];
rz(0.050267537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9208357) q[2];
sx q[2];
rz(-1.7831384) q[2];
sx q[2];
rz(2.9662568) q[2];
rz(-2.1389424) q[3];
sx q[3];
rz(-2.9084539) q[3];
sx q[3];
rz(1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.47068448) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(0.15057527) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(-0.19828037) q[2];
sx q[2];
rz(-1.5783327) q[2];
sx q[2];
rz(-0.21680149) q[2];
rz(-0.2507052) q[3];
sx q[3];
rz(-1.7964994) q[3];
sx q[3];
rz(-2.9356706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
