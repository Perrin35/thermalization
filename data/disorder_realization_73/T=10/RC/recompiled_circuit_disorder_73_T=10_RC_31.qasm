OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(3.6448195) q[0];
sx q[0];
rz(10.148944) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(0.78483265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9871702) q[0];
sx q[0];
rz(-1.3599456) q[0];
sx q[0];
rz(-0.2984557) q[0];
x q[1];
rz(-2.9169693) q[2];
sx q[2];
rz(-2.7135239) q[2];
sx q[2];
rz(-0.12878865) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19476328) q[1];
sx q[1];
rz(-1.0461055) q[1];
sx q[1];
rz(-2.8835433) q[1];
x q[2];
rz(-1.9853398) q[3];
sx q[3];
rz(-1.539955) q[3];
sx q[3];
rz(-2.8950092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.589754) q[2];
sx q[2];
rz(-1.4171615) q[2];
sx q[2];
rz(0.067967728) q[2];
rz(-0.12456482) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(-1.7547866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92000604) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(-2.8979229) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(1.3557281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6512017) q[0];
sx q[0];
rz(-0.25679195) q[0];
sx q[0];
rz(-1.6780361) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7849814) q[2];
sx q[2];
rz(-1.2043673) q[2];
sx q[2];
rz(-2.8474142) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1214952) q[1];
sx q[1];
rz(-2.0165682) q[1];
sx q[1];
rz(-2.9989468) q[1];
x q[2];
rz(1.5242819) q[3];
sx q[3];
rz(-2.1624613) q[3];
sx q[3];
rz(-3.068963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0624861) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-2.8919354) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8963985) q[0];
sx q[0];
rz(-1.9165374) q[0];
sx q[0];
rz(2.2431592) q[0];
rz(-1.8067182) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(-1.2737087) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35667426) q[0];
sx q[0];
rz(-1.2151562) q[0];
sx q[0];
rz(0.26892923) q[0];
x q[1];
rz(1.6689698) q[2];
sx q[2];
rz(-1.2766826) q[2];
sx q[2];
rz(-1.3054747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.18759218) q[1];
sx q[1];
rz(-0.91916579) q[1];
sx q[1];
rz(1.1263532) q[1];
rz(-pi) q[2];
rz(2.0687194) q[3];
sx q[3];
rz(-2.2306799) q[3];
sx q[3];
rz(-2.4592196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0818103) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(-2.188142) q[2];
rz(3.1070784) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(-2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8614486) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(2.8934073) q[0];
rz(-2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(3.0674556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.733487) q[0];
sx q[0];
rz(-1.8251849) q[0];
sx q[0];
rz(1.0803726) q[0];
rz(-pi) q[1];
rz(1.9937236) q[2];
sx q[2];
rz(-2.4175156) q[2];
sx q[2];
rz(-2.1507182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6756235) q[1];
sx q[1];
rz(-0.33756653) q[1];
sx q[1];
rz(1.2077431) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3965963) q[3];
sx q[3];
rz(-1.8074236) q[3];
sx q[3];
rz(-3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8903824) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(3.1029491) q[2];
rz(-2.1679227) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(1.779153) q[0];
rz(-0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.978925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29279941) q[0];
sx q[0];
rz(-1.6822364) q[0];
sx q[0];
rz(3.1116027) q[0];
x q[1];
rz(1.0983724) q[2];
sx q[2];
rz(-2.782151) q[2];
sx q[2];
rz(0.098509468) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8414383) q[1];
sx q[1];
rz(-1.9342124) q[1];
sx q[1];
rz(-0.20519786) q[1];
rz(-pi) q[2];
rz(-1.053247) q[3];
sx q[3];
rz(-1.9540817) q[3];
sx q[3];
rz(2.7305207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(0.14222063) q[2];
rz(2.2375315) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(0.18946762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.6739155) q[0];
rz(2.5698075) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(-2.8335559) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338617) q[0];
sx q[0];
rz(-2.9836285) q[0];
sx q[0];
rz(-0.64361848) q[0];
x q[1];
rz(0.96254827) q[2];
sx q[2];
rz(-1.0481917) q[2];
sx q[2];
rz(-2.7163497) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8923556) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(-1.874079) q[1];
x q[2];
rz(-0.24758731) q[3];
sx q[3];
rz(-0.35841225) q[3];
sx q[3];
rz(-0.7022411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3623111) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(-2.4273196) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(-2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(-0.1299783) q[0];
rz(-3.1107483) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-2.470509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111611) q[0];
sx q[0];
rz(-1.5357619) q[0];
sx q[0];
rz(3.0879211) q[0];
rz(3.1110711) q[2];
sx q[2];
rz(-1.3866716) q[2];
sx q[2];
rz(2.7437999) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5741201) q[1];
sx q[1];
rz(-1.6735055) q[1];
sx q[1];
rz(-1.1979539) q[1];
rz(-0.89550771) q[3];
sx q[3];
rz(-1.5338147) q[3];
sx q[3];
rz(1.0469588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.122763) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(3.1130062) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(-1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(2.7291765) q[0];
rz(1.4498129) q[1];
sx q[1];
rz(-1.342536) q[1];
sx q[1];
rz(-1.1669881) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1056571) q[0];
sx q[0];
rz(-1.3577537) q[0];
sx q[0];
rz(-2.1980594) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86161676) q[2];
sx q[2];
rz(-1.3590727) q[2];
sx q[2];
rz(-0.35702969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2003157) q[1];
sx q[1];
rz(-1.3822767) q[1];
sx q[1];
rz(1.6016866) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9917631) q[3];
sx q[3];
rz(-2.095788) q[3];
sx q[3];
rz(-2.6168952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.940544) q[2];
sx q[2];
rz(-0.87934914) q[2];
sx q[2];
rz(2.9525625) q[2];
rz(2.9947301) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(-1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(0.92700672) q[0];
rz(-1.758763) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(-1.4896726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4179174) q[0];
sx q[0];
rz(-2.0360887) q[0];
sx q[0];
rz(-2.377541) q[0];
x q[1];
rz(1.0987368) q[2];
sx q[2];
rz(-1.615881) q[2];
sx q[2];
rz(2.2301607) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5786963) q[1];
sx q[1];
rz(-1.7637196) q[1];
sx q[1];
rz(-1.7460515) q[1];
rz(-2.6418266) q[3];
sx q[3];
rz(-1.0904113) q[3];
sx q[3];
rz(1.2479316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(-1.4302953) q[2];
rz(0.57724214) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(-1.757471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5532613) q[0];
sx q[0];
rz(-1.7734779) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-2.9945701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518236) q[0];
sx q[0];
rz(-0.71634403) q[0];
sx q[0];
rz(2.1728974) q[0];
rz(-pi) q[1];
rz(-2.8617919) q[2];
sx q[2];
rz(-0.24926148) q[2];
sx q[2];
rz(1.3485497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30565572) q[1];
sx q[1];
rz(-0.61969295) q[1];
sx q[1];
rz(1.5020919) q[1];
rz(-pi) q[2];
rz(2.9045194) q[3];
sx q[3];
rz(-0.58963886) q[3];
sx q[3];
rz(1.3414563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0599351) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(0.74404136) q[2];
rz(-2.384281) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1159146) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(-0.81746447) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(-1.6384009) q[2];
sx q[2];
rz(-2.0220145) q[2];
sx q[2];
rz(1.8215712) q[2];
rz(3.01132) q[3];
sx q[3];
rz(-1.0126922) q[3];
sx q[3];
rz(-1.0425413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
