OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1571932) q[0];
sx q[0];
rz(-2.0318883) q[0];
sx q[0];
rz(-1.00534) q[0];
rz(0.73468626) q[1];
sx q[1];
rz(3.4975657) q[1];
sx q[1];
rz(11.674292) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9514314) q[0];
sx q[0];
rz(-2.3888458) q[0];
sx q[0];
rz(-2.997082) q[0];
x q[1];
rz(2.1744035) q[2];
sx q[2];
rz(-1.4481067) q[2];
sx q[2];
rz(-1.762378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0301117) q[1];
sx q[1];
rz(-1.6684384) q[1];
sx q[1];
rz(0.48508118) q[1];
rz(-pi) q[2];
rz(2.0964811) q[3];
sx q[3];
rz(-2.1108004) q[3];
sx q[3];
rz(-1.857855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72508183) q[2];
sx q[2];
rz(-1.9981013) q[2];
sx q[2];
rz(1.4408646) q[2];
rz(0.95669389) q[3];
sx q[3];
rz(-1.1172833) q[3];
sx q[3];
rz(2.8982437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6956534) q[0];
sx q[0];
rz(-0.39458269) q[0];
sx q[0];
rz(1.12895) q[0];
rz(-2.8961862) q[1];
sx q[1];
rz(-1.3134198) q[1];
sx q[1];
rz(2.479877) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.322987) q[0];
sx q[0];
rz(-1.855369) q[0];
sx q[0];
rz(2.5045082) q[0];
rz(-pi) q[1];
rz(1.7956338) q[2];
sx q[2];
rz(-1.944146) q[2];
sx q[2];
rz(-0.38885798) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0894818) q[1];
sx q[1];
rz(-0.27325892) q[1];
sx q[1];
rz(-0.11788003) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2798987) q[3];
sx q[3];
rz(-1.4778294) q[3];
sx q[3];
rz(-2.7806843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5976065) q[2];
sx q[2];
rz(-2.643955) q[2];
sx q[2];
rz(-1.0428766) q[2];
rz(-1.6889702) q[3];
sx q[3];
rz(-0.85113227) q[3];
sx q[3];
rz(-2.0378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76688981) q[0];
sx q[0];
rz(-2.4536528) q[0];
sx q[0];
rz(1.3091298) q[0];
rz(-0.4492999) q[1];
sx q[1];
rz(-2.4520912) q[1];
sx q[1];
rz(1.7151394) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5912522) q[0];
sx q[0];
rz(-0.36446291) q[0];
sx q[0];
rz(-0.82247199) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9932991) q[2];
sx q[2];
rz(-1.5543756) q[2];
sx q[2];
rz(2.8016479) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.215308) q[1];
sx q[1];
rz(-2.4325772) q[1];
sx q[1];
rz(3.0224817) q[1];
rz(1.2481232) q[3];
sx q[3];
rz(-2.8827658) q[3];
sx q[3];
rz(0.62773529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7766777) q[2];
sx q[2];
rz(-1.1740351) q[2];
sx q[2];
rz(-0.066369973) q[2];
rz(2.5868609) q[3];
sx q[3];
rz(-1.4812255) q[3];
sx q[3];
rz(1.6248645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1944815) q[0];
sx q[0];
rz(-0.28488657) q[0];
sx q[0];
rz(-2.3696005) q[0];
rz(-2.4783065) q[1];
sx q[1];
rz(-2.3978077) q[1];
sx q[1];
rz(0.1276806) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4467112) q[0];
sx q[0];
rz(-2.237473) q[0];
sx q[0];
rz(1.5645909) q[0];
x q[1];
rz(1.1947317) q[2];
sx q[2];
rz(-2.0916307) q[2];
sx q[2];
rz(-1.8686881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7230715) q[1];
sx q[1];
rz(-1.9772031) q[1];
sx q[1];
rz(-1.0999098) q[1];
rz(0.43032448) q[3];
sx q[3];
rz(-1.4269623) q[3];
sx q[3];
rz(0.50912428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53106236) q[2];
sx q[2];
rz(-2.1805111) q[2];
sx q[2];
rz(2.5430211) q[2];
rz(-2.5564204) q[3];
sx q[3];
rz(-1.3734615) q[3];
sx q[3];
rz(0.055805512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7912306) q[0];
sx q[0];
rz(-1.5086011) q[0];
sx q[0];
rz(0.97306657) q[0];
rz(-2.9452501) q[1];
sx q[1];
rz(-2.0025496) q[1];
sx q[1];
rz(1.6729209) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40590433) q[0];
sx q[0];
rz(-0.87474147) q[0];
sx q[0];
rz(-0.52175867) q[0];
rz(-pi) q[1];
rz(-2.009654) q[2];
sx q[2];
rz(-0.54897749) q[2];
sx q[2];
rz(-1.0145607) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.10188516) q[1];
sx q[1];
rz(-2.7927924) q[1];
sx q[1];
rz(1.0451911) q[1];
rz(-pi) q[2];
rz(1.4782952) q[3];
sx q[3];
rz(-2.6008462) q[3];
sx q[3];
rz(-2.1860022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68957442) q[2];
sx q[2];
rz(-1.3268027) q[2];
sx q[2];
rz(-0.16239521) q[2];
rz(1.1299805) q[3];
sx q[3];
rz(-0.88310784) q[3];
sx q[3];
rz(-1.1348772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3910386) q[0];
sx q[0];
rz(-2.6226608) q[0];
sx q[0];
rz(-0.050405141) q[0];
rz(-2.9734036) q[1];
sx q[1];
rz(-1.5610118) q[1];
sx q[1];
rz(0.87356299) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024603) q[0];
sx q[0];
rz(-0.96193571) q[0];
sx q[0];
rz(1.2723421) q[0];
rz(-pi) q[1];
rz(-0.033848714) q[2];
sx q[2];
rz(-1.7122972) q[2];
sx q[2];
rz(1.7514829) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.37994036) q[1];
sx q[1];
rz(-0.72387513) q[1];
sx q[1];
rz(2.4340802) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29804067) q[3];
sx q[3];
rz(-1.7665157) q[3];
sx q[3];
rz(2.6258371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1217338) q[2];
sx q[2];
rz(-0.21616082) q[2];
sx q[2];
rz(0.38062322) q[2];
rz(2.5753283) q[3];
sx q[3];
rz(-1.4004613) q[3];
sx q[3];
rz(2.3316021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64067632) q[0];
sx q[0];
rz(-2.930142) q[0];
sx q[0];
rz(2.7389615) q[0];
rz(1.3817878) q[1];
sx q[1];
rz(-0.80843061) q[1];
sx q[1];
rz(-3.0852539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9531247) q[0];
sx q[0];
rz(-1.896965) q[0];
sx q[0];
rz(1.1475546) q[0];
rz(-pi) q[1];
rz(2.0414786) q[2];
sx q[2];
rz(-1.5342481) q[2];
sx q[2];
rz(1.728003) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64526865) q[1];
sx q[1];
rz(-2.5773415) q[1];
sx q[1];
rz(-1.5936046) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67752892) q[3];
sx q[3];
rz(-1.9871681) q[3];
sx q[3];
rz(0.29034055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9075883) q[2];
sx q[2];
rz(-1.3434709) q[2];
sx q[2];
rz(2.6410356) q[2];
rz(-0.060700011) q[3];
sx q[3];
rz(-1.3541636) q[3];
sx q[3];
rz(-0.48262706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.56977) q[0];
sx q[0];
rz(-1.7800542) q[0];
sx q[0];
rz(1.7806336) q[0];
rz(-0.743615) q[1];
sx q[1];
rz(-1.7230325) q[1];
sx q[1];
rz(0.13882151) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4512135) q[0];
sx q[0];
rz(-1.0745966) q[0];
sx q[0];
rz(-0.98753937) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93075625) q[2];
sx q[2];
rz(-0.70022196) q[2];
sx q[2];
rz(-2.8559358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.19473361) q[1];
sx q[1];
rz(-1.7725962) q[1];
sx q[1];
rz(-1.2386333) q[1];
x q[2];
rz(1.0984135) q[3];
sx q[3];
rz(-0.43398991) q[3];
sx q[3];
rz(-1.7349617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6737785) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(0.71411258) q[2];
rz(2.1821187) q[3];
sx q[3];
rz(-1.6474612) q[3];
sx q[3];
rz(-0.97904557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4176843) q[0];
sx q[0];
rz(-2.4651616) q[0];
sx q[0];
rz(-0.12829256) q[0];
rz(2.7543606) q[1];
sx q[1];
rz(-0.59806824) q[1];
sx q[1];
rz(-1.3234352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6559728) q[0];
sx q[0];
rz(-1.6682297) q[0];
sx q[0];
rz(-1.3809588) q[0];
rz(-pi) q[1];
rz(-0.88487423) q[2];
sx q[2];
rz(-1.0909683) q[2];
sx q[2];
rz(2.4201833) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7625757) q[1];
sx q[1];
rz(-1.5847995) q[1];
sx q[1];
rz(1.4268616) q[1];
rz(1.3714482) q[3];
sx q[3];
rz(-2.2888765) q[3];
sx q[3];
rz(1.1465286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6748176) q[2];
sx q[2];
rz(-1.6206436) q[2];
sx q[2];
rz(-2.9368994) q[2];
rz(1.2734867) q[3];
sx q[3];
rz(-0.73015648) q[3];
sx q[3];
rz(-1.5654806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765163) q[0];
sx q[0];
rz(-0.62224329) q[0];
sx q[0];
rz(2.9283071) q[0];
rz(-2.9215096) q[1];
sx q[1];
rz(-1.2525696) q[1];
sx q[1];
rz(-1.3594886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8187788) q[0];
sx q[0];
rz(-1.5548163) q[0];
sx q[0];
rz(3.131991) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.67138) q[2];
sx q[2];
rz(-0.5222975) q[2];
sx q[2];
rz(3.0015415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90725856) q[1];
sx q[1];
rz(-1.2266907) q[1];
sx q[1];
rz(-2.1109778) q[1];
x q[2];
rz(-0.45234802) q[3];
sx q[3];
rz(-0.95146639) q[3];
sx q[3];
rz(-2.1296219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6243817) q[2];
sx q[2];
rz(-1.5959847) q[2];
sx q[2];
rz(2.8037996) q[2];
rz(1.7289915) q[3];
sx q[3];
rz(-1.9905041) q[3];
sx q[3];
rz(0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74078858) q[0];
sx q[0];
rz(-0.82397006) q[0];
sx q[0];
rz(-1.9116221) q[0];
rz(0.38372718) q[1];
sx q[1];
rz(-1.4839254) q[1];
sx q[1];
rz(0.76212777) q[1];
rz(2.4356213) q[2];
sx q[2];
rz(-2.6298475) q[2];
sx q[2];
rz(2.0012326) q[2];
rz(-0.2445515) q[3];
sx q[3];
rz(-1.0724194) q[3];
sx q[3];
rz(-2.7735143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
