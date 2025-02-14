OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98439944) q[0];
sx q[0];
rz(-1.1097044) q[0];
sx q[0];
rz(-2.1362526) q[0];
rz(0.73468626) q[1];
sx q[1];
rz(-2.7856196) q[1];
sx q[1];
rz(2.2495143) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9933424) q[0];
sx q[0];
rz(-2.3138232) q[0];
sx q[0];
rz(-1.4367075) q[0];
x q[1];
rz(1.7847127) q[2];
sx q[2];
rz(-2.5271673) q[2];
sx q[2];
rz(0.015903552) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2764276) q[1];
sx q[1];
rz(-0.49404198) q[1];
sx q[1];
rz(0.20706351) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69697505) q[3];
sx q[3];
rz(-2.4066952) q[3];
sx q[3];
rz(-2.1293726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4165108) q[2];
sx q[2];
rz(-1.9981013) q[2];
sx q[2];
rz(-1.4408646) q[2];
rz(-0.95669389) q[3];
sx q[3];
rz(-1.1172833) q[3];
sx q[3];
rz(-2.8982437) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44593921) q[0];
sx q[0];
rz(-0.39458269) q[0];
sx q[0];
rz(-2.0126427) q[0];
rz(0.24540643) q[1];
sx q[1];
rz(-1.3134198) q[1];
sx q[1];
rz(2.479877) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6889383) q[0];
sx q[0];
rz(-2.1784349) q[0];
sx q[0];
rz(1.9197965) q[0];
rz(-pi) q[1];
rz(0.38208802) q[2];
sx q[2];
rz(-1.3616865) q[2];
sx q[2];
rz(1.2651625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1744888) q[1];
sx q[1];
rz(-1.8421108) q[1];
sx q[1];
rz(-1.6037462) q[1];
rz(-pi) q[2];
rz(0.097019086) q[3];
sx q[3];
rz(-1.8604014) q[3];
sx q[3];
rz(1.1821018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5439862) q[2];
sx q[2];
rz(-0.49763766) q[2];
sx q[2];
rz(1.0428766) q[2];
rz(1.6889702) q[3];
sx q[3];
rz(-2.2904604) q[3];
sx q[3];
rz(-2.0378621) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747028) q[0];
sx q[0];
rz(-2.4536528) q[0];
sx q[0];
rz(-1.8324628) q[0];
rz(2.6922928) q[1];
sx q[1];
rz(-2.4520912) q[1];
sx q[1];
rz(-1.4264533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4065789) q[0];
sx q[0];
rz(-1.8157703) q[0];
sx q[0];
rz(1.2981775) q[0];
x q[1];
rz(-3.1235891) q[2];
sx q[2];
rz(-1.1483542) q[2];
sx q[2];
rz(1.2382335) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.215308) q[1];
sx q[1];
rz(-0.70901543) q[1];
sx q[1];
rz(0.11911094) q[1];
rz(-pi) q[2];
rz(-1.2481232) q[3];
sx q[3];
rz(-0.25882684) q[3];
sx q[3];
rz(0.62773529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36491498) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9471112) q[0];
sx q[0];
rz(-2.8567061) q[0];
sx q[0];
rz(-0.77199212) q[0];
rz(2.4783065) q[1];
sx q[1];
rz(-2.3978077) q[1];
sx q[1];
rz(-0.1276806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6948815) q[0];
sx q[0];
rz(-2.237473) q[0];
sx q[0];
rz(1.5645909) q[0];
x q[1];
rz(-2.588932) q[2];
sx q[2];
rz(-1.895012) q[2];
sx q[2];
rz(-0.49190258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7230715) q[1];
sx q[1];
rz(-1.1643895) q[1];
sx q[1];
rz(-1.0999098) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4127619) q[3];
sx q[3];
rz(-1.1452066) q[3];
sx q[3];
rz(-1.1273718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6105303) q[2];
sx q[2];
rz(-2.1805111) q[2];
sx q[2];
rz(-2.5430211) q[2];
rz(0.5851723) q[3];
sx q[3];
rz(-1.7681311) q[3];
sx q[3];
rz(3.0857871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-2.1685261) q[0];
rz(2.9452501) q[1];
sx q[1];
rz(-2.0025496) q[1];
sx q[1];
rz(-1.6729209) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7356883) q[0];
sx q[0];
rz(-2.2668512) q[0];
sx q[0];
rz(2.619834) q[0];
rz(-2.8873047) q[2];
sx q[2];
rz(-1.0788222) q[2];
sx q[2];
rz(2.6300583) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65485146) q[1];
sx q[1];
rz(-1.2706725) q[1];
sx q[1];
rz(-2.9611118) q[1];
x q[2];
rz(-2.1096538) q[3];
sx q[3];
rz(-1.5232289) q[3];
sx q[3];
rz(0.69456929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68957442) q[2];
sx q[2];
rz(-1.3268027) q[2];
sx q[2];
rz(-0.16239521) q[2];
rz(-1.1299805) q[3];
sx q[3];
rz(-0.88310784) q[3];
sx q[3];
rz(-2.0067154) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75055403) q[0];
sx q[0];
rz(-0.51893187) q[0];
sx q[0];
rz(0.050405141) q[0];
rz(-0.16818908) q[1];
sx q[1];
rz(-1.5610118) q[1];
sx q[1];
rz(2.2680297) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03913232) q[0];
sx q[0];
rz(-2.1796569) q[0];
sx q[0];
rz(-1.2723421) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4292154) q[2];
sx q[2];
rz(-1.537286) q[2];
sx q[2];
rz(-0.17591116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.62100121) q[1];
sx q[1];
rz(-1.1258003) q[1];
sx q[1];
rz(0.59153647) q[1];
x q[2];
rz(-1.3662947) q[3];
sx q[3];
rz(-1.8629774) q[3];
sx q[3];
rz(0.99537163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1217338) q[2];
sx q[2];
rz(-0.21616082) q[2];
sx q[2];
rz(2.7609694) q[2];
rz(2.5753283) q[3];
sx q[3];
rz(-1.7411313) q[3];
sx q[3];
rz(-2.3316021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009163) q[0];
sx q[0];
rz(-0.2114507) q[0];
sx q[0];
rz(2.7389615) q[0];
rz(1.3817878) q[1];
sx q[1];
rz(-2.333162) q[1];
sx q[1];
rz(-0.056338739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23898174) q[0];
sx q[0];
rz(-1.9704116) q[0];
sx q[0];
rz(-0.35524551) q[0];
rz(3.10059) q[2];
sx q[2];
rz(-1.1004538) q[2];
sx q[2];
rz(0.13861632) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4693335) q[1];
sx q[1];
rz(-2.1348828) q[1];
sx q[1];
rz(3.1271598) q[1];
rz(2.5272156) q[3];
sx q[3];
rz(-2.3639332) q[3];
sx q[3];
rz(-1.3953502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9075883) q[2];
sx q[2];
rz(-1.3434709) q[2];
sx q[2];
rz(-0.50055707) q[2];
rz(0.060700011) q[3];
sx q[3];
rz(-1.7874291) q[3];
sx q[3];
rz(2.6589656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57182264) q[0];
sx q[0];
rz(-1.7800542) q[0];
sx q[0];
rz(-1.7806336) q[0];
rz(2.3979777) q[1];
sx q[1];
rz(-1.7230325) q[1];
sx q[1];
rz(-3.0027711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69037914) q[0];
sx q[0];
rz(-2.0669961) q[0];
sx q[0];
rz(-0.98753937) q[0];
rz(0.46625579) q[2];
sx q[2];
rz(-1.0276349) q[2];
sx q[2];
rz(-2.0838497) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.946859) q[1];
sx q[1];
rz(-1.3689965) q[1];
sx q[1];
rz(1.9029593) q[1];
rz(-1.9622063) q[3];
sx q[3];
rz(-1.7633121) q[3];
sx q[3];
rz(-0.26993902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6737785) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(0.71411258) q[2];
rz(-2.1821187) q[3];
sx q[3];
rz(-1.4941314) q[3];
sx q[3];
rz(-0.97904557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7239083) q[0];
sx q[0];
rz(-2.4651616) q[0];
sx q[0];
rz(-0.12829256) q[0];
rz(0.3872321) q[1];
sx q[1];
rz(-2.5435244) q[1];
sx q[1];
rz(1.8181575) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.58777) q[0];
sx q[0];
rz(-2.9284796) q[0];
sx q[0];
rz(2.0487259) q[0];
rz(-pi) q[1];
rz(2.2567184) q[2];
sx q[2];
rz(-1.0909683) q[2];
sx q[2];
rz(2.4201833) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9518428) q[1];
sx q[1];
rz(-1.4268759) q[1];
sx q[1];
rz(-0.014149498) q[1];
rz(1.3714482) q[3];
sx q[3];
rz(-2.2888765) q[3];
sx q[3];
rz(1.1465286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6748176) q[2];
sx q[2];
rz(-1.6206436) q[2];
sx q[2];
rz(0.20469323) q[2];
rz(1.2734867) q[3];
sx q[3];
rz(-0.73015648) q[3];
sx q[3];
rz(-1.5654806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650763) q[0];
sx q[0];
rz(-2.5193494) q[0];
sx q[0];
rz(2.9283071) q[0];
rz(0.22008303) q[1];
sx q[1];
rz(-1.2525696) q[1];
sx q[1];
rz(-1.3594886) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78173554) q[0];
sx q[0];
rz(-3.1229501) q[0];
sx q[0];
rz(-2.1117979) q[0];
rz(0.47021265) q[2];
sx q[2];
rz(-2.6192952) q[2];
sx q[2];
rz(0.14005113) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9656187) q[1];
sx q[1];
rz(-2.5104264) q[1];
sx q[1];
rz(2.1793773) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1207934) q[3];
sx q[3];
rz(-0.74899835) q[3];
sx q[3];
rz(-1.7089546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.517211) q[2];
sx q[2];
rz(-1.5456079) q[2];
sx q[2];
rz(-2.8037996) q[2];
rz(1.4126011) q[3];
sx q[3];
rz(-1.1510886) q[3];
sx q[3];
rz(0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4008041) q[0];
sx q[0];
rz(-2.3176226) q[0];
sx q[0];
rz(1.2299706) q[0];
rz(0.38372718) q[1];
sx q[1];
rz(-1.4839254) q[1];
sx q[1];
rz(0.76212777) q[1];
rz(2.7376851) q[2];
sx q[2];
rz(-1.8941034) q[2];
sx q[2];
rz(-0.20878172) q[2];
rz(-2.8970412) q[3];
sx q[3];
rz(-2.0691732) q[3];
sx q[3];
rz(0.3680784) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
