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
rz(0.62975878) q[0];
sx q[0];
rz(3.44343) q[0];
sx q[0];
rz(11.259196) q[0];
rz(2.6131926) q[1];
sx q[1];
rz(-2.3101248) q[1];
sx q[1];
rz(2.6822579) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6469) q[0];
sx q[0];
rz(-1.5216869) q[0];
sx q[0];
rz(1.8270593) q[0];
x q[1];
rz(-0.95593837) q[2];
sx q[2];
rz(-0.3345851) q[2];
sx q[2];
rz(-2.0841887) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6428549) q[1];
sx q[1];
rz(-2.5044522) q[1];
sx q[1];
rz(0.29344198) q[1];
rz(-1.5422056) q[3];
sx q[3];
rz(-2.4725879) q[3];
sx q[3];
rz(-1.0278078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1186195) q[2];
sx q[2];
rz(-1.3381253) q[2];
sx q[2];
rz(-0.55707923) q[2];
rz(-2.6856375) q[3];
sx q[3];
rz(-2.3796701) q[3];
sx q[3];
rz(0.60321155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93422455) q[0];
sx q[0];
rz(-2.4226483) q[0];
sx q[0];
rz(1.7508605) q[0];
rz(1.3181744) q[1];
sx q[1];
rz(-1.7134106) q[1];
sx q[1];
rz(0.21600977) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073450449) q[0];
sx q[0];
rz(-1.2552069) q[0];
sx q[0];
rz(0.95629779) q[0];
rz(-pi) q[1];
rz(2.87848) q[2];
sx q[2];
rz(-0.43242681) q[2];
sx q[2];
rz(2.1869834) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5358754) q[1];
sx q[1];
rz(-0.89369666) q[1];
sx q[1];
rz(2.2435921) q[1];
rz(-1.029929) q[3];
sx q[3];
rz(-1.7591624) q[3];
sx q[3];
rz(2.4595367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7184489) q[2];
sx q[2];
rz(-0.37280145) q[2];
sx q[2];
rz(0.049840363) q[2];
rz(-2.5981264) q[3];
sx q[3];
rz(-1.1246357) q[3];
sx q[3];
rz(-1.0928833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.8922888) q[0];
sx q[0];
rz(-1.1019305) q[0];
sx q[0];
rz(-2.1732543) q[0];
rz(3.1210506) q[1];
sx q[1];
rz(-1.3803218) q[1];
sx q[1];
rz(-1.4302018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50043369) q[0];
sx q[0];
rz(-2.4999188) q[0];
sx q[0];
rz(-0.037317567) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5780952) q[2];
sx q[2];
rz(-1.2545956) q[2];
sx q[2];
rz(-1.344156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80950981) q[1];
sx q[1];
rz(-0.44594279) q[1];
sx q[1];
rz(2.6757702) q[1];
rz(-pi) q[2];
rz(-0.084094989) q[3];
sx q[3];
rz(-2.7635592) q[3];
sx q[3];
rz(1.7769906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4358383) q[2];
sx q[2];
rz(-0.64283723) q[2];
sx q[2];
rz(-1.359681) q[2];
rz(-0.5082353) q[3];
sx q[3];
rz(-1.619092) q[3];
sx q[3];
rz(-1.9967509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.085676) q[0];
sx q[0];
rz(-0.29656947) q[0];
sx q[0];
rz(1.0067518) q[0];
rz(-2.309917) q[1];
sx q[1];
rz(-0.74200231) q[1];
sx q[1];
rz(1.1078018) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2058973) q[0];
sx q[0];
rz(-0.95608053) q[0];
sx q[0];
rz(0.45386916) q[0];
rz(-pi) q[1];
rz(-0.57834703) q[2];
sx q[2];
rz(-2.7014534) q[2];
sx q[2];
rz(-2.4348093) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0084787) q[1];
sx q[1];
rz(-1.263887) q[1];
sx q[1];
rz(3.1208193) q[1];
rz(-pi) q[2];
rz(0.44509586) q[3];
sx q[3];
rz(-1.3161049) q[3];
sx q[3];
rz(1.2459918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6882249) q[2];
sx q[2];
rz(-1.0658762) q[2];
sx q[2];
rz(-0.31309703) q[2];
rz(-0.39737663) q[3];
sx q[3];
rz(-1.4704967) q[3];
sx q[3];
rz(0.1853005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6322286) q[0];
sx q[0];
rz(-0.86357421) q[0];
sx q[0];
rz(0.92556992) q[0];
rz(0.46982345) q[1];
sx q[1];
rz(-1.3146105) q[1];
sx q[1];
rz(2.125461) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3811263) q[0];
sx q[0];
rz(-2.010049) q[0];
sx q[0];
rz(-0.45760568) q[0];
rz(0.92570289) q[2];
sx q[2];
rz(-2.7113554) q[2];
sx q[2];
rz(2.1003342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.713607) q[1];
sx q[1];
rz(-1.813289) q[1];
sx q[1];
rz(-2.9338783) q[1];
rz(-pi) q[2];
rz(-0.031948745) q[3];
sx q[3];
rz(-1.2476139) q[3];
sx q[3];
rz(-1.972071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4887345) q[2];
sx q[2];
rz(-0.32565871) q[2];
sx q[2];
rz(2.6540836) q[2];
rz(-0.89753914) q[3];
sx q[3];
rz(-1.2582015) q[3];
sx q[3];
rz(-1.0645617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94051903) q[0];
sx q[0];
rz(-0.93795332) q[0];
sx q[0];
rz(0.67365375) q[0];
rz(1.9081839) q[1];
sx q[1];
rz(-2.1234832) q[1];
sx q[1];
rz(-0.76603755) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8476794) q[0];
sx q[0];
rz(-1.840217) q[0];
sx q[0];
rz(-1.040017) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9658017) q[2];
sx q[2];
rz(-1.5002709) q[2];
sx q[2];
rz(-2.8688237) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5140347) q[1];
sx q[1];
rz(-0.50495353) q[1];
sx q[1];
rz(-1.2673753) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7286639) q[3];
sx q[3];
rz(-2.8097162) q[3];
sx q[3];
rz(-1.6433034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.96887863) q[2];
sx q[2];
rz(-1.6438899) q[2];
sx q[2];
rz(-2.3806351) q[2];
rz(0.40192762) q[3];
sx q[3];
rz(-0.26508489) q[3];
sx q[3];
rz(0.5886122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6589979) q[0];
sx q[0];
rz(-2.3880279) q[0];
sx q[0];
rz(1.6931417) q[0];
rz(-2.6407369) q[1];
sx q[1];
rz(-1.294699) q[1];
sx q[1];
rz(0.81333152) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7194448) q[0];
sx q[0];
rz(-1.1851743) q[0];
sx q[0];
rz(-1.7305806) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7649037) q[2];
sx q[2];
rz(-2.5849205) q[2];
sx q[2];
rz(2.106032) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2709368) q[1];
sx q[1];
rz(-1.7084048) q[1];
sx q[1];
rz(-0.44160053) q[1];
x q[2];
rz(1.1239979) q[3];
sx q[3];
rz(-0.13544336) q[3];
sx q[3];
rz(-2.9970084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27998754) q[2];
sx q[2];
rz(-0.6051175) q[2];
sx q[2];
rz(-2.81847) q[2];
rz(1.4019639) q[3];
sx q[3];
rz(-1.2819382) q[3];
sx q[3];
rz(1.7590947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0793656) q[0];
sx q[0];
rz(-1.0492188) q[0];
sx q[0];
rz(-2.3754689) q[0];
rz(2.3221723) q[1];
sx q[1];
rz(-2.1111646) q[1];
sx q[1];
rz(1.6105509) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19056828) q[0];
sx q[0];
rz(-1.6389209) q[0];
sx q[0];
rz(-2.9124252) q[0];
rz(2.4862637) q[2];
sx q[2];
rz(-1.7458143) q[2];
sx q[2];
rz(-1.9536215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4447915) q[1];
sx q[1];
rz(-0.91233048) q[1];
sx q[1];
rz(2.3347003) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2977912) q[3];
sx q[3];
rz(-0.42359023) q[3];
sx q[3];
rz(-0.18720489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7649585) q[2];
sx q[2];
rz(-2.9464293) q[2];
sx q[2];
rz(-0.39806077) q[2];
rz(0.6063439) q[3];
sx q[3];
rz(-1.5667934) q[3];
sx q[3];
rz(0.91355598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9048318) q[0];
sx q[0];
rz(-2.8396711) q[0];
sx q[0];
rz(-2.6171369) q[0];
rz(0.66871387) q[1];
sx q[1];
rz(-1.5293744) q[1];
sx q[1];
rz(-1.761577) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76261824) q[0];
sx q[0];
rz(-1.6143823) q[0];
sx q[0];
rz(2.144078) q[0];
rz(-2.8369806) q[2];
sx q[2];
rz(-1.1122983) q[2];
sx q[2];
rz(-0.079041399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7126969) q[1];
sx q[1];
rz(-1.98723) q[1];
sx q[1];
rz(-0.45559366) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59137592) q[3];
sx q[3];
rz(-2.5246361) q[3];
sx q[3];
rz(-2.5643831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1606719) q[2];
sx q[2];
rz(-2.394684) q[2];
sx q[2];
rz(2.6778632) q[2];
rz(-1.4240228) q[3];
sx q[3];
rz(-1.0878891) q[3];
sx q[3];
rz(0.033528479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5793107) q[0];
sx q[0];
rz(-0.71826851) q[0];
sx q[0];
rz(-2.723519) q[0];
rz(0.75830165) q[1];
sx q[1];
rz(-0.98379358) q[1];
sx q[1];
rz(-1.2850579) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4649974) q[0];
sx q[0];
rz(-0.9228068) q[0];
sx q[0];
rz(1.6654439) q[0];
x q[1];
rz(0.34512102) q[2];
sx q[2];
rz(-2.1639544) q[2];
sx q[2];
rz(-0.56559222) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27800769) q[1];
sx q[1];
rz(-1.7628432) q[1];
sx q[1];
rz(-0.553756) q[1];
rz(-pi) q[2];
rz(-2.140652) q[3];
sx q[3];
rz(-0.39467282) q[3];
sx q[3];
rz(2.3230524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16984223) q[2];
sx q[2];
rz(-2.1368133) q[2];
sx q[2];
rz(2.8954835) q[2];
rz(-1.2717815) q[3];
sx q[3];
rz(-1.566794) q[3];
sx q[3];
rz(-1.7790214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9825738) q[0];
sx q[0];
rz(-1.6652501) q[0];
sx q[0];
rz(-0.2022947) q[0];
rz(-0.62494878) q[1];
sx q[1];
rz(-2.2849871) q[1];
sx q[1];
rz(-2.7402592) q[1];
rz(-2.305793) q[2];
sx q[2];
rz(-2.8636572) q[2];
sx q[2];
rz(2.5000769) q[2];
rz(-1.2558595) q[3];
sx q[3];
rz(-1.7618014) q[3];
sx q[3];
rz(-0.8509064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
