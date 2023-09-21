OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(3.9324023) q[0];
sx q[0];
rz(12.232236) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5028635) q[0];
sx q[0];
rz(-1.6067061) q[0];
sx q[0];
rz(2.3907651) q[0];
rz(-pi) q[1];
rz(-3.0604612) q[2];
sx q[2];
rz(-0.46617026) q[2];
sx q[2];
rz(2.1612338) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7114746) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(-0.75573604) q[1];
rz(1.3487885) q[3];
sx q[3];
rz(-1.7553925) q[3];
sx q[3];
rz(0.27664646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-2.5640326) q[2];
rz(-1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-0.66453385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98786551) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(-2.7541449) q[0];
rz(-0.93915835) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.739025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57099045) q[0];
sx q[0];
rz(-1.6002858) q[0];
sx q[0];
rz(1.5667314) q[0];
rz(-pi) q[1];
rz(2.8005373) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(-1.7797433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.95784159) q[1];
sx q[1];
rz(-1.9722003) q[1];
sx q[1];
rz(-0.5387696) q[1];
rz(-pi) q[2];
rz(-1.6658695) q[3];
sx q[3];
rz(-2.2007696) q[3];
sx q[3];
rz(2.7271987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7188321) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(0.20673949) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(-2.3247705) q[3];
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
rz(1.7725672) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(-1.4136219) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51940489) q[0];
sx q[0];
rz(-1.4675958) q[0];
sx q[0];
rz(-1.2445356) q[0];
rz(-pi) q[1];
rz(-2.8917519) q[2];
sx q[2];
rz(-0.65953883) q[2];
sx q[2];
rz(1.2741054) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.229278) q[1];
sx q[1];
rz(-2.4114128) q[1];
sx q[1];
rz(0.0095403949) q[1];
rz(-pi) q[2];
rz(2.9455455) q[3];
sx q[3];
rz(-1.184433) q[3];
sx q[3];
rz(2.8783609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5473189) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(0.55580124) q[2];
rz(-0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926086) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(-1.697631) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(-2.8881883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69617535) q[0];
sx q[0];
rz(-1.6735958) q[0];
sx q[0];
rz(-1.0284543) q[0];
x q[1];
rz(-0.01732145) q[2];
sx q[2];
rz(-2.0565196) q[2];
sx q[2];
rz(2.3522365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.464444) q[1];
sx q[1];
rz(-0.943999) q[1];
sx q[1];
rz(-1.4518751) q[1];
rz(-2.4921791) q[3];
sx q[3];
rz(-1.6522437) q[3];
sx q[3];
rz(-2.7922975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(3.0786247) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(2.6834992) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8040745) q[0];
sx q[0];
rz(-2.6989557) q[0];
sx q[0];
rz(3.1361561) q[0];
rz(2.1796218) q[2];
sx q[2];
rz(-1.95032) q[2];
sx q[2];
rz(-2.6517207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.749436) q[1];
sx q[1];
rz(-1.6877618) q[1];
sx q[1];
rz(-0.82526308) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0714508) q[3];
sx q[3];
rz(-1.6541012) q[3];
sx q[3];
rz(1.7617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(-0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5620419) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(3.0793076) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.8744291) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29612449) q[0];
sx q[0];
rz(-0.86588973) q[0];
sx q[0];
rz(1.0494997) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3976372) q[2];
sx q[2];
rz(-1.6774872) q[2];
sx q[2];
rz(2.8731186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8333203) q[1];
sx q[1];
rz(-2.7566524) q[1];
sx q[1];
rz(-1.0233378) q[1];
rz(-2.0503067) q[3];
sx q[3];
rz(-1.3538085) q[3];
sx q[3];
rz(-0.12245164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(0.70872712) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-1.0725347) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(-2.0297208) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5522963) q[0];
sx q[0];
rz(-2.6483905) q[0];
sx q[0];
rz(1.769355) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6107437) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(1.2209148) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0791196) q[1];
sx q[1];
rz(-1.2812496) q[1];
sx q[1];
rz(1.4909966) q[1];
rz(-pi) q[2];
rz(0.38970077) q[3];
sx q[3];
rz(-2.0412363) q[3];
sx q[3];
rz(-2.5999992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.973935) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325571) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-2.2413975) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-2.3866167) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46144339) q[0];
sx q[0];
rz(-1.4247243) q[0];
sx q[0];
rz(2.2938674) q[0];
x q[1];
rz(2.1515498) q[2];
sx q[2];
rz(-2.2120737) q[2];
sx q[2];
rz(-2.5255447) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4432115) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(2.5618782) q[1];
rz(-0.98324361) q[3];
sx q[3];
rz(-1.1508815) q[3];
sx q[3];
rz(2.0411232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(-1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144142) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(-0.019502217) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7050539) q[0];
sx q[0];
rz(-2.9199613) q[0];
sx q[0];
rz(1.1987232) q[0];
rz(-pi) q[1];
rz(-2.9901864) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(0.99934794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94432482) q[1];
sx q[1];
rz(-1.7515344) q[1];
sx q[1];
rz(2.1987869) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2583624) q[3];
sx q[3];
rz(-1.7282681) q[3];
sx q[3];
rz(-1.0851932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.6000115) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-0.48450255) q[0];
rz(1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(1.9932995) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023037993) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(-1.5492357) q[0];
x q[1];
rz(2.8414367) q[2];
sx q[2];
rz(-0.65320063) q[2];
sx q[2];
rz(2.4094827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9705829) q[1];
sx q[1];
rz(-0.90921558) q[1];
sx q[1];
rz(2.9243102) q[1];
rz(-0.34928068) q[3];
sx q[3];
rz(-1.5216856) q[3];
sx q[3];
rz(1.9625361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(-3.0129516) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0317595) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-2.1784492) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(0.66394304) q[2];
sx q[2];
rz(-1.8891816) q[2];
sx q[2];
rz(1.328697) q[2];
rz(2.1309489) q[3];
sx q[3];
rz(-1.602136) q[3];
sx q[3];
rz(-1.7020561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];