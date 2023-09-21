OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(3.7319558) q[0];
sx q[0];
rz(9.0537602) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(-1.7655656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57173079) q[0];
sx q[0];
rz(-0.96643448) q[0];
sx q[0];
rz(-0.5423003) q[0];
rz(-pi) q[1];
rz(0.55144989) q[2];
sx q[2];
rz(-2.3460238) q[2];
sx q[2];
rz(-1.431682) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10378097) q[1];
sx q[1];
rz(-1.4191287) q[1];
sx q[1];
rz(1.2396461) q[1];
rz(-pi) q[2];
rz(1.6269496) q[3];
sx q[3];
rz(-2.3893642) q[3];
sx q[3];
rz(-0.55263954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.084289) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(-0.75254285) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.9616615) q[0];
rz(0.99769366) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-0.72431272) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99740072) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(2.2749167) q[0];
rz(-2.9421147) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(2.2897838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1867379) q[1];
sx q[1];
rz(-1.9222944) q[1];
sx q[1];
rz(0.21912205) q[1];
rz(-pi) q[2];
rz(-0.23726666) q[3];
sx q[3];
rz(-1.7627343) q[3];
sx q[3];
rz(0.81258472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90536845) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(2.779707) q[2];
rz(-3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(-1.3954337) q[0];
rz(-0.46229258) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.2190855) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1299382) q[0];
sx q[0];
rz(-1.219795) q[0];
sx q[0];
rz(-2.8610693) q[0];
rz(-pi) q[1];
rz(-2.2723324) q[2];
sx q[2];
rz(-1.7995036) q[2];
sx q[2];
rz(-1.9321835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29875007) q[1];
sx q[1];
rz(-1.0423653) q[1];
sx q[1];
rz(-2.8612086) q[1];
rz(-pi) q[2];
rz(-1.5145281) q[3];
sx q[3];
rz(-2.8731822) q[3];
sx q[3];
rz(-1.0750107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(-3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179984) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(-0.82675654) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0543538) q[0];
sx q[0];
rz(-1.0687437) q[0];
sx q[0];
rz(-2.0421844) q[0];
rz(2.3739359) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(-1.8061639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3863694) q[1];
sx q[1];
rz(-1.1896903) q[1];
sx q[1];
rz(1.8783046) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8825674) q[3];
sx q[3];
rz(-0.87083737) q[3];
sx q[3];
rz(1.242897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.15239079) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82692659) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(0.28516969) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7171705) q[0];
sx q[0];
rz(-1.5283913) q[0];
sx q[0];
rz(2.0457343) q[0];
x q[1];
rz(-1.9589013) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(2.2270122) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3030745) q[1];
sx q[1];
rz(-0.31305602) q[1];
sx q[1];
rz(-1.5721553) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4739591) q[3];
sx q[3];
rz(-1.1708461) q[3];
sx q[3];
rz(-1.6638343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(-2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.3670115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0317504) q[0];
sx q[0];
rz(-1.5378008) q[0];
sx q[0];
rz(-1.5026291) q[0];
x q[1];
rz(-0.579367) q[2];
sx q[2];
rz(-2.2420792) q[2];
sx q[2];
rz(-2.6210149) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64703343) q[1];
sx q[1];
rz(-1.4899947) q[1];
sx q[1];
rz(-1.8606436) q[1];
x q[2];
rz(-0.35297024) q[3];
sx q[3];
rz(-1.6680696) q[3];
sx q[3];
rz(-2.4091165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-2.7381251) q[2];
rz(-0.48163313) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.1706932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6252977) q[0];
sx q[0];
rz(-1.5646311) q[0];
sx q[0];
rz(0.023029285) q[0];
x q[1];
rz(-0.41202338) q[2];
sx q[2];
rz(-0.12672666) q[2];
sx q[2];
rz(2.3840981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2611744) q[1];
sx q[1];
rz(-2.2194127) q[1];
sx q[1];
rz(0.76678126) q[1];
x q[2];
rz(0.4865173) q[3];
sx q[3];
rz(-2.5518637) q[3];
sx q[3];
rz(0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0968904) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24191813) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(0.29712594) q[0];
rz(-1.3946474) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(-0.64613211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0028249) q[0];
sx q[0];
rz(-1.7768304) q[0];
sx q[0];
rz(-0.10319184) q[0];
rz(-pi) q[1];
rz(1.4023151) q[2];
sx q[2];
rz(-2.0373166) q[2];
sx q[2];
rz(2.434935) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7039973) q[1];
sx q[1];
rz(-2.5181209) q[1];
sx q[1];
rz(2.5295579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5313247) q[3];
sx q[3];
rz(-2.2329997) q[3];
sx q[3];
rz(2.3074647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.064676553) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(0.45483744) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(-2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(-3.033175) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71589564) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(-0.9066559) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8234532) q[2];
sx q[2];
rz(-0.18441072) q[2];
sx q[2];
rz(2.2968963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3474116) q[1];
sx q[1];
rz(-1.3550183) q[1];
sx q[1];
rz(-2.9706035) q[1];
rz(-pi) q[2];
rz(-2.9547144) q[3];
sx q[3];
rz(-2.2664321) q[3];
sx q[3];
rz(-2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0124399) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6091992) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(-0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(-3.0864339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.178135) q[0];
sx q[0];
rz(-1.9473416) q[0];
sx q[0];
rz(-1.5012653) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.360838) q[2];
sx q[2];
rz(-2.1205489) q[2];
sx q[2];
rz(1.9090261) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.060505796) q[1];
sx q[1];
rz(-0.85716893) q[1];
sx q[1];
rz(0.75018261) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0268029) q[3];
sx q[3];
rz(-1.6037233) q[3];
sx q[3];
rz(-1.1820716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(-0.49017635) q[2];
rz(-3.0040719) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162083) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(2.8293777) q[2];
sx q[2];
rz(-1.6630465) q[2];
sx q[2];
rz(1.9065471) q[2];
rz(2.901554) q[3];
sx q[3];
rz(-1.6143027) q[3];
sx q[3];
rz(-1.3054813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
