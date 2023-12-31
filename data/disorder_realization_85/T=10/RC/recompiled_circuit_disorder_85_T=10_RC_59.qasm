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
rz(-2.5512295) q[0];
sx q[0];
rz(-0.37101775) q[0];
rz(2.7603005) q[1];
sx q[1];
rz(-2.5420904) q[1];
sx q[1];
rz(-1.376027) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3289514) q[0];
sx q[0];
rz(-1.1321804) q[0];
sx q[0];
rz(-2.2493275) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71557578) q[2];
sx q[2];
rz(-1.1872429) q[2];
sx q[2];
rz(-2.5959612) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4151167) q[1];
sx q[1];
rz(-1.8980025) q[1];
sx q[1];
rz(2.9813558) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81935482) q[3];
sx q[3];
rz(-1.5324394) q[3];
sx q[3];
rz(1.0591782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(2.1526745) q[2];
rz(-0.75254285) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(2.4172799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82942078) q[0];
sx q[0];
rz(-2.2342626) q[0];
sx q[0];
rz(-2.7396766) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7472277) q[2];
sx q[2];
rz(-0.21557237) q[2];
sx q[2];
rz(2.8087316) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7608632) q[1];
sx q[1];
rz(-0.41178307) q[1];
sx q[1];
rz(2.1058583) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69085391) q[3];
sx q[3];
rz(-2.8375531) q[3];
sx q[3];
rz(-3.0512878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(-0.36188564) q[2];
rz(-0.13606717) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.3954337) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.2190855) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31419471) q[0];
sx q[0];
rz(-2.6959246) q[0];
sx q[0];
rz(-0.92339869) q[0];
x q[1];
rz(2.8457882) q[2];
sx q[2];
rz(-2.25053) q[2];
sx q[2];
rz(-2.5909397) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8176986) q[1];
sx q[1];
rz(-2.5496799) q[1];
sx q[1];
rz(2.0134258) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3027906) q[3];
sx q[3];
rz(-1.5857113) q[3];
sx q[3];
rz(-0.55004317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(-2.5727663) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(-3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-2.9179984) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(2.3148361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2430192) q[0];
sx q[0];
rz(-1.9802226) q[0];
sx q[0];
rz(0.55222521) q[0];
rz(-0.13917285) q[2];
sx q[2];
rz(-2.3690802) q[2];
sx q[2];
rz(-0.13538361) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.066847853) q[1];
sx q[1];
rz(-1.8555992) q[1];
sx q[1];
rz(0.39798255) q[1];
x q[2];
rz(1.8825674) q[3];
sx q[3];
rz(-2.2707553) q[3];
sx q[3];
rz(1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9892019) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(-1.4962083) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3146661) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(2.856423) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171705) q[0];
sx q[0];
rz(-1.5283913) q[0];
sx q[0];
rz(1.0958584) q[0];
rz(-pi) q[1];
rz(-1.1826913) q[2];
sx q[2];
rz(-2.1974807) q[2];
sx q[2];
rz(2.2270122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83708977) q[1];
sx q[1];
rz(-1.883852) q[1];
sx q[1];
rz(-0.00043991107) q[1];
x q[2];
rz(0.40163715) q[3];
sx q[3];
rz(-1.65997) q[3];
sx q[3];
rz(-0.13084403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-2.6021393) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-2.6873798) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(-2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(1.7745811) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088928662) q[0];
sx q[0];
rz(-3.0658709) q[0];
sx q[0];
rz(-1.1195539) q[0];
x q[1];
rz(-0.579367) q[2];
sx q[2];
rz(-2.2420792) q[2];
sx q[2];
rz(-2.6210149) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4821266) q[1];
sx q[1];
rz(-2.8409993) q[1];
sx q[1];
rz(1.2947004) q[1];
rz(2.7886224) q[3];
sx q[3];
rz(-1.4735231) q[3];
sx q[3];
rz(2.4091165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8453025) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(0.48163313) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24213174) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(-1.1706932) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51629492) q[0];
sx q[0];
rz(-1.5769616) q[0];
sx q[0];
rz(3.1185634) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0253719) q[2];
sx q[2];
rz(-1.6214317) q[2];
sx q[2];
rz(-1.2223787) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.891174) q[1];
sx q[1];
rz(-0.95953566) q[1];
sx q[1];
rz(0.82959081) q[1];
rz(-pi) q[2];
rz(0.53415926) q[3];
sx q[3];
rz(-1.8338406) q[3];
sx q[3];
rz(-1.1155038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0447023) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(2.5308385) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8996745) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-0.64613211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32983366) q[0];
sx q[0];
rz(-0.23010294) q[0];
sx q[0];
rz(2.0287201) q[0];
rz(0.32142873) q[2];
sx q[2];
rz(-2.6476963) q[2];
sx q[2];
rz(2.0733881) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99107332) q[1];
sx q[1];
rz(-1.0725613) q[1];
sx q[1];
rz(-1.9626161) q[1];
rz(-pi) q[2];
rz(0.81031462) q[3];
sx q[3];
rz(-2.0397566) q[3];
sx q[3];
rz(-1.9988434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0769161) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(0.45483744) q[2];
rz(0.70139766) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(-1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(0.79750693) q[0];
rz(2.6240255) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(-3.033175) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3603044) q[0];
sx q[0];
rz(-0.66626781) q[0];
sx q[0];
rz(-1.6643307) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3181395) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(-0.84469634) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3280231) q[1];
sx q[1];
rz(-1.4038101) q[1];
sx q[1];
rz(-1.7896673) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9547144) q[3];
sx q[3];
rz(-0.87516057) q[3];
sx q[3];
rz(2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0124399) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(-2.7774096) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(-0.055158786) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3653152) q[0];
sx q[0];
rz(-0.3826097) q[0];
sx q[0];
rz(2.9676653) q[0];
x q[1];
rz(-2.8137384) q[2];
sx q[2];
rz(-0.58460669) q[2];
sx q[2];
rz(2.2962928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0835417) q[1];
sx q[1];
rz(-1.0293048) q[1];
sx q[1];
rz(2.4400649) q[1];
x q[2];
rz(1.6454562) q[3];
sx q[3];
rz(-2.6844822) q[3];
sx q[3];
rz(0.32170579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(-2.6514163) q[2];
rz(0.13752078) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.4162083) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(-2.9329119) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(-1.4738884) q[2];
sx q[2];
rz(-1.8816392) q[2];
sx q[2];
rz(0.30602602) q[2];
rz(-2.901554) q[3];
sx q[3];
rz(-1.5272899) q[3];
sx q[3];
rz(1.8361113) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
