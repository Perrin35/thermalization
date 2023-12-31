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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8986172) q[0];
sx q[0];
rz(-0.78865047) q[0];
sx q[0];
rz(-2.2126161) q[0];
rz(-pi) q[1];
rz(-0.71557578) q[2];
sx q[2];
rz(-1.9543497) q[2];
sx q[2];
rz(-0.54563145) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7264759) q[1];
sx q[1];
rz(-1.2435902) q[1];
sx q[1];
rz(-2.9813558) q[1];
rz(-0.81935482) q[3];
sx q[3];
rz(-1.5324394) q[3];
sx q[3];
rz(2.0824144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0573037) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(0.99769366) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-0.72431272) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2253101) q[0];
sx q[0];
rz(-0.75964576) q[0];
sx q[0];
rz(2.0347974) q[0];
rz(-pi) q[1];
rz(-1.486859) q[2];
sx q[2];
rz(-1.3720023) q[2];
sx q[2];
rz(-2.4059911) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.38072941) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(-2.1058583) q[1];
x q[2];
rz(0.23726666) q[3];
sx q[3];
rz(-1.3788584) q[3];
sx q[3];
rz(-2.3290079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(0.13606717) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(-3.0959685) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(1.3954337) q[0];
rz(-2.6793001) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(1.9225072) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8273979) q[0];
sx q[0];
rz(-0.44566804) q[0];
sx q[0];
rz(2.218194) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9169541) q[2];
sx q[2];
rz(-0.73181728) q[2];
sx q[2];
rz(3.0423622) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-2.5496799) q[1];
sx q[1];
rz(-2.0134258) q[1];
rz(-pi) q[2];
rz(-1.6270646) q[3];
sx q[3];
rz(-0.26841044) q[3];
sx q[3];
rz(-1.0750107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.4455618) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(-3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(2.4261684) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(2.3148361) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137705) q[0];
sx q[0];
rz(-0.6745406) q[0];
sx q[0];
rz(-0.69112372) q[0];
x q[1];
rz(2.3739359) q[2];
sx q[2];
rz(-1.4738238) q[2];
sx q[2];
rz(-1.3354288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0930867) q[1];
sx q[1];
rz(-0.48492453) q[1];
sx q[1];
rz(2.4946458) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7923146) q[3];
sx q[3];
rz(-2.3861902) q[3];
sx q[3];
rz(0.7790156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.3789122) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(-2.0157053) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(-0.28516969) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42442214) q[0];
sx q[0];
rz(-1.6132014) q[0];
sx q[0];
rz(1.0958584) q[0];
rz(1.1826913) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(2.2270122) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4080216) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(1.2577406) q[1];
x q[2];
rz(1.4739591) q[3];
sx q[3];
rz(-1.9707465) q[3];
sx q[3];
rz(-1.6638343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(0.53945333) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-0.15701292) q[0];
rz(-0.69333386) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.7745811) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088928662) q[0];
sx q[0];
rz(-0.0757218) q[0];
sx q[0];
rz(2.0220387) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5622257) q[2];
sx q[2];
rz(-2.2420792) q[2];
sx q[2];
rz(-0.52057779) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94783084) q[1];
sx q[1];
rz(-1.281922) q[1];
sx q[1];
rz(0.08430251) q[1];
x q[2];
rz(2.7886224) q[3];
sx q[3];
rz(-1.4735231) q[3];
sx q[3];
rz(-0.73247611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8453025) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(-0.40346754) q[2];
rz(0.48163313) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(-2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(-0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.1706932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0543594) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(1.5646294) q[0];
rz(-3.0253719) q[2];
sx q[2];
rz(-1.6214317) q[2];
sx q[2];
rz(1.2223787) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16329855) q[1];
sx q[1];
rz(-0.98493176) q[1];
sx q[1];
rz(-2.3818124) q[1];
rz(-pi) q[2];
rz(-0.4865173) q[3];
sx q[3];
rz(-2.5518637) q[3];
sx q[3];
rz(-0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0968904) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(-0.61075413) q[2];
rz(0.47510535) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(-0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-0.29712594) q[0];
rz(1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532115) q[0];
sx q[0];
rz(-1.4697945) q[0];
sx q[0];
rz(-1.7779011) q[0];
rz(-2.8201639) q[2];
sx q[2];
rz(-0.49389631) q[2];
sx q[2];
rz(-2.0733881) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1505193) q[1];
sx q[1];
rz(-1.0725613) q[1];
sx q[1];
rz(-1.9626161) q[1];
x q[2];
rz(-0.81031462) q[3];
sx q[3];
rz(-2.0397566) q[3];
sx q[3];
rz(-1.1427493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0769161) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(0.45483744) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(-1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4421473) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(0.79750693) q[0];
rz(2.6240255) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(-3.033175) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78128821) q[0];
sx q[0];
rz(-0.66626781) q[0];
sx q[0];
rz(-1.6643307) q[0];
rz(-pi) q[1];
x q[1];
rz(0.046594521) q[2];
sx q[2];
rz(-1.7492883) q[2];
sx q[2];
rz(0.5878693) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11549599) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(2.2309169) q[1];
rz(1.3518203) q[3];
sx q[3];
rz(-0.7162381) q[3];
sx q[3];
rz(2.4718474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(2.8016395) q[2];
rz(-0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(3.0864339) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3653152) q[0];
sx q[0];
rz(-0.3826097) q[0];
sx q[0];
rz(2.9676653) q[0];
rz(-pi) q[1];
rz(0.32785428) q[2];
sx q[2];
rz(-2.556986) q[2];
sx q[2];
rz(0.84529982) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0810869) q[1];
sx q[1];
rz(-0.85716893) q[1];
sx q[1];
rz(-0.75018261) q[1];
x q[2];
rz(1.4961365) q[3];
sx q[3];
rz(-2.6844822) q[3];
sx q[3];
rz(2.8198869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(-0.49017635) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162083) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(1.6677042) q[2];
sx q[2];
rz(-1.8816392) q[2];
sx q[2];
rz(0.30602602) q[2];
rz(-0.18110885) q[3];
sx q[3];
rz(-2.8977179) q[3];
sx q[3];
rz(0.44117622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
