OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.4829798) q[0];
sx q[0];
rz(3.5456181) q[0];
sx q[0];
rz(9.0344949) q[0];
rz(3.1080988) q[1];
sx q[1];
rz(-0.53116763) q[1];
sx q[1];
rz(2.9386428) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34863198) q[0];
sx q[0];
rz(-1.9237776) q[0];
sx q[0];
rz(0.99162942) q[0];
x q[1];
rz(-1.5076981) q[2];
sx q[2];
rz(-0.089091688) q[2];
sx q[2];
rz(2.2697946) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6773078) q[1];
sx q[1];
rz(-0.57972902) q[1];
sx q[1];
rz(2.7955758) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1139939) q[3];
sx q[3];
rz(-0.8445794) q[3];
sx q[3];
rz(-0.82311326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94886327) q[2];
sx q[2];
rz(-2.296083) q[2];
sx q[2];
rz(1.3884937) q[2];
rz(3.0697611) q[3];
sx q[3];
rz(-2.6303232) q[3];
sx q[3];
rz(-0.82740074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097505957) q[0];
sx q[0];
rz(-0.16479099) q[0];
sx q[0];
rz(2.6440115) q[0];
rz(-1.3961821) q[1];
sx q[1];
rz(-2.1131682) q[1];
sx q[1];
rz(2.5713249) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3511757) q[0];
sx q[0];
rz(-2.6714258) q[0];
sx q[0];
rz(-0.29542653) q[0];
rz(-0.79945081) q[2];
sx q[2];
rz(-0.28944566) q[2];
sx q[2];
rz(-1.4674526) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5757338) q[1];
sx q[1];
rz(-1.5354146) q[1];
sx q[1];
rz(2.4206764) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2482613) q[3];
sx q[3];
rz(-2.3260197) q[3];
sx q[3];
rz(1.4295242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1121062) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(2.4278329) q[2];
rz(-1.7589689) q[3];
sx q[3];
rz(-0.52240038) q[3];
sx q[3];
rz(-2.6944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5722028) q[0];
sx q[0];
rz(-2.1694006) q[0];
sx q[0];
rz(-0.33706459) q[0];
rz(0.49346787) q[1];
sx q[1];
rz(-2.4512873) q[1];
sx q[1];
rz(0.92672551) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0686058) q[0];
sx q[0];
rz(-1.301501) q[0];
sx q[0];
rz(-2.8728743) q[0];
rz(-0.91442911) q[2];
sx q[2];
rz(-0.85322748) q[2];
sx q[2];
rz(0.85384439) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34636075) q[1];
sx q[1];
rz(-1.7344001) q[1];
sx q[1];
rz(-1.9925881) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95229228) q[3];
sx q[3];
rz(-2.3626948) q[3];
sx q[3];
rz(-2.2924045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9530764) q[2];
sx q[2];
rz(-2.2761554) q[2];
sx q[2];
rz(-0.68149978) q[2];
rz(1.513688) q[3];
sx q[3];
rz(-1.3368006) q[3];
sx q[3];
rz(-0.17076913) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3156768) q[0];
sx q[0];
rz(-3.0382394) q[0];
sx q[0];
rz(2.090825) q[0];
rz(0.61327618) q[1];
sx q[1];
rz(-0.79137099) q[1];
sx q[1];
rz(1.9176066) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036788615) q[0];
sx q[0];
rz(-2.4147644) q[0];
sx q[0];
rz(2.6340371) q[0];
rz(-0.54398016) q[2];
sx q[2];
rz(-2.1270985) q[2];
sx q[2];
rz(-2.408319) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5889086) q[1];
sx q[1];
rz(-0.97637227) q[1];
sx q[1];
rz(-1.5556704) q[1];
x q[2];
rz(-0.59215109) q[3];
sx q[3];
rz(-1.4363406) q[3];
sx q[3];
rz(2.8176702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0918538) q[2];
sx q[2];
rz(-0.78511304) q[2];
sx q[2];
rz(-2.4741057) q[2];
rz(1.5664172) q[3];
sx q[3];
rz(-2.9508041) q[3];
sx q[3];
rz(-0.13458399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62395537) q[0];
sx q[0];
rz(-1.0409545) q[0];
sx q[0];
rz(2.5868296) q[0];
rz(1.293921) q[1];
sx q[1];
rz(-0.41176739) q[1];
sx q[1];
rz(-2.256934) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22065565) q[0];
sx q[0];
rz(-1.0405128) q[0];
sx q[0];
rz(-0.028593731) q[0];
x q[1];
rz(-1.1049461) q[2];
sx q[2];
rz(-2.2171582) q[2];
sx q[2];
rz(2.4798648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5371478) q[1];
sx q[1];
rz(-1.7824518) q[1];
sx q[1];
rz(2.0404633) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8507666) q[3];
sx q[3];
rz(-1.9970702) q[3];
sx q[3];
rz(-2.4170162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1137997) q[2];
sx q[2];
rz(-2.5774559) q[2];
sx q[2];
rz(2.0417058) q[2];
rz(-2.9187628) q[3];
sx q[3];
rz(-1.9375216) q[3];
sx q[3];
rz(0.13404624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564388) q[0];
sx q[0];
rz(-2.8414861) q[0];
sx q[0];
rz(-1.1837748) q[0];
rz(1.3166332) q[1];
sx q[1];
rz(-2.1778409) q[1];
sx q[1];
rz(-2.1060941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1610704) q[0];
sx q[0];
rz(-2.5455406) q[0];
sx q[0];
rz(-0.97575326) q[0];
rz(-pi) q[1];
rz(-0.99656099) q[2];
sx q[2];
rz(-2.445526) q[2];
sx q[2];
rz(2.0526759) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4343513) q[1];
sx q[1];
rz(-1.3907258) q[1];
sx q[1];
rz(1.7456013) q[1];
rz(1.6430969) q[3];
sx q[3];
rz(-2.3359155) q[3];
sx q[3];
rz(1.3209381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3731132) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(-0.38144544) q[2];
rz(-1.2162195) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(-0.60429627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(2.4878047) q[0];
sx q[0];
rz(-0.31823802) q[0];
sx q[0];
rz(-2.8741264) q[0];
rz(-1.5221315) q[1];
sx q[1];
rz(-0.61264241) q[1];
sx q[1];
rz(-2.6681275) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2758613) q[0];
sx q[0];
rz(-1.9124219) q[0];
sx q[0];
rz(-2.7462358) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47961819) q[2];
sx q[2];
rz(-1.4066182) q[2];
sx q[2];
rz(-0.93725433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50187868) q[1];
sx q[1];
rz(-0.25020975) q[1];
sx q[1];
rz(-2.6544957) q[1];
rz(-1.7340937) q[3];
sx q[3];
rz(-1.9400404) q[3];
sx q[3];
rz(2.6881517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9240616) q[2];
sx q[2];
rz(-1.6465829) q[2];
sx q[2];
rz(1.7688497) q[2];
rz(0.30630201) q[3];
sx q[3];
rz(-1.0270216) q[3];
sx q[3];
rz(0.33941227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.821625) q[0];
sx q[0];
rz(-0.29682934) q[0];
sx q[0];
rz(2.6676275) q[0];
rz(-2.3560246) q[1];
sx q[1];
rz(-0.57890099) q[1];
sx q[1];
rz(-0.89326352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2300917) q[0];
sx q[0];
rz(-2.8727838) q[0];
sx q[0];
rz(-1.5066766) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39126663) q[2];
sx q[2];
rz(-2.3377468) q[2];
sx q[2];
rz(0.93144691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6345659) q[1];
sx q[1];
rz(-2.9302672) q[1];
sx q[1];
rz(2.9853249) q[1];
x q[2];
rz(1.8056554) q[3];
sx q[3];
rz(-2.7357172) q[3];
sx q[3];
rz(1.0535976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6713509) q[2];
sx q[2];
rz(-1.525815) q[2];
sx q[2];
rz(-0.5980171) q[2];
rz(2.9371069) q[3];
sx q[3];
rz(-3.0308767) q[3];
sx q[3];
rz(2.428875) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0996899) q[0];
sx q[0];
rz(-0.99869096) q[0];
sx q[0];
rz(2.4424851) q[0];
rz(2.7514669) q[1];
sx q[1];
rz(-2.4603619) q[1];
sx q[1];
rz(2.1652538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49347827) q[0];
sx q[0];
rz(-2.8472487) q[0];
sx q[0];
rz(-1.0674547) q[0];
x q[1];
rz(1.8985229) q[2];
sx q[2];
rz(-1.9877745) q[2];
sx q[2];
rz(2.8627739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55667952) q[1];
sx q[1];
rz(-2.0784857) q[1];
sx q[1];
rz(-0.084225144) q[1];
rz(-pi) q[2];
rz(-0.41918079) q[3];
sx q[3];
rz(-0.59801918) q[3];
sx q[3];
rz(2.647915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1583027) q[2];
sx q[2];
rz(-2.172564) q[2];
sx q[2];
rz(0.14652531) q[2];
rz(2.8807785) q[3];
sx q[3];
rz(-1.1365889) q[3];
sx q[3];
rz(-0.28723106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2448267) q[0];
sx q[0];
rz(-0.34938669) q[0];
sx q[0];
rz(-0.31325999) q[0];
rz(-0.80097711) q[1];
sx q[1];
rz(-1.5104834) q[1];
sx q[1];
rz(2.7371178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069720751) q[0];
sx q[0];
rz(-1.3861462) q[0];
sx q[0];
rz(1.429256) q[0];
rz(1.7687665) q[2];
sx q[2];
rz(-2.1083197) q[2];
sx q[2];
rz(-0.24453577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69533379) q[1];
sx q[1];
rz(-2.0768169) q[1];
sx q[1];
rz(-2.395588) q[1];
x q[2];
rz(3.0048921) q[3];
sx q[3];
rz(-1.2368844) q[3];
sx q[3];
rz(-2.2025618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6910088) q[2];
sx q[2];
rz(-2.6477224) q[2];
sx q[2];
rz(0.79130006) q[2];
rz(2.6386236) q[3];
sx q[3];
rz(-1.0478323) q[3];
sx q[3];
rz(-0.80016971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7293411) q[0];
sx q[0];
rz(-1.3164192) q[0];
sx q[0];
rz(-1.3715716) q[0];
rz(0.62660632) q[1];
sx q[1];
rz(-1.659844) q[1];
sx q[1];
rz(-1.0214092) q[1];
rz(1.7328429) q[2];
sx q[2];
rz(-1.2037983) q[2];
sx q[2];
rz(-1.498602) q[2];
rz(-3.088763) q[3];
sx q[3];
rz(-0.23374791) q[3];
sx q[3];
rz(1.4240627) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
