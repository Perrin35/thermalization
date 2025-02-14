OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1825778) q[0];
sx q[0];
rz(-0.20354095) q[0];
sx q[0];
rz(1.8939053) q[0];
rz(-1.7893451) q[1];
sx q[1];
rz(-2.4440553) q[1];
sx q[1];
rz(0.69812671) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20839918) q[0];
sx q[0];
rz(-1.0675479) q[0];
sx q[0];
rz(1.1017975) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1280367) q[2];
sx q[2];
rz(-0.45326158) q[2];
sx q[2];
rz(1.0176941) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0586529) q[1];
sx q[1];
rz(-1.8720683) q[1];
sx q[1];
rz(2.3430587) q[1];
rz(-pi) q[2];
rz(1.5145958) q[3];
sx q[3];
rz(-0.25543419) q[3];
sx q[3];
rz(-1.5843624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2932566) q[2];
sx q[2];
rz(-1.2245327) q[2];
sx q[2];
rz(-2.2621034) q[2];
rz(0.73421156) q[3];
sx q[3];
rz(-1.1314393) q[3];
sx q[3];
rz(0.82936275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64049983) q[0];
sx q[0];
rz(-1.1859897) q[0];
sx q[0];
rz(-0.99047852) q[0];
rz(-2.1462006) q[1];
sx q[1];
rz(-1.4442911) q[1];
sx q[1];
rz(2.676414) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93056193) q[0];
sx q[0];
rz(-0.4586691) q[0];
sx q[0];
rz(-0.42210292) q[0];
rz(-2.0700109) q[2];
sx q[2];
rz(-2.4518161) q[2];
sx q[2];
rz(-0.76329714) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8930186) q[1];
sx q[1];
rz(-2.0573528) q[1];
sx q[1];
rz(-2.3793329) q[1];
rz(1.7725132) q[3];
sx q[3];
rz(-1.3139825) q[3];
sx q[3];
rz(-1.4434003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5140932) q[2];
sx q[2];
rz(-1.2343312) q[2];
sx q[2];
rz(0.17847432) q[2];
rz(-2.576339) q[3];
sx q[3];
rz(-1.7491128) q[3];
sx q[3];
rz(0.55725151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8239215) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(2.4406261) q[0];
rz(2.7340381) q[1];
sx q[1];
rz(-1.9302255) q[1];
sx q[1];
rz(2.0661381) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70301818) q[0];
sx q[0];
rz(-2.10963) q[0];
sx q[0];
rz(1.1807022) q[0];
x q[1];
rz(-2.2268217) q[2];
sx q[2];
rz(-1.390995) q[2];
sx q[2];
rz(0.54043661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4891995) q[1];
sx q[1];
rz(-1.7013738) q[1];
sx q[1];
rz(1.9741535) q[1];
rz(-pi) q[2];
x q[2];
rz(1.386679) q[3];
sx q[3];
rz(-0.26911165) q[3];
sx q[3];
rz(-1.5003913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0740697) q[2];
sx q[2];
rz(-1.8068376) q[2];
sx q[2];
rz(1.734181) q[2];
rz(-0.47809005) q[3];
sx q[3];
rz(-2.7985088) q[3];
sx q[3];
rz(-0.34496719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5673229) q[0];
sx q[0];
rz(-1.3936717) q[0];
sx q[0];
rz(-2.0950914) q[0];
rz(-2.8872755) q[1];
sx q[1];
rz(-2.3260702) q[1];
sx q[1];
rz(-2.8291124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1069665) q[0];
sx q[0];
rz(-2.7101332) q[0];
sx q[0];
rz(1.8851243) q[0];
rz(2.0489659) q[2];
sx q[2];
rz(-2.1743283) q[2];
sx q[2];
rz(2.9750533) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29280832) q[1];
sx q[1];
rz(-2.7393746) q[1];
sx q[1];
rz(2.4821698) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.11361) q[3];
sx q[3];
rz(-0.58251721) q[3];
sx q[3];
rz(-1.5912671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8326524) q[2];
sx q[2];
rz(-0.71844429) q[2];
sx q[2];
rz(-0.18860513) q[2];
rz(-2.9077933) q[3];
sx q[3];
rz(-2.2457687) q[3];
sx q[3];
rz(2.5578267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.69498649) q[0];
sx q[0];
rz(-1.3124895) q[0];
sx q[0];
rz(-3.1332698) q[0];
rz(0.43909973) q[1];
sx q[1];
rz(-1.7726026) q[1];
sx q[1];
rz(1.2459374) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74730342) q[0];
sx q[0];
rz(-0.4957333) q[0];
sx q[0];
rz(0.41252211) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8077085) q[2];
sx q[2];
rz(-1.4258175) q[2];
sx q[2];
rz(-2.0857834) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.627927) q[1];
sx q[1];
rz(-0.49066077) q[1];
sx q[1];
rz(-2.2918244) q[1];
rz(-pi) q[2];
rz(1.2034791) q[3];
sx q[3];
rz(-2.4532336) q[3];
sx q[3];
rz(1.7662774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40450725) q[2];
sx q[2];
rz(-0.046701996) q[2];
sx q[2];
rz(1.6339634) q[2];
rz(0.59219939) q[3];
sx q[3];
rz(-1.3785572) q[3];
sx q[3];
rz(-0.73970214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020788766) q[0];
sx q[0];
rz(-1.692166) q[0];
sx q[0];
rz(-2.8756496) q[0];
rz(1.9256915) q[1];
sx q[1];
rz(-2.3561056) q[1];
sx q[1];
rz(-1.9507834) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94534369) q[0];
sx q[0];
rz(-0.73305819) q[0];
sx q[0];
rz(-1.5668012) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1235775) q[2];
sx q[2];
rz(-1.6016804) q[2];
sx q[2];
rz(-0.72252233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9126351) q[1];
sx q[1];
rz(-1.5785145) q[1];
sx q[1];
rz(-2.387729) q[1];
x q[2];
rz(-0.27482185) q[3];
sx q[3];
rz(-1.6259818) q[3];
sx q[3];
rz(1.7589057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8094981) q[2];
sx q[2];
rz(-1.3846493) q[2];
sx q[2];
rz(2.9936227) q[2];
rz(0.45251265) q[3];
sx q[3];
rz(-1.0871525) q[3];
sx q[3];
rz(-2.7361659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0026523503) q[0];
sx q[0];
rz(-1.8754706) q[0];
sx q[0];
rz(-1.3462322) q[0];
rz(-3.1345308) q[1];
sx q[1];
rz(-0.66771737) q[1];
sx q[1];
rz(2.6341569) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1088828) q[0];
sx q[0];
rz(-1.0531972) q[0];
sx q[0];
rz(-0.93982307) q[0];
rz(-2.9865737) q[2];
sx q[2];
rz(-1.9117725) q[2];
sx q[2];
rz(2.2237958) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40312156) q[1];
sx q[1];
rz(-2.7475121) q[1];
sx q[1];
rz(3.1172063) q[1];
x q[2];
rz(1.5370338) q[3];
sx q[3];
rz(-0.27245263) q[3];
sx q[3];
rz(-1.1794943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4355882) q[2];
sx q[2];
rz(-1.0554577) q[2];
sx q[2];
rz(-0.17534193) q[2];
rz(-0.38241479) q[3];
sx q[3];
rz(-1.9491111) q[3];
sx q[3];
rz(-2.6800938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14313702) q[0];
sx q[0];
rz(-2.1303506) q[0];
sx q[0];
rz(1.829041) q[0];
rz(0.10874272) q[1];
sx q[1];
rz(-0.60832945) q[1];
sx q[1];
rz(2.463602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10763845) q[0];
sx q[0];
rz(-2.0642515) q[0];
sx q[0];
rz(-3.1163868) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74900519) q[2];
sx q[2];
rz(-1.5856447) q[2];
sx q[2];
rz(-0.9204364) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99876632) q[1];
sx q[1];
rz(-1.8590392) q[1];
sx q[1];
rz(-2.114813) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59321065) q[3];
sx q[3];
rz(-2.3311391) q[3];
sx q[3];
rz(2.5814501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8678191) q[2];
sx q[2];
rz(-1.5673693) q[2];
sx q[2];
rz(-1.1567953) q[2];
rz(0.50446883) q[3];
sx q[3];
rz(-1.0320832) q[3];
sx q[3];
rz(-0.57197905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93541637) q[0];
sx q[0];
rz(-1.6125866) q[0];
sx q[0];
rz(-2.5231498) q[0];
rz(-0.84856021) q[1];
sx q[1];
rz(-2.6237374) q[1];
sx q[1];
rz(-2.3458164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9002118) q[0];
sx q[0];
rz(-2.6707895) q[0];
sx q[0];
rz(1.5174675) q[0];
rz(0.82473849) q[2];
sx q[2];
rz(-1.6055664) q[2];
sx q[2];
rz(0.29037133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.578131) q[1];
sx q[1];
rz(-1.8860779) q[1];
sx q[1];
rz(0.28750802) q[1];
rz(0.14988092) q[3];
sx q[3];
rz(-1.9871885) q[3];
sx q[3];
rz(-0.021377953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.15647469) q[2];
sx q[2];
rz(-1.5072344) q[2];
sx q[2];
rz(-1.1953243) q[2];
rz(-0.34504238) q[3];
sx q[3];
rz(-2.2796567) q[3];
sx q[3];
rz(-0.8663469) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4320375) q[0];
sx q[0];
rz(-2.9032752) q[0];
sx q[0];
rz(0.077614345) q[0];
rz(-1.2331102) q[1];
sx q[1];
rz(-2.1014919) q[1];
sx q[1];
rz(-0.79201039) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38538489) q[0];
sx q[0];
rz(-2.1127708) q[0];
sx q[0];
rz(0.78800027) q[0];
rz(-pi) q[1];
rz(1.4155238) q[2];
sx q[2];
rz(-1.4129935) q[2];
sx q[2];
rz(-1.4580245) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35490552) q[1];
sx q[1];
rz(-0.49017683) q[1];
sx q[1];
rz(-1.2398943) q[1];
rz(-0.11802063) q[3];
sx q[3];
rz(-2.3343918) q[3];
sx q[3];
rz(-2.0995875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9218257) q[2];
sx q[2];
rz(-2.2302088) q[2];
sx q[2];
rz(1.0731953) q[2];
rz(-2.8940767) q[3];
sx q[3];
rz(-1.211834) q[3];
sx q[3];
rz(3.1246576) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76219227) q[0];
sx q[0];
rz(-1.301855) q[0];
sx q[0];
rz(-0.72574885) q[0];
rz(-0.87860592) q[1];
sx q[1];
rz(-0.85522063) q[1];
sx q[1];
rz(-1.9488889) q[1];
rz(-2.939777) q[2];
sx q[2];
rz(-2.0463662) q[2];
sx q[2];
rz(-3.0007888) q[2];
rz(0.056817354) q[3];
sx q[3];
rz(-2.29689) q[3];
sx q[3];
rz(2.5162027) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
