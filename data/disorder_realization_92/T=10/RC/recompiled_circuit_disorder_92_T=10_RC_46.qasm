OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5151514) q[0];
sx q[0];
rz(-0.03349537) q[0];
sx q[0];
rz(1.7749696) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(4.7935901) q[1];
sx q[1];
rz(10.556769) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6080947) q[0];
sx q[0];
rz(-0.84871263) q[0];
sx q[0];
rz(0.46790926) q[0];
x q[1];
rz(-1.7366473) q[2];
sx q[2];
rz(-1.7300786) q[2];
sx q[2];
rz(0.35370358) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1177897) q[1];
sx q[1];
rz(-2.1315247) q[1];
sx q[1];
rz(2.6610713) q[1];
x q[2];
rz(-1.8559998) q[3];
sx q[3];
rz(-1.0343026) q[3];
sx q[3];
rz(-2.5453018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.91036096) q[2];
sx q[2];
rz(-1.3277162) q[2];
sx q[2];
rz(1.1532016) q[2];
rz(0.48405805) q[3];
sx q[3];
rz(-2.3279133) q[3];
sx q[3];
rz(-2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83950481) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(-0.50305811) q[0];
rz(-1.5867651) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(-2.9876626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4707697) q[0];
sx q[0];
rz(-1.5564859) q[0];
sx q[0];
rz(-1.0679246) q[0];
rz(-pi) q[1];
rz(-1.1459848) q[2];
sx q[2];
rz(-0.69485352) q[2];
sx q[2];
rz(-1.4171464) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.82984) q[1];
sx q[1];
rz(-2.0318188) q[1];
sx q[1];
rz(-0.47383576) q[1];
rz(-pi) q[2];
rz(-0.50977317) q[3];
sx q[3];
rz(-1.5700649) q[3];
sx q[3];
rz(-0.36611205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8872035) q[2];
sx q[2];
rz(-2.7837191) q[2];
sx q[2];
rz(1.3228234) q[2];
rz(1.4860738) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(-2.4310908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.1948497) q[0];
sx q[0];
rz(-1.1179593) q[0];
sx q[0];
rz(2.2316566) q[0];
rz(0.77727708) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(-2.1562703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1811718) q[0];
sx q[0];
rz(-0.67502484) q[0];
sx q[0];
rz(-1.8280562) q[0];
x q[1];
rz(1.9820205) q[2];
sx q[2];
rz(-2.4422438) q[2];
sx q[2];
rz(-1.5027836) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4640376) q[1];
sx q[1];
rz(-2.8863393) q[1];
sx q[1];
rz(2.9222701) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3722234) q[3];
sx q[3];
rz(-2.3017075) q[3];
sx q[3];
rz(-1.4589256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.862792) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(-1.8939691) q[2];
rz(-0.4425846) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(-2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.9001532) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-2.0181657) q[0];
rz(2.5627047) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(-1.3935864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0279562) q[0];
sx q[0];
rz(-0.96158577) q[0];
sx q[0];
rz(-1.05818) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4030928) q[2];
sx q[2];
rz(-1.1412914) q[2];
sx q[2];
rz(1.6657366) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.046603831) q[1];
sx q[1];
rz(-2.2294083) q[1];
sx q[1];
rz(2.1556426) q[1];
x q[2];
rz(2.9080503) q[3];
sx q[3];
rz(-0.84905784) q[3];
sx q[3];
rz(-2.8535709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1018155) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(0.56143108) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(2.5533365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(1.2258688) q[0];
rz(1.6745802) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(1.3668758) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5850692) q[0];
sx q[0];
rz(-1.9559304) q[0];
sx q[0];
rz(-1.3655846) q[0];
x q[1];
rz(-1.5405802) q[2];
sx q[2];
rz(-2.317791) q[2];
sx q[2];
rz(-3.1291762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.67147672) q[1];
sx q[1];
rz(-1.6504382) q[1];
sx q[1];
rz(0.86298841) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59646888) q[3];
sx q[3];
rz(-0.58613741) q[3];
sx q[3];
rz(-0.49367192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86429578) q[2];
sx q[2];
rz(-2.0531451) q[2];
sx q[2];
rz(-0.40536353) q[2];
rz(0.46164414) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(1.595114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6500403) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(2.6382085) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(0.65178451) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0359905) q[0];
sx q[0];
rz(-2.833617) q[0];
sx q[0];
rz(2.4686747) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0405365) q[2];
sx q[2];
rz(-0.82436845) q[2];
sx q[2];
rz(-1.9161759) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1106907) q[1];
sx q[1];
rz(-1.4363524) q[1];
sx q[1];
rz(3.0893374) q[1];
x q[2];
rz(1.6739453) q[3];
sx q[3];
rz(-1.8522989) q[3];
sx q[3];
rz(-0.75479773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6289604) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(-0.39166489) q[2];
rz(2.9351249) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(-2.4519043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7713292) q[0];
sx q[0];
rz(-1.4498793) q[0];
sx q[0];
rz(2.1719334) q[0];
rz(-0.58352739) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-0.13024174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1657432) q[0];
sx q[0];
rz(-1.9231057) q[0];
sx q[0];
rz(1.6238814) q[0];
rz(-0.98353705) q[2];
sx q[2];
rz(-2.1526255) q[2];
sx q[2];
rz(-2.8772417) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.817037) q[1];
sx q[1];
rz(-0.55700028) q[1];
sx q[1];
rz(-0.42028285) q[1];
x q[2];
rz(0.98826615) q[3];
sx q[3];
rz(-1.2417955) q[3];
sx q[3];
rz(-2.695801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6234201) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(2.0557892) q[2];
rz(-3.0893677) q[3];
sx q[3];
rz(-1.4445392) q[3];
sx q[3];
rz(1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6770342) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(-1.1664671) q[0];
rz(-2.4160066) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(-0.74277791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9184473) q[0];
sx q[0];
rz(-0.76876516) q[0];
sx q[0];
rz(-0.2785152) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5119945) q[2];
sx q[2];
rz(-2.3592735) q[2];
sx q[2];
rz(0.19537374) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2837977) q[1];
sx q[1];
rz(-2.2504914) q[1];
sx q[1];
rz(0.53417511) q[1];
x q[2];
rz(-0.86117427) q[3];
sx q[3];
rz(-2.3434601) q[3];
sx q[3];
rz(0.29336151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5143738) q[2];
sx q[2];
rz(-1.6122931) q[2];
sx q[2];
rz(0.26088866) q[2];
rz(-2.0339113) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(-1.0816983) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35659197) q[0];
sx q[0];
rz(-0.78105015) q[0];
sx q[0];
rz(-2.2055431) q[0];
rz(-3.127457) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(-2.4900808) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0148841) q[0];
sx q[0];
rz(-1.5944905) q[0];
sx q[0];
rz(3.1372877) q[0];
rz(-pi) q[1];
rz(-0.75479836) q[2];
sx q[2];
rz(-2.3279466) q[2];
sx q[2];
rz(-2.1458643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.028101746) q[1];
sx q[1];
rz(-2.7999702) q[1];
sx q[1];
rz(1.2042868) q[1];
rz(-1.2048079) q[3];
sx q[3];
rz(-1.1567187) q[3];
sx q[3];
rz(0.28596349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8462048) q[2];
sx q[2];
rz(-2.442895) q[2];
sx q[2];
rz(0.52337581) q[2];
rz(0.30803099) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(1.3380922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.4073407) q[0];
sx q[0];
rz(-2.9946406) q[0];
sx q[0];
rz(-1.403632) q[0];
rz(2.667528) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(-1.3778936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5869851) q[0];
sx q[0];
rz(-2.4025612) q[0];
sx q[0];
rz(-0.87297312) q[0];
rz(-pi) q[1];
rz(-2.2950298) q[2];
sx q[2];
rz(-0.91869527) q[2];
sx q[2];
rz(-2.6400499) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4442515) q[1];
sx q[1];
rz(-1.8129983) q[1];
sx q[1];
rz(0.40476207) q[1];
rz(-pi) q[2];
x q[2];
rz(2.58425) q[3];
sx q[3];
rz(-1.7094269) q[3];
sx q[3];
rz(-2.0687452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3906117) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.6101458) q[2];
rz(-1.6837998) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(-2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67125852) q[0];
sx q[0];
rz(-2.8699734) q[0];
sx q[0];
rz(1.097453) q[0];
rz(0.63256565) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(-0.40614265) q[2];
sx q[2];
rz(-2.4293025) q[2];
sx q[2];
rz(0.14355125) q[2];
rz(-0.74606568) q[3];
sx q[3];
rz(-0.88701556) q[3];
sx q[3];
rz(-1.5942667) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];