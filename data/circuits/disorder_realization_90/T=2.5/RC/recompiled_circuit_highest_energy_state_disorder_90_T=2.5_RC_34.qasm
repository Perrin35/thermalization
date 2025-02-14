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
rz(-2.3174602) q[0];
sx q[0];
rz(-1.1622518) q[0];
sx q[0];
rz(2.9293769) q[0];
rz(-2.4556887) q[1];
sx q[1];
rz(-0.75690126) q[1];
sx q[1];
rz(-0.25605717) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60535359) q[0];
sx q[0];
rz(-2.1383717) q[0];
sx q[0];
rz(2.6044659) q[0];
rz(-pi) q[1];
rz(0.08698626) q[2];
sx q[2];
rz(-1.0115919) q[2];
sx q[2];
rz(2.5285697) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5642173) q[1];
sx q[1];
rz(-1.7322473) q[1];
sx q[1];
rz(-0.20906469) q[1];
rz(-2.831826) q[3];
sx q[3];
rz(-1.0599146) q[3];
sx q[3];
rz(2.8155308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.889692) q[2];
sx q[2];
rz(-2.8163741) q[2];
sx q[2];
rz(-1.088781) q[2];
rz(-0.47993663) q[3];
sx q[3];
rz(-1.9710541) q[3];
sx q[3];
rz(-1.1133194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8460409) q[0];
sx q[0];
rz(-0.22838455) q[0];
sx q[0];
rz(2.8725655) q[0];
rz(-2.524952) q[1];
sx q[1];
rz(-2.4863305) q[1];
sx q[1];
rz(0.90589595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30681102) q[0];
sx q[0];
rz(-2.3521345) q[0];
sx q[0];
rz(-0.10460268) q[0];
x q[1];
rz(-0.53585662) q[2];
sx q[2];
rz(-1.1108221) q[2];
sx q[2];
rz(-0.64086781) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91859282) q[1];
sx q[1];
rz(-1.4529374) q[1];
sx q[1];
rz(-2.8840114) q[1];
x q[2];
rz(-1.9578827) q[3];
sx q[3];
rz(-2.0969694) q[3];
sx q[3];
rz(-2.5691606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1632094) q[2];
sx q[2];
rz(-2.3613112) q[2];
sx q[2];
rz(0.015446375) q[2];
rz(-0.77203006) q[3];
sx q[3];
rz(-2.6431712) q[3];
sx q[3];
rz(1.8933403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60927272) q[0];
sx q[0];
rz(-2.5072704) q[0];
sx q[0];
rz(2.3442205) q[0];
rz(0.66870052) q[1];
sx q[1];
rz(-1.0182321) q[1];
sx q[1];
rz(2.5858137) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41125339) q[0];
sx q[0];
rz(-0.78124917) q[0];
sx q[0];
rz(-2.8528105) q[0];
rz(1.749281) q[2];
sx q[2];
rz(-0.67214078) q[2];
sx q[2];
rz(-1.3591421) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2712471) q[1];
sx q[1];
rz(-1.5588798) q[1];
sx q[1];
rz(0.72626497) q[1];
x q[2];
rz(-0.45826073) q[3];
sx q[3];
rz(-2.3790559) q[3];
sx q[3];
rz(-2.6242695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7734311) q[2];
sx q[2];
rz(-1.1661466) q[2];
sx q[2];
rz(-1.9318761) q[2];
rz(-3.0220253) q[3];
sx q[3];
rz(-0.60496324) q[3];
sx q[3];
rz(-0.97924489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6794353) q[0];
sx q[0];
rz(-1.698751) q[0];
sx q[0];
rz(0.47732842) q[0];
rz(0.10074549) q[1];
sx q[1];
rz(-1.0557749) q[1];
sx q[1];
rz(-0.13564067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0207708) q[0];
sx q[0];
rz(-0.77540708) q[0];
sx q[0];
rz(0.94803079) q[0];
x q[1];
rz(1.863998) q[2];
sx q[2];
rz(-2.0005156) q[2];
sx q[2];
rz(1.7079084) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1592334) q[1];
sx q[1];
rz(-1.9894588) q[1];
sx q[1];
rz(1.406698) q[1];
rz(-pi) q[2];
rz(0.78285406) q[3];
sx q[3];
rz(-1.4042036) q[3];
sx q[3];
rz(1.596791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9292844) q[2];
sx q[2];
rz(-1.4886651) q[2];
sx q[2];
rz(-2.6279214) q[2];
rz(-0.16289991) q[3];
sx q[3];
rz(-2.8859911) q[3];
sx q[3];
rz(-1.1027078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5279919) q[0];
sx q[0];
rz(-0.010471424) q[0];
sx q[0];
rz(-0.42798671) q[0];
rz(1.9635268) q[1];
sx q[1];
rz(-1.6798429) q[1];
sx q[1];
rz(1.4129432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75668272) q[0];
sx q[0];
rz(-1.9506978) q[0];
sx q[0];
rz(-1.0713473) q[0];
rz(-3.1052633) q[2];
sx q[2];
rz(-1.6128572) q[2];
sx q[2];
rz(-1.0147699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.08623422) q[1];
sx q[1];
rz(-1.5155122) q[1];
sx q[1];
rz(0.30200821) q[1];
rz(-pi) q[2];
rz(-0.67107646) q[3];
sx q[3];
rz(-0.6830712) q[3];
sx q[3];
rz(2.4623614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54521769) q[2];
sx q[2];
rz(-1.7866106) q[2];
sx q[2];
rz(0.4893111) q[2];
rz(-0.69928586) q[3];
sx q[3];
rz(-2.0977061) q[3];
sx q[3];
rz(-1.6635241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74772239) q[0];
sx q[0];
rz(-1.007653) q[0];
sx q[0];
rz(-1.7973768) q[0];
rz(-0.6952374) q[1];
sx q[1];
rz(-1.8724915) q[1];
sx q[1];
rz(-0.42323798) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3812034) q[0];
sx q[0];
rz(-1.4013327) q[0];
sx q[0];
rz(-1.5080835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6793785) q[2];
sx q[2];
rz(-0.75410289) q[2];
sx q[2];
rz(-1.6414409) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5343439) q[1];
sx q[1];
rz(-0.75466506) q[1];
sx q[1];
rz(0.21790803) q[1];
rz(2.9598992) q[3];
sx q[3];
rz(-1.5156997) q[3];
sx q[3];
rz(0.7167992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77373475) q[2];
sx q[2];
rz(-0.23907842) q[2];
sx q[2];
rz(1.9284922) q[2];
rz(-0.086094543) q[3];
sx q[3];
rz(-1.8910858) q[3];
sx q[3];
rz(-1.7566173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4940779) q[0];
sx q[0];
rz(-0.32951847) q[0];
sx q[0];
rz(-2.9373017) q[0];
rz(-0.34945166) q[1];
sx q[1];
rz(-1.5710257) q[1];
sx q[1];
rz(1.516516) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7766315) q[0];
sx q[0];
rz(-1.6775515) q[0];
sx q[0];
rz(-2.0364709) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9728685) q[2];
sx q[2];
rz(-1.8622145) q[2];
sx q[2];
rz(0.60044392) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8015281) q[1];
sx q[1];
rz(-1.6710647) q[1];
sx q[1];
rz(3.0951436) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99269331) q[3];
sx q[3];
rz(-1.5823964) q[3];
sx q[3];
rz(1.8314401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1427631) q[2];
sx q[2];
rz(-0.62935722) q[2];
sx q[2];
rz(-2.2289185) q[2];
rz(-3.1257889) q[3];
sx q[3];
rz(-2.0533357) q[3];
sx q[3];
rz(0.44939941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939746) q[0];
sx q[0];
rz(-2.0517218) q[0];
sx q[0];
rz(0.81899482) q[0];
rz(0.35797572) q[1];
sx q[1];
rz(-2.0134108) q[1];
sx q[1];
rz(0.63041675) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6651745) q[0];
sx q[0];
rz(-0.074284241) q[0];
sx q[0];
rz(-1.7887542) q[0];
rz(2.6825601) q[2];
sx q[2];
rz(-0.9046281) q[2];
sx q[2];
rz(-0.55202548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9220029) q[1];
sx q[1];
rz(-2.1685984) q[1];
sx q[1];
rz(-2.647318) q[1];
rz(-pi) q[2];
rz(-2.1958293) q[3];
sx q[3];
rz(-1.9595567) q[3];
sx q[3];
rz(0.22076535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7716052) q[2];
sx q[2];
rz(-2.6321754) q[2];
sx q[2];
rz(2.7484757) q[2];
rz(-2.3734132) q[3];
sx q[3];
rz(-2.3438175) q[3];
sx q[3];
rz(1.2099077) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2792252) q[0];
sx q[0];
rz(-2.6287862) q[0];
sx q[0];
rz(-0.38238907) q[0];
rz(0.44876107) q[1];
sx q[1];
rz(-0.02350137) q[1];
sx q[1];
rz(2.1429367) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6210839) q[0];
sx q[0];
rz(-1.3338085) q[0];
sx q[0];
rz(-2.468796) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.08491718) q[2];
sx q[2];
rz(-1.907427) q[2];
sx q[2];
rz(-1.4648598) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4289394) q[1];
sx q[1];
rz(-1.3725159) q[1];
sx q[1];
rz(2.8444931) q[1];
x q[2];
rz(2.4527384) q[3];
sx q[3];
rz(-1.5176532) q[3];
sx q[3];
rz(-2.1518858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0465595) q[2];
sx q[2];
rz(-2.3456489) q[2];
sx q[2];
rz(-2.156192) q[2];
rz(2.4424148) q[3];
sx q[3];
rz(-2.2907084) q[3];
sx q[3];
rz(-0.053475577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92626524) q[0];
sx q[0];
rz(-0.43020058) q[0];
sx q[0];
rz(2.3660124) q[0];
rz(2.8714478) q[1];
sx q[1];
rz(-1.6855449) q[1];
sx q[1];
rz(0.57033479) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71513373) q[0];
sx q[0];
rz(-1.1021779) q[0];
sx q[0];
rz(0.56632407) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.774046) q[2];
sx q[2];
rz(-2.654944) q[2];
sx q[2];
rz(-1.1543076) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66460863) q[1];
sx q[1];
rz(-0.23566569) q[1];
sx q[1];
rz(0.13500555) q[1];
x q[2];
rz(-1.4943284) q[3];
sx q[3];
rz(-0.88275331) q[3];
sx q[3];
rz(-2.9043759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21738416) q[2];
sx q[2];
rz(-1.1576098) q[2];
sx q[2];
rz(2.6317934) q[2];
rz(-0.23586759) q[3];
sx q[3];
rz(-0.1736975) q[3];
sx q[3];
rz(-2.164446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8248642) q[0];
sx q[0];
rz(-1.7289799) q[0];
sx q[0];
rz(-1.254542) q[0];
rz(-2.6084689) q[1];
sx q[1];
rz(-1.6785379) q[1];
sx q[1];
rz(1.0482845) q[1];
rz(1.6021846) q[2];
sx q[2];
rz(-2.7161408) q[2];
sx q[2];
rz(-1.3795992) q[2];
rz(0.20797603) q[3];
sx q[3];
rz(-2.6317876) q[3];
sx q[3];
rz(-0.77710487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
