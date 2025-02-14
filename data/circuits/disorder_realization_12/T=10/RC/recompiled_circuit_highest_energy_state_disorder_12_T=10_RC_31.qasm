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
rz(2.501261) q[0];
sx q[0];
rz(-2.2202272) q[0];
sx q[0];
rz(1.5322354) q[0];
rz(0.20707239) q[1];
sx q[1];
rz(4.3011811) q[1];
sx q[1];
rz(11.094697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4764413) q[0];
sx q[0];
rz(-1.8855321) q[0];
sx q[0];
rz(1.9010309) q[0];
x q[1];
rz(0.81469131) q[2];
sx q[2];
rz(-1.0202131) q[2];
sx q[2];
rz(1.4514635) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67480278) q[1];
sx q[1];
rz(-0.25898749) q[1];
sx q[1];
rz(-1.966443) q[1];
rz(-pi) q[2];
rz(-0.26443414) q[3];
sx q[3];
rz(-0.61226058) q[3];
sx q[3];
rz(-2.0443288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6472935) q[2];
sx q[2];
rz(-0.99227253) q[2];
sx q[2];
rz(-0.54599071) q[2];
rz(0.93849385) q[3];
sx q[3];
rz(-0.5135082) q[3];
sx q[3];
rz(-2.688496) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12478011) q[0];
sx q[0];
rz(-1.1174959) q[0];
sx q[0];
rz(2.5002531) q[0];
rz(-1.1891787) q[1];
sx q[1];
rz(-2.7180505) q[1];
sx q[1];
rz(2.0372527) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46342349) q[0];
sx q[0];
rz(-2.1974653) q[0];
sx q[0];
rz(-1.3889484) q[0];
rz(-1.0510089) q[2];
sx q[2];
rz(-0.6001282) q[2];
sx q[2];
rz(2.5243197) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84766372) q[1];
sx q[1];
rz(-2.1184455) q[1];
sx q[1];
rz(2.6002162) q[1];
rz(-1.6112174) q[3];
sx q[3];
rz(-1.6580402) q[3];
sx q[3];
rz(0.62007346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7426593) q[2];
sx q[2];
rz(-1.089596) q[2];
sx q[2];
rz(0.99417865) q[2];
rz(1.5587156) q[3];
sx q[3];
rz(-2.4623058) q[3];
sx q[3];
rz(2.5905632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0068552103) q[0];
sx q[0];
rz(-2.0515433) q[0];
sx q[0];
rz(2.5110974) q[0];
rz(-2.9948803) q[1];
sx q[1];
rz(-1.881733) q[1];
sx q[1];
rz(-2.2781118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1711581) q[0];
sx q[0];
rz(-2.2169431) q[0];
sx q[0];
rz(-2.7617654) q[0];
rz(-1.2892154) q[2];
sx q[2];
rz(-1.3669434) q[2];
sx q[2];
rz(-1.2512164) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5456704) q[1];
sx q[1];
rz(-1.1356259) q[1];
sx q[1];
rz(0.11637139) q[1];
x q[2];
rz(2.9894838) q[3];
sx q[3];
rz(-2.8242691) q[3];
sx q[3];
rz(-0.0025686669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33519393) q[2];
sx q[2];
rz(-1.6890182) q[2];
sx q[2];
rz(-2.2550968) q[2];
rz(2.7116306) q[3];
sx q[3];
rz(-2.6468247) q[3];
sx q[3];
rz(2.2494242) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48915136) q[0];
sx q[0];
rz(-2.8471071) q[0];
sx q[0];
rz(2.6926706) q[0];
rz(-0.65662518) q[1];
sx q[1];
rz(-1.7078327) q[1];
sx q[1];
rz(0.33033672) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8632148) q[0];
sx q[0];
rz(-1.2513573) q[0];
sx q[0];
rz(-2.9055975) q[0];
rz(-1.0895715) q[2];
sx q[2];
rz(-2.6767566) q[2];
sx q[2];
rz(-0.46457738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9724839) q[1];
sx q[1];
rz(-2.0044998) q[1];
sx q[1];
rz(0.32167158) q[1];
rz(-3.1325601) q[3];
sx q[3];
rz(-1.6918283) q[3];
sx q[3];
rz(-1.5980522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0535447) q[2];
sx q[2];
rz(-1.6036754) q[2];
sx q[2];
rz(-0.21928445) q[2];
rz(0.80026475) q[3];
sx q[3];
rz(-2.181668) q[3];
sx q[3];
rz(2.3805591) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4735755) q[0];
sx q[0];
rz(-2.0125084) q[0];
sx q[0];
rz(-1.4833204) q[0];
rz(-2.5738916) q[1];
sx q[1];
rz(-1.9300108) q[1];
sx q[1];
rz(-1.6426881) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2788102) q[0];
sx q[0];
rz(-1.4238951) q[0];
sx q[0];
rz(0.53333078) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0437356) q[2];
sx q[2];
rz(-1.5817533) q[2];
sx q[2];
rz(-1.1681854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3659542) q[1];
sx q[1];
rz(-1.3729551) q[1];
sx q[1];
rz(-0.61590804) q[1];
x q[2];
rz(-0.12171774) q[3];
sx q[3];
rz(-1.7608402) q[3];
sx q[3];
rz(0.42480436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.5119001) q[2];
sx q[2];
rz(-2.3834159) q[2];
sx q[2];
rz(0.32621113) q[2];
rz(-0.1263667) q[3];
sx q[3];
rz(-2.5765403) q[3];
sx q[3];
rz(2.8700184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7997953) q[0];
sx q[0];
rz(-1.3038776) q[0];
sx q[0];
rz(-0.52249348) q[0];
rz(2.3937468) q[1];
sx q[1];
rz(-1.5380842) q[1];
sx q[1];
rz(-1.1572256) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3078454) q[0];
sx q[0];
rz(-2.5933449) q[0];
sx q[0];
rz(-0.98497106) q[0];
x q[1];
rz(2.0339478) q[2];
sx q[2];
rz(-1.9310631) q[2];
sx q[2];
rz(2.2778794) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.044837601) q[1];
sx q[1];
rz(-1.6628578) q[1];
sx q[1];
rz(-1.7141059) q[1];
rz(-pi) q[2];
rz(-1.4943407) q[3];
sx q[3];
rz(-1.8969733) q[3];
sx q[3];
rz(-0.13274176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.010765643) q[2];
sx q[2];
rz(-0.5126493) q[2];
sx q[2];
rz(-1.1831247) q[2];
rz(-0.042923953) q[3];
sx q[3];
rz(-1.5019417) q[3];
sx q[3];
rz(-2.8973798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6965028) q[0];
sx q[0];
rz(-0.88560167) q[0];
sx q[0];
rz(-0.67286277) q[0];
rz(-0.54975763) q[1];
sx q[1];
rz(-1.8472698) q[1];
sx q[1];
rz(-1.6514282) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8880418) q[0];
sx q[0];
rz(-1.7490938) q[0];
sx q[0];
rz(-0.56714296) q[0];
x q[1];
rz(2.6308833) q[2];
sx q[2];
rz(-1.129653) q[2];
sx q[2];
rz(-0.06076755) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0800327) q[1];
sx q[1];
rz(-1.4747689) q[1];
sx q[1];
rz(-2.705234) q[1];
x q[2];
rz(2.1792322) q[3];
sx q[3];
rz(-0.058789805) q[3];
sx q[3];
rz(-0.99942452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.043639) q[2];
sx q[2];
rz(-0.54936469) q[2];
sx q[2];
rz(-1.1489541) q[2];
rz(-0.39135459) q[3];
sx q[3];
rz(-1.9214182) q[3];
sx q[3];
rz(2.2933188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5963762) q[0];
sx q[0];
rz(-1.2059809) q[0];
sx q[0];
rz(0.69860953) q[0];
rz(-0.85083234) q[1];
sx q[1];
rz(-2.5406676) q[1];
sx q[1];
rz(2.1934379) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20336172) q[0];
sx q[0];
rz(-1.6530418) q[0];
sx q[0];
rz(-2.4191678) q[0];
rz(-pi) q[1];
rz(-1.9004563) q[2];
sx q[2];
rz(-1.5200117) q[2];
sx q[2];
rz(0.18731681) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8906456) q[1];
sx q[1];
rz(-2.1275824) q[1];
sx q[1];
rz(2.6614038) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3513673) q[3];
sx q[3];
rz(-0.91013113) q[3];
sx q[3];
rz(-0.017221551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4725388) q[2];
sx q[2];
rz(-1.6927745) q[2];
sx q[2];
rz(-1.6834458) q[2];
rz(-2.0504047) q[3];
sx q[3];
rz(-0.89717054) q[3];
sx q[3];
rz(-2.9567961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75337306) q[0];
sx q[0];
rz(-0.18167697) q[0];
sx q[0];
rz(1.8411807) q[0];
rz(-0.45937195) q[1];
sx q[1];
rz(-1.6875024) q[1];
sx q[1];
rz(2.557911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.822264) q[0];
sx q[0];
rz(-2.2177025) q[0];
sx q[0];
rz(-2.1859593) q[0];
x q[1];
rz(-1.5351717) q[2];
sx q[2];
rz(-1.5322663) q[2];
sx q[2];
rz(2.5690479) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4286268) q[1];
sx q[1];
rz(-1.2199063) q[1];
sx q[1];
rz(2.116973) q[1];
rz(-1.0631869) q[3];
sx q[3];
rz(-0.48590966) q[3];
sx q[3];
rz(-2.8775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61648458) q[2];
sx q[2];
rz(-2.8664092) q[2];
sx q[2];
rz(-0.52213651) q[2];
rz(2.3790242) q[3];
sx q[3];
rz(-1.6430166) q[3];
sx q[3];
rz(-2.7975119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6793215) q[0];
sx q[0];
rz(-1.2079879) q[0];
sx q[0];
rz(2.5479877) q[0];
rz(0.69315928) q[1];
sx q[1];
rz(-2.3489372) q[1];
sx q[1];
rz(2.3902182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7369738) q[0];
sx q[0];
rz(-1.829379) q[0];
sx q[0];
rz(0.41249852) q[0];
rz(-0.61155106) q[2];
sx q[2];
rz(-2.1566506) q[2];
sx q[2];
rz(-0.73918804) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33022949) q[1];
sx q[1];
rz(-0.6682446) q[1];
sx q[1];
rz(0.10295566) q[1];
rz(-pi) q[2];
rz(0.42507024) q[3];
sx q[3];
rz(-0.18501013) q[3];
sx q[3];
rz(1.6744705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1410602) q[2];
sx q[2];
rz(-2.4324721) q[2];
sx q[2];
rz(0.023690311) q[2];
rz(-2.0247816) q[3];
sx q[3];
rz(-1.6938554) q[3];
sx q[3];
rz(-2.007133) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58917191) q[0];
sx q[0];
rz(-1.9235274) q[0];
sx q[0];
rz(-2.0148475) q[0];
rz(-1.8213656) q[1];
sx q[1];
rz(-1.8658493) q[1];
sx q[1];
rz(1.0125926) q[1];
rz(-2.7350551) q[2];
sx q[2];
rz(-2.5589522) q[2];
sx q[2];
rz(2.027579) q[2];
rz(0.22707012) q[3];
sx q[3];
rz(-1.3454701) q[3];
sx q[3];
rz(-3.0440192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
