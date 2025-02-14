OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.00081113022) q[0];
sx q[0];
rz(-0.82511628) q[0];
sx q[0];
rz(-2.7197279) q[0];
rz(2.3946664) q[1];
sx q[1];
rz(-0.59023017) q[1];
sx q[1];
rz(-2.3205369) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29971545) q[0];
sx q[0];
rz(-1.4869191) q[0];
sx q[0];
rz(-1.6158577) q[0];
rz(1.3652322) q[2];
sx q[2];
rz(-1.7920307) q[2];
sx q[2];
rz(1.4874171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3154308) q[1];
sx q[1];
rz(-1.9524442) q[1];
sx q[1];
rz(0.30604575) q[1];
rz(-pi) q[2];
rz(3.0077137) q[3];
sx q[3];
rz(-1.9124096) q[3];
sx q[3];
rz(-2.435129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39516285) q[2];
sx q[2];
rz(-0.91150993) q[2];
sx q[2];
rz(2.5968623) q[2];
rz(-1.4970477) q[3];
sx q[3];
rz(-1.112273) q[3];
sx q[3];
rz(1.8018319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4561975) q[0];
sx q[0];
rz(-0.74390262) q[0];
sx q[0];
rz(-2.5681382) q[0];
rz(3.0990797) q[1];
sx q[1];
rz(-1.3923693) q[1];
sx q[1];
rz(-1.858985) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79270482) q[0];
sx q[0];
rz(-0.12785873) q[0];
sx q[0];
rz(-0.16221817) q[0];
x q[1];
rz(-2.4290724) q[2];
sx q[2];
rz(-2.7606568) q[2];
sx q[2];
rz(-2.088758) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.01067082) q[1];
sx q[1];
rz(-2.0951443) q[1];
sx q[1];
rz(1.4078543) q[1];
rz(-pi) q[2];
rz(-2.1759791) q[3];
sx q[3];
rz(-2.11016) q[3];
sx q[3];
rz(2.2210647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87806988) q[2];
sx q[2];
rz(-1.1493378) q[2];
sx q[2];
rz(-0.44352201) q[2];
rz(1.589132) q[3];
sx q[3];
rz(-2.7510721) q[3];
sx q[3];
rz(0.21903567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1294915) q[0];
sx q[0];
rz(-0.87915593) q[0];
sx q[0];
rz(0.40170676) q[0];
rz(-1.4178287) q[1];
sx q[1];
rz(-0.44770733) q[1];
sx q[1];
rz(-2.0254501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4075466) q[0];
sx q[0];
rz(-0.73483682) q[0];
sx q[0];
rz(-0.46291344) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1537519) q[2];
sx q[2];
rz(-2.2567686) q[2];
sx q[2];
rz(-1.5123715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9520719) q[1];
sx q[1];
rz(-2.4671989) q[1];
sx q[1];
rz(0.86878573) q[1];
rz(-pi) q[2];
rz(2.5930034) q[3];
sx q[3];
rz(-1.635713) q[3];
sx q[3];
rz(1.3260076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9125354) q[2];
sx q[2];
rz(-2.8450862) q[2];
sx q[2];
rz(2.1336446) q[2];
rz(1.5063162) q[3];
sx q[3];
rz(-1.9624036) q[3];
sx q[3];
rz(0.87872163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85663831) q[0];
sx q[0];
rz(-3.1143739) q[0];
sx q[0];
rz(-2.304049) q[0];
rz(-1.8800927) q[1];
sx q[1];
rz(-1.7439758) q[1];
sx q[1];
rz(1.8998442) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1364131) q[0];
sx q[0];
rz(-1.2997928) q[0];
sx q[0];
rz(-0.47309978) q[0];
rz(-2.1793887) q[2];
sx q[2];
rz(-1.7458916) q[2];
sx q[2];
rz(2.1661557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95835244) q[1];
sx q[1];
rz(-1.1179525) q[1];
sx q[1];
rz(2.0264105) q[1];
x q[2];
rz(2.1513274) q[3];
sx q[3];
rz(-1.255569) q[3];
sx q[3];
rz(0.48993708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9813098) q[2];
sx q[2];
rz(-2.1872988) q[2];
sx q[2];
rz(-1.0120288) q[2];
rz(-2.1435598) q[3];
sx q[3];
rz(-2.4055552) q[3];
sx q[3];
rz(-1.6691104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0618458) q[0];
sx q[0];
rz(-2.340402) q[0];
sx q[0];
rz(2.8629942) q[0];
rz(0.31696907) q[1];
sx q[1];
rz(-1.7520889) q[1];
sx q[1];
rz(0.46151361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6742548) q[0];
sx q[0];
rz(-1.6188356) q[0];
sx q[0];
rz(1.4018628) q[0];
rz(1.0217083) q[2];
sx q[2];
rz(-2.4442992) q[2];
sx q[2];
rz(1.0233699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0089142) q[1];
sx q[1];
rz(-1.3323093) q[1];
sx q[1];
rz(-3.0452252) q[1];
x q[2];
rz(-1.0255542) q[3];
sx q[3];
rz(-1.488113) q[3];
sx q[3];
rz(2.290639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3658112) q[2];
sx q[2];
rz(-1.8041939) q[2];
sx q[2];
rz(-2.865045) q[2];
rz(0.60025275) q[3];
sx q[3];
rz(-0.7163896) q[3];
sx q[3];
rz(2.6071809) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14715956) q[0];
sx q[0];
rz(-2.5614547) q[0];
sx q[0];
rz(0.68122) q[0];
rz(0.45411202) q[1];
sx q[1];
rz(-1.157016) q[1];
sx q[1];
rz(-0.62201321) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21888622) q[0];
sx q[0];
rz(-1.6657296) q[0];
sx q[0];
rz(-2.9430998) q[0];
rz(-pi) q[1];
rz(-0.76592457) q[2];
sx q[2];
rz(-0.799338) q[2];
sx q[2];
rz(0.23169416) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4677538) q[1];
sx q[1];
rz(-2.2430393) q[1];
sx q[1];
rz(-0.52947395) q[1];
rz(-pi) q[2];
rz(-2.2174453) q[3];
sx q[3];
rz(-2.7139276) q[3];
sx q[3];
rz(-0.88811737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5155718) q[2];
sx q[2];
rz(-0.61352789) q[2];
sx q[2];
rz(2.386509) q[2];
rz(0.10609047) q[3];
sx q[3];
rz(-1.2664436) q[3];
sx q[3];
rz(-0.46060002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0945011) q[0];
sx q[0];
rz(-2.5300808) q[0];
sx q[0];
rz(3.0754572) q[0];
rz(0.051008929) q[1];
sx q[1];
rz(-0.74791932) q[1];
sx q[1];
rz(-1.8775108) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8569778) q[0];
sx q[0];
rz(-0.98941411) q[0];
sx q[0];
rz(-2.7063497) q[0];
rz(-0.71526975) q[2];
sx q[2];
rz(-1.8085305) q[2];
sx q[2];
rz(-2.0572853) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3466929) q[1];
sx q[1];
rz(-2.3155619) q[1];
sx q[1];
rz(-0.54462437) q[1];
rz(-0.73406666) q[3];
sx q[3];
rz(-2.5841004) q[3];
sx q[3];
rz(-2.8898784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2833726) q[2];
sx q[2];
rz(-2.0141352) q[2];
sx q[2];
rz(-0.14632012) q[2];
rz(2.537651) q[3];
sx q[3];
rz(-2.5950409) q[3];
sx q[3];
rz(-1.4124136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42501763) q[0];
sx q[0];
rz(-1.9442433) q[0];
sx q[0];
rz(3.0498411) q[0];
rz(-3.0692406) q[1];
sx q[1];
rz(-1.3183343) q[1];
sx q[1];
rz(-0.66973698) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2429072) q[0];
sx q[0];
rz(-1.7062441) q[0];
sx q[0];
rz(2.8464635) q[0];
x q[1];
rz(2.8477741) q[2];
sx q[2];
rz(-2.71596) q[2];
sx q[2];
rz(0.64070492) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9409166) q[1];
sx q[1];
rz(-1.7155316) q[1];
sx q[1];
rz(-0.30818589) q[1];
rz(-2.0132695) q[3];
sx q[3];
rz(-1.4391581) q[3];
sx q[3];
rz(-1.4030171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79121315) q[2];
sx q[2];
rz(-2.4243441) q[2];
sx q[2];
rz(0.72163248) q[2];
rz(2.3676938) q[3];
sx q[3];
rz(-0.98267233) q[3];
sx q[3];
rz(-0.53988808) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4197107) q[0];
sx q[0];
rz(-1.9419436) q[0];
sx q[0];
rz(2.6019959) q[0];
rz(-2.3609912) q[1];
sx q[1];
rz(-1.2465979) q[1];
sx q[1];
rz(2.7194729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0151505) q[0];
sx q[0];
rz(-1.0202304) q[0];
sx q[0];
rz(0.88811876) q[0];
rz(-2.6494041) q[2];
sx q[2];
rz(-1.5180615) q[2];
sx q[2];
rz(-1.9143788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.064116) q[1];
sx q[1];
rz(-1.6088652) q[1];
sx q[1];
rz(-1.0910973) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5255181) q[3];
sx q[3];
rz(-2.1589734) q[3];
sx q[3];
rz(-2.8578293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.098794) q[2];
sx q[2];
rz(-2.1197987) q[2];
sx q[2];
rz(-2.5928024) q[2];
rz(0.15268606) q[3];
sx q[3];
rz(-1.0278253) q[3];
sx q[3];
rz(2.4299183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3720836) q[0];
sx q[0];
rz(-2.7476269) q[0];
sx q[0];
rz(2.5373996) q[0];
rz(1.4671885) q[1];
sx q[1];
rz(-1.7908275) q[1];
sx q[1];
rz(1.9942572) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7580099) q[0];
sx q[0];
rz(-0.95057633) q[0];
sx q[0];
rz(-1.0538641) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34530039) q[2];
sx q[2];
rz(-0.7658813) q[2];
sx q[2];
rz(-1.6102546) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5529788) q[1];
sx q[1];
rz(-2.0535894) q[1];
sx q[1];
rz(-1.5534414) q[1];
x q[2];
rz(2.2828034) q[3];
sx q[3];
rz(-0.50903532) q[3];
sx q[3];
rz(1.7767752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6818105) q[2];
sx q[2];
rz(-1.1638389) q[2];
sx q[2];
rz(-0.60010827) q[2];
rz(2.4188304) q[3];
sx q[3];
rz(-2.1434651) q[3];
sx q[3];
rz(0.37989894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.4031274) q[0];
sx q[0];
rz(-1.6275788) q[0];
sx q[0];
rz(1.7049261) q[0];
rz(1.5557095) q[1];
sx q[1];
rz(-1.7441505) q[1];
sx q[1];
rz(2.2081262) q[1];
rz(2.6513097) q[2];
sx q[2];
rz(-0.96302196) q[2];
sx q[2];
rz(-1.8204126) q[2];
rz(-2.4248471) q[3];
sx q[3];
rz(-2.8195753) q[3];
sx q[3];
rz(-1.4618518) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
