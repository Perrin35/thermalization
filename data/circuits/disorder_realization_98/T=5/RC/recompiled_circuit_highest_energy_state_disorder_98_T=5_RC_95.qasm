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
rz(-0.67796081) q[0];
sx q[0];
rz(-0.27422658) q[0];
sx q[0];
rz(1.8507313) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(1.4996127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9737273) q[0];
sx q[0];
rz(-1.8059547) q[0];
sx q[0];
rz(-0.27882378) q[0];
x q[1];
rz(-2.847509) q[2];
sx q[2];
rz(-0.80162797) q[2];
sx q[2];
rz(-1.5844567) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5747524) q[1];
sx q[1];
rz(-0.87927239) q[1];
sx q[1];
rz(-1.1471466) q[1];
x q[2];
rz(0.94120195) q[3];
sx q[3];
rz(-1.968002) q[3];
sx q[3];
rz(-1.8091699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9700254) q[2];
sx q[2];
rz(-1.0415404) q[2];
sx q[2];
rz(-2.2402666) q[2];
rz(2.6775635) q[3];
sx q[3];
rz(-1.2814859) q[3];
sx q[3];
rz(0.06981167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0932015) q[0];
sx q[0];
rz(-0.036660107) q[0];
sx q[0];
rz(-1.8638336) q[0];
rz(0.094206421) q[1];
sx q[1];
rz(-0.46368805) q[1];
sx q[1];
rz(-1.6069848) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2814935) q[0];
sx q[0];
rz(-1.5380385) q[0];
sx q[0];
rz(0.97086914) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7530551) q[2];
sx q[2];
rz(-1.6267952) q[2];
sx q[2];
rz(0.1268498) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.822545) q[1];
sx q[1];
rz(-1.3101868) q[1];
sx q[1];
rz(-0.50995639) q[1];
rz(-1.9877795) q[3];
sx q[3];
rz(-1.9440042) q[3];
sx q[3];
rz(2.5905379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4216807) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(-2.9812532) q[2];
rz(-0.73873377) q[3];
sx q[3];
rz(-0.84638798) q[3];
sx q[3];
rz(1.4275449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.5147603) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(0.50977388) q[0];
rz(1.7104507) q[1];
sx q[1];
rz(-1.1809843) q[1];
sx q[1];
rz(-0.74657718) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11543848) q[0];
sx q[0];
rz(-2.5887449) q[0];
sx q[0];
rz(1.5325969) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4116729) q[2];
sx q[2];
rz(-1.0906719) q[2];
sx q[2];
rz(-1.4381222) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1521897) q[1];
sx q[1];
rz(-2.2446989) q[1];
sx q[1];
rz(2.6722355) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88991092) q[3];
sx q[3];
rz(-1.2632252) q[3];
sx q[3];
rz(2.6079082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15472445) q[2];
sx q[2];
rz(-2.8747323) q[2];
sx q[2];
rz(1.9179087) q[2];
rz(-1.0736505) q[3];
sx q[3];
rz(-1.4370388) q[3];
sx q[3];
rz(0.8980208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.8680442) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(-2.7622188) q[0];
rz(0.1768449) q[1];
sx q[1];
rz(-1.6601446) q[1];
sx q[1];
rz(-0.79536974) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5569755) q[0];
sx q[0];
rz(-1.4780227) q[0];
sx q[0];
rz(1.2014821) q[0];
rz(-1.8537117) q[2];
sx q[2];
rz(-2.3672315) q[2];
sx q[2];
rz(-1.1224358) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6713555) q[1];
sx q[1];
rz(-2.0908326) q[1];
sx q[1];
rz(-2.8440471) q[1];
x q[2];
rz(-1.5136817) q[3];
sx q[3];
rz(-1.4453672) q[3];
sx q[3];
rz(2.7386087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7103601) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(1.360652) q[2];
rz(-2.1121173) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(-1.3220538) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4904356) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(-1.5929476) q[0];
rz(-2.7410638) q[1];
sx q[1];
rz(-1.6533886) q[1];
sx q[1];
rz(-1.7318447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5142877) q[0];
sx q[0];
rz(-0.76111932) q[0];
sx q[0];
rz(-2.850015) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39938853) q[2];
sx q[2];
rz(-2.7542973) q[2];
sx q[2];
rz(-0.072810955) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61297902) q[1];
sx q[1];
rz(-0.76957146) q[1];
sx q[1];
rz(2.2320094) q[1];
x q[2];
rz(-1.665776) q[3];
sx q[3];
rz(-2.5174369) q[3];
sx q[3];
rz(0.74172663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3132402) q[2];
sx q[2];
rz(-1.4408377) q[2];
sx q[2];
rz(2.7102615) q[2];
rz(3.0865772) q[3];
sx q[3];
rz(-2.8300245) q[3];
sx q[3];
rz(0.55606786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738275) q[0];
sx q[0];
rz(-0.25280935) q[0];
sx q[0];
rz(1.025169) q[0];
rz(-0.45285666) q[1];
sx q[1];
rz(-2.5621474) q[1];
sx q[1];
rz(2.2377009) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4898281) q[0];
sx q[0];
rz(-1.8790885) q[0];
sx q[0];
rz(0.90109189) q[0];
x q[1];
rz(-2.507693) q[2];
sx q[2];
rz(-1.0154775) q[2];
sx q[2];
rz(-2.7125396) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0597983) q[1];
sx q[1];
rz(-0.79080338) q[1];
sx q[1];
rz(0.38500144) q[1];
rz(1.4056095) q[3];
sx q[3];
rz(-1.4017229) q[3];
sx q[3];
rz(-1.0195093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.30770939) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(0.21635381) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-2.4560865) q[3];
sx q[3];
rz(-1.640813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-1.1067156) q[0];
sx q[0];
rz(-0.69881451) q[0];
rz(-1.1837333) q[1];
sx q[1];
rz(-1.7203169) q[1];
sx q[1];
rz(0.85174495) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158902) q[0];
sx q[0];
rz(-0.53759241) q[0];
sx q[0];
rz(-1.6449152) q[0];
x q[1];
rz(2.5642407) q[2];
sx q[2];
rz(-1.5492348) q[2];
sx q[2];
rz(-0.93761629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5761583) q[1];
sx q[1];
rz(-1.9322188) q[1];
sx q[1];
rz(-2.4225077) q[1];
rz(-pi) q[2];
rz(0.64166358) q[3];
sx q[3];
rz(-1.9681276) q[3];
sx q[3];
rz(-1.404055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8375887) q[2];
sx q[2];
rz(-1.3112661) q[2];
sx q[2];
rz(-2.6336929) q[2];
rz(1.0287644) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(0.60104162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69865882) q[0];
sx q[0];
rz(-1.3154987) q[0];
sx q[0];
rz(-0.0096631924) q[0];
rz(1.6161605) q[1];
sx q[1];
rz(-1.7639672) q[1];
sx q[1];
rz(-0.38472167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98823901) q[0];
sx q[0];
rz(-1.8073449) q[0];
sx q[0];
rz(-0.1519712) q[0];
rz(-0.20236774) q[2];
sx q[2];
rz(-1.9500537) q[2];
sx q[2];
rz(2.8059562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5009969) q[1];
sx q[1];
rz(-1.847155) q[1];
sx q[1];
rz(-0.21548693) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58458768) q[3];
sx q[3];
rz(-0.86170022) q[3];
sx q[3];
rz(2.5392836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5821417) q[2];
sx q[2];
rz(-2.1820575) q[2];
sx q[2];
rz(3.1257296) q[2];
rz(-1.9150241) q[3];
sx q[3];
rz(-0.80569402) q[3];
sx q[3];
rz(-0.62172833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7568307) q[0];
sx q[0];
rz(-2.4767196) q[0];
sx q[0];
rz(-1.0913947) q[0];
rz(-0.81820828) q[1];
sx q[1];
rz(-1.810377) q[1];
sx q[1];
rz(2.4726726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92907897) q[0];
sx q[0];
rz(-2.2336279) q[0];
sx q[0];
rz(2.3431091) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0537655) q[2];
sx q[2];
rz(-0.7730128) q[2];
sx q[2];
rz(0.20394606) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3260169) q[1];
sx q[1];
rz(-1.6417802) q[1];
sx q[1];
rz(-2.1287588) q[1];
x q[2];
rz(2.8112765) q[3];
sx q[3];
rz(-1.013375) q[3];
sx q[3];
rz(-1.6176318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6341256) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(-2.1779306) q[2];
rz(-1.4108747) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(-1.3655519) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1249579) q[0];
sx q[0];
rz(-2.5732915) q[0];
sx q[0];
rz(1.4547263) q[0];
rz(1.1085054) q[1];
sx q[1];
rz(-1.6051555) q[1];
sx q[1];
rz(-0.32807168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1301918) q[0];
sx q[0];
rz(-1.6284124) q[0];
sx q[0];
rz(-1.4559697) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56888442) q[2];
sx q[2];
rz(-2.9279104) q[2];
sx q[2];
rz(1.3485707) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83392622) q[1];
sx q[1];
rz(-0.95912479) q[1];
sx q[1];
rz(2.0425379) q[1];
rz(-pi) q[2];
rz(2.1222018) q[3];
sx q[3];
rz(-0.20272045) q[3];
sx q[3];
rz(-0.77714506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46128094) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(-2.7247562) q[2];
rz(-0.69878116) q[3];
sx q[3];
rz(-0.67658934) q[3];
sx q[3];
rz(0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1998491) q[0];
sx q[0];
rz(-1.9244292) q[0];
sx q[0];
rz(1.695965) q[0];
rz(0.70882123) q[1];
sx q[1];
rz(-0.91239057) q[1];
sx q[1];
rz(-0.053587996) q[1];
rz(2.2075072) q[2];
sx q[2];
rz(-1.4082505) q[2];
sx q[2];
rz(0.55278964) q[2];
rz(0.96434595) q[3];
sx q[3];
rz(-2.3026005) q[3];
sx q[3];
rz(2.8540924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
