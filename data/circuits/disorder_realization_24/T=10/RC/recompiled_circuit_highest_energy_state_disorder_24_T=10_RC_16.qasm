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
rz(-0.10676323) q[0];
sx q[0];
rz(-1.8101298) q[0];
sx q[0];
rz(-1.3504299) q[0];
rz(0.30836937) q[1];
sx q[1];
rz(-0.24873939) q[1];
sx q[1];
rz(2.1623936) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0672507) q[0];
sx q[0];
rz(-0.76267159) q[0];
sx q[0];
rz(3.0347155) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1264231) q[2];
sx q[2];
rz(-0.70894402) q[2];
sx q[2];
rz(2.5042748) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82047082) q[1];
sx q[1];
rz(-2.781457) q[1];
sx q[1];
rz(2.9481357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37110801) q[3];
sx q[3];
rz(-1.4083997) q[3];
sx q[3];
rz(-1.1134256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2778492) q[2];
sx q[2];
rz(-1.2110445) q[2];
sx q[2];
rz(-0.63300526) q[2];
rz(0.64166075) q[3];
sx q[3];
rz(-2.6810665) q[3];
sx q[3];
rz(2.3708926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4577456) q[0];
sx q[0];
rz(-0.96413833) q[0];
sx q[0];
rz(0.82474166) q[0];
rz(-1.224158) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(-2.8214084) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.510957) q[0];
sx q[0];
rz(-0.63586006) q[0];
sx q[0];
rz(-2.4102734) q[0];
rz(-pi) q[1];
rz(1.0616706) q[2];
sx q[2];
rz(-1.3855644) q[2];
sx q[2];
rz(2.9502466) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0517438) q[1];
sx q[1];
rz(-1.4898572) q[1];
sx q[1];
rz(-1.593394) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3125962) q[3];
sx q[3];
rz(-2.9123983) q[3];
sx q[3];
rz(2.3865139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1227293) q[2];
sx q[2];
rz(-1.9247232) q[2];
sx q[2];
rz(1.2790206) q[2];
rz(-0.095712885) q[3];
sx q[3];
rz(-2.0666104) q[3];
sx q[3];
rz(1.8224645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9886446) q[0];
sx q[0];
rz(-2.4929292) q[0];
sx q[0];
rz(1.0308107) q[0];
rz(-0.55054322) q[1];
sx q[1];
rz(-0.85854733) q[1];
sx q[1];
rz(-1.2446838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7787832) q[0];
sx q[0];
rz(-0.16805695) q[0];
sx q[0];
rz(2.8166978) q[0];
x q[1];
rz(1.8049525) q[2];
sx q[2];
rz(-2.5114254) q[2];
sx q[2];
rz(2.5182398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0708187) q[1];
sx q[1];
rz(-0.70291513) q[1];
sx q[1];
rz(1.2513112) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9774417) q[3];
sx q[3];
rz(-2.3465524) q[3];
sx q[3];
rz(2.1926978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9024258) q[2];
sx q[2];
rz(-2.5776165) q[2];
sx q[2];
rz(1.172056) q[2];
rz(-2.3160882) q[3];
sx q[3];
rz(-1.2647311) q[3];
sx q[3];
rz(0.66506213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.0498407) q[0];
sx q[0];
rz(-1.6850152) q[0];
sx q[0];
rz(0.41393429) q[0];
rz(-0.44231689) q[1];
sx q[1];
rz(-1.0041753) q[1];
sx q[1];
rz(0.54534674) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6028311) q[0];
sx q[0];
rz(-1.6206546) q[0];
sx q[0];
rz(-2.8661348) q[0];
rz(-0.77540437) q[2];
sx q[2];
rz(-1.0808965) q[2];
sx q[2];
rz(-2.6012914) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8607747) q[1];
sx q[1];
rz(-2.0866261) q[1];
sx q[1];
rz(2.8223553) q[1];
x q[2];
rz(-2.993587) q[3];
sx q[3];
rz(-2.0750351) q[3];
sx q[3];
rz(-1.809169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44117323) q[2];
sx q[2];
rz(-1.8808695) q[2];
sx q[2];
rz(0.1612266) q[2];
rz(1.553933) q[3];
sx q[3];
rz(-2.5416538) q[3];
sx q[3];
rz(2.5509295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081838354) q[0];
sx q[0];
rz(-1.8738382) q[0];
sx q[0];
rz(-0.72342122) q[0];
rz(1.3149892) q[1];
sx q[1];
rz(-1.2774717) q[1];
sx q[1];
rz(1.0378729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4525131) q[0];
sx q[0];
rz(-0.85001365) q[0];
sx q[0];
rz(-1.6759273) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5332301) q[2];
sx q[2];
rz(-0.82348292) q[2];
sx q[2];
rz(-2.0840816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8602058) q[1];
sx q[1];
rz(-1.1347927) q[1];
sx q[1];
rz(-1.866775) q[1];
rz(-2.9014189) q[3];
sx q[3];
rz(-0.55910149) q[3];
sx q[3];
rz(0.19776519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7448685) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(0.57207668) q[2];
rz(0.62503302) q[3];
sx q[3];
rz(-2.1376938) q[3];
sx q[3];
rz(1.9151789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201651) q[0];
sx q[0];
rz(-0.033997424) q[0];
sx q[0];
rz(-2.6183364) q[0];
rz(1.8190067) q[1];
sx q[1];
rz(-1.6736284) q[1];
sx q[1];
rz(-0.62017131) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70411982) q[0];
sx q[0];
rz(-1.0123555) q[0];
sx q[0];
rz(1.062514) q[0];
rz(-pi) q[1];
rz(2.6729726) q[2];
sx q[2];
rz(-1.5053362) q[2];
sx q[2];
rz(0.58131389) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.090724548) q[1];
sx q[1];
rz(-1.4830989) q[1];
sx q[1];
rz(-0.94322738) q[1];
x q[2];
rz(-2.9240588) q[3];
sx q[3];
rz(-1.399588) q[3];
sx q[3];
rz(-1.8678766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4124477) q[2];
sx q[2];
rz(-0.77249384) q[2];
sx q[2];
rz(1.8577417) q[2];
rz(-2.1733952) q[3];
sx q[3];
rz(-1.3788393) q[3];
sx q[3];
rz(1.5865954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17484434) q[0];
sx q[0];
rz(-2.0059678) q[0];
sx q[0];
rz(1.9298166) q[0];
rz(-2.1719596) q[1];
sx q[1];
rz(-1.8200834) q[1];
sx q[1];
rz(-1.8671573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6912708) q[0];
sx q[0];
rz(-0.88435882) q[0];
sx q[0];
rz(-0.19569389) q[0];
rz(-0.78584255) q[2];
sx q[2];
rz(-1.5034901) q[2];
sx q[2];
rz(2.6493361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3434747) q[1];
sx q[1];
rz(-0.96303446) q[1];
sx q[1];
rz(-2.7704984) q[1];
rz(-pi) q[2];
rz(-0.59341615) q[3];
sx q[3];
rz(-0.53916603) q[3];
sx q[3];
rz(-0.19684686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4659287) q[2];
sx q[2];
rz(-1.482654) q[2];
sx q[2];
rz(2.8181804) q[2];
rz(3.1392426) q[3];
sx q[3];
rz(-1.6833143) q[3];
sx q[3];
rz(0.6215483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49331409) q[0];
sx q[0];
rz(-1.6596153) q[0];
sx q[0];
rz(1.3508654) q[0];
rz(1.3510652) q[1];
sx q[1];
rz(-1.8565145) q[1];
sx q[1];
rz(-1.8538808) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13437102) q[0];
sx q[0];
rz(-0.30728455) q[0];
sx q[0];
rz(1.968712) q[0];
rz(0.28952629) q[2];
sx q[2];
rz(-2.2643201) q[2];
sx q[2];
rz(0.62563459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2877174) q[1];
sx q[1];
rz(-2.7734904) q[1];
sx q[1];
rz(2.6516857) q[1];
rz(-pi) q[2];
rz(1.5519556) q[3];
sx q[3];
rz(-1.2219489) q[3];
sx q[3];
rz(2.2486389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52435851) q[2];
sx q[2];
rz(-1.989216) q[2];
sx q[2];
rz(-2.1481245) q[2];
rz(-1.0348381) q[3];
sx q[3];
rz(-2.5351758) q[3];
sx q[3];
rz(0.16916999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.33191037) q[0];
sx q[0];
rz(-1.1746291) q[0];
sx q[0];
rz(1.5274973) q[0];
rz(2.2612259) q[1];
sx q[1];
rz(-2.7208734) q[1];
sx q[1];
rz(-0.31164935) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0535705) q[0];
sx q[0];
rz(-2.7854475) q[0];
sx q[0];
rz(2.1689264) q[0];
rz(-pi) q[1];
rz(1.3394322) q[2];
sx q[2];
rz(-1.5548493) q[2];
sx q[2];
rz(2.0248264) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2276409) q[1];
sx q[1];
rz(-2.3415509) q[1];
sx q[1];
rz(2.6197919) q[1];
rz(-pi) q[2];
rz(-0.33297379) q[3];
sx q[3];
rz(-2.1978488) q[3];
sx q[3];
rz(-2.0505333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46939048) q[2];
sx q[2];
rz(-2.2180874) q[2];
sx q[2];
rz(-1.3908609) q[2];
rz(1.5270799) q[3];
sx q[3];
rz(-2.0938087) q[3];
sx q[3];
rz(2.0670149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63956082) q[0];
sx q[0];
rz(-1.9968411) q[0];
sx q[0];
rz(2.5290053) q[0];
rz(2.1514429) q[1];
sx q[1];
rz(-1.1724816) q[1];
sx q[1];
rz(1.490907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4717141) q[0];
sx q[0];
rz(-1.151514) q[0];
sx q[0];
rz(-0.74434481) q[0];
rz(-pi) q[1];
rz(-2.2305626) q[2];
sx q[2];
rz(-1.7740031) q[2];
sx q[2];
rz(-0.045762941) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4348169) q[1];
sx q[1];
rz(-1.5770327) q[1];
sx q[1];
rz(1.4093136) q[1];
rz(-2.8464523) q[3];
sx q[3];
rz(-0.71986976) q[3];
sx q[3];
rz(0.9269549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.426173) q[2];
sx q[2];
rz(-2.3121068) q[2];
sx q[2];
rz(3.1066011) q[2];
rz(1.4404826) q[3];
sx q[3];
rz(-2.248843) q[3];
sx q[3];
rz(2.6006202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69915019) q[0];
sx q[0];
rz(-2.2461666) q[0];
sx q[0];
rz(2.9641892) q[0];
rz(1.9433446) q[1];
sx q[1];
rz(-1.6234963) q[1];
sx q[1];
rz(-2.1877098) q[1];
rz(-2.4631029) q[2];
sx q[2];
rz(-2.0526921) q[2];
sx q[2];
rz(-1.0865281) q[2];
rz(1.1097601) q[3];
sx q[3];
rz(-2.2363935) q[3];
sx q[3];
rz(1.4105894) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
