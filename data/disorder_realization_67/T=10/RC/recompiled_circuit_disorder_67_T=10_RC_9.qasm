OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(-1.0597205) q[0];
sx q[0];
rz(0.73097316) q[0];
rz(1.641474) q[1];
sx q[1];
rz(5.2483622) q[1];
sx q[1];
rz(11.622826) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1424926) q[0];
sx q[0];
rz(-1.5038135) q[0];
sx q[0];
rz(1.5357369) q[0];
rz(-1.4859096) q[2];
sx q[2];
rz(-0.92157084) q[2];
sx q[2];
rz(-2.6390586) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.700625) q[1];
sx q[1];
rz(-1.1807627) q[1];
sx q[1];
rz(-2.7729211) q[1];
rz(0.30480095) q[3];
sx q[3];
rz(-1.7497352) q[3];
sx q[3];
rz(-1.6182871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8865108) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(1.8908267) q[2];
rz(1.7154153) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(-0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.132906) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(0.59894484) q[0];
rz(-1.8006181) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(0.96639955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67993977) q[0];
sx q[0];
rz(-2.0113809) q[0];
sx q[0];
rz(2.5011714) q[0];
rz(-0.62218372) q[2];
sx q[2];
rz(-1.6112279) q[2];
sx q[2];
rz(-1.4036904) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2981373) q[1];
sx q[1];
rz(-1.6846488) q[1];
sx q[1];
rz(1.650412) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3923799) q[3];
sx q[3];
rz(-1.2470761) q[3];
sx q[3];
rz(0.92326984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0559343) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-0.30109626) q[2];
rz(-1.1931233) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(1.1863856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(3.0369204) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(-1.954129) q[0];
rz(1.2359515) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.8240066) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353969) q[0];
sx q[0];
rz(-2.4164696) q[0];
sx q[0];
rz(-2.8185185) q[0];
rz(0.74430978) q[2];
sx q[2];
rz(-0.28713206) q[2];
sx q[2];
rz(-0.14936514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8220362) q[1];
sx q[1];
rz(-1.6720547) q[1];
sx q[1];
rz(-0.42145573) q[1];
rz(-pi) q[2];
rz(2.0882323) q[3];
sx q[3];
rz(-1.468717) q[3];
sx q[3];
rz(-1.7774338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(2.9349566) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(-2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77465039) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(-2.2553717) q[0];
rz(2.1318502) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(-1.9151691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7851631) q[0];
sx q[0];
rz(-2.4117081) q[0];
sx q[0];
rz(-1.6701783) q[0];
x q[1];
rz(-2.6817276) q[2];
sx q[2];
rz(-0.93732873) q[2];
sx q[2];
rz(-1.8749274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7734079) q[1];
sx q[1];
rz(-0.28973026) q[1];
sx q[1];
rz(0.75399953) q[1];
rz(0.63663441) q[3];
sx q[3];
rz(-2.5599179) q[3];
sx q[3];
rz(-3.0326774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6440789) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(0.99299661) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-0.48373568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(1.1859878) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-2.5591154) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7459481) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(-0.086918513) q[0];
rz(0.14898665) q[2];
sx q[2];
rz(-2.2431231) q[2];
sx q[2];
rz(2.7186) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.40618784) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(0.71449844) q[1];
x q[2];
rz(-1.6378239) q[3];
sx q[3];
rz(-1.1653295) q[3];
sx q[3];
rz(-2.1342579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(-0.081710903) q[2];
rz(-0.47406667) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(-1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.896647) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(-0.0078049302) q[0];
rz(1.4004978) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(2.0369464) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1305599) q[0];
sx q[0];
rz(-1.7337807) q[0];
sx q[0];
rz(-2.4736604) q[0];
x q[1];
rz(2.5014624) q[2];
sx q[2];
rz(-1.374561) q[2];
sx q[2];
rz(-0.093984691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2075344) q[1];
sx q[1];
rz(-0.96734069) q[1];
sx q[1];
rz(2.250309) q[1];
rz(2.036318) q[3];
sx q[3];
rz(-1.421531) q[3];
sx q[3];
rz(-2.535459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(0.55523038) q[2];
rz(-2.9688719) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.5482607) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5884488) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(3.0798262) q[0];
rz(-0.24208367) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(-1.1118836) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2568946) q[0];
sx q[0];
rz(-1.9848794) q[0];
sx q[0];
rz(2.3143682) q[0];
rz(-pi) q[1];
x q[1];
rz(1.493181) q[2];
sx q[2];
rz(-1.5570119) q[2];
sx q[2];
rz(0.26253653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.249914) q[1];
sx q[1];
rz(-0.87712446) q[1];
sx q[1];
rz(1.9338884) q[1];
rz(-pi) q[2];
rz(-1.816733) q[3];
sx q[3];
rz(-1.2517559) q[3];
sx q[3];
rz(-2.6815573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3322488) q[2];
sx q[2];
rz(-0.76433864) q[2];
sx q[2];
rz(0.81364441) q[2];
rz(-1.404473) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(-2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5161045) q[0];
sx q[0];
rz(-1.3502716) q[0];
sx q[0];
rz(-1.7161436) q[0];
rz(1.6199934) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(2.5040748) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.92256) q[0];
sx q[0];
rz(-1.0244644) q[0];
sx q[0];
rz(2.4541897) q[0];
x q[1];
rz(1.2136739) q[2];
sx q[2];
rz(-0.87840688) q[2];
sx q[2];
rz(1.1536319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1453104) q[1];
sx q[1];
rz(-1.9687555) q[1];
sx q[1];
rz(0.63509649) q[1];
rz(-pi) q[2];
rz(1.4333126) q[3];
sx q[3];
rz(-2.1372037) q[3];
sx q[3];
rz(0.24368225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1311538) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(-1.1361702) q[2];
rz(-1.4853959) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(1.9320528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5159601) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(-0.82538429) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(-2.1264145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4207941) q[0];
sx q[0];
rz(-1.7013229) q[0];
sx q[0];
rz(-0.41850787) q[0];
x q[1];
rz(-1.7245618) q[2];
sx q[2];
rz(-0.92369881) q[2];
sx q[2];
rz(-0.53714067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79267348) q[1];
sx q[1];
rz(-1.9209849) q[1];
sx q[1];
rz(-2.7677571) q[1];
x q[2];
rz(0.34992643) q[3];
sx q[3];
rz(-2.3438128) q[3];
sx q[3];
rz(1.6328904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.72145808) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(0.96735111) q[2];
rz(-1.5445276) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6185146) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(-1.7369695) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4182189) q[0];
sx q[0];
rz(-1.2105816) q[0];
sx q[0];
rz(-0.13803137) q[0];
rz(-pi) q[1];
rz(1.6925473) q[2];
sx q[2];
rz(-1.2687614) q[2];
sx q[2];
rz(-1.2388602) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2056634) q[1];
sx q[1];
rz(-2.231039) q[1];
sx q[1];
rz(-0.6026938) q[1];
rz(2.9910562) q[3];
sx q[3];
rz(-1.5794465) q[3];
sx q[3];
rz(2.0899525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.848032) q[2];
sx q[2];
rz(-2.3876987) q[2];
sx q[2];
rz(-0.30612293) q[2];
rz(2.7434769) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33481471) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.1595935) q[1];
sx q[1];
rz(-1.0212785) q[1];
sx q[1];
rz(-2.8550128) q[1];
rz(0.26328662) q[2];
sx q[2];
rz(-1.2821715) q[2];
sx q[2];
rz(0.53496219) q[2];
rz(1.0767827) q[3];
sx q[3];
rz(-1.8643338) q[3];
sx q[3];
rz(-2.1716933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
