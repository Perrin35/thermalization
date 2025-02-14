OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8706239) q[0];
sx q[0];
rz(3.6977036) q[0];
sx q[0];
rz(7.2365427) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(1.035773) q[1];
sx q[1];
rz(8.4761578) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9069288) q[0];
sx q[0];
rz(-1.0415823) q[0];
sx q[0];
rz(2.0509023) q[0];
rz(-pi) q[1];
rz(-0.23748246) q[2];
sx q[2];
rz(-1.5680702) q[2];
sx q[2];
rz(-0.021955519) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.030101209) q[1];
sx q[1];
rz(-1.1426569) q[1];
sx q[1];
rz(0.74049048) q[1];
rz(2.1236739) q[3];
sx q[3];
rz(-1.8686668) q[3];
sx q[3];
rz(1.2978993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76401508) q[2];
sx q[2];
rz(-3.0674051) q[2];
sx q[2];
rz(0.78262502) q[2];
rz(-3.022656) q[3];
sx q[3];
rz(-2.1075893) q[3];
sx q[3];
rz(-2.9940166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19621944) q[0];
sx q[0];
rz(-1.1213028) q[0];
sx q[0];
rz(-0.25783208) q[0];
rz(0.091015426) q[1];
sx q[1];
rz(-2.0423404) q[1];
sx q[1];
rz(-1.6450504) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2473619) q[0];
sx q[0];
rz(-1.6391488) q[0];
sx q[0];
rz(-1.0175045) q[0];
rz(-1.2256669) q[2];
sx q[2];
rz(-1.218443) q[2];
sx q[2];
rz(2.8690668) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0219258) q[1];
sx q[1];
rz(-2.7283784) q[1];
sx q[1];
rz(1.4959123) q[1];
rz(-pi) q[2];
rz(1.9807545) q[3];
sx q[3];
rz(-2.1968968) q[3];
sx q[3];
rz(-0.61625851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0195134) q[2];
sx q[2];
rz(-1.2747108) q[2];
sx q[2];
rz(2.7369734) q[2];
rz(2.1520065) q[3];
sx q[3];
rz(-1.3000969) q[3];
sx q[3];
rz(2.5185481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0290381) q[0];
sx q[0];
rz(-0.14415388) q[0];
sx q[0];
rz(-2.3032904) q[0];
rz(2.2932032) q[1];
sx q[1];
rz(-1.9219857) q[1];
sx q[1];
rz(-0.99260509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557195) q[0];
sx q[0];
rz(-0.8536754) q[0];
sx q[0];
rz(-1.3846057) q[0];
rz(-2.8072629) q[2];
sx q[2];
rz(-1.624776) q[2];
sx q[2];
rz(0.24654972) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.015945399) q[1];
sx q[1];
rz(-1.1635259) q[1];
sx q[1];
rz(-1.8512878) q[1];
rz(-pi) q[2];
rz(-1.6818524) q[3];
sx q[3];
rz(-1.8594311) q[3];
sx q[3];
rz(1.9185324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8695716) q[2];
sx q[2];
rz(-1.1676936) q[2];
sx q[2];
rz(-1.1248379) q[2];
rz(-2.520842) q[3];
sx q[3];
rz(-0.97095942) q[3];
sx q[3];
rz(3.0264405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426185) q[0];
sx q[0];
rz(-1.9816575) q[0];
sx q[0];
rz(2.9677891) q[0];
rz(2.4558892) q[1];
sx q[1];
rz(-1.4861636) q[1];
sx q[1];
rz(-0.83820835) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3371171) q[0];
sx q[0];
rz(-2.4675779) q[0];
sx q[0];
rz(0.78916691) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.051866) q[2];
sx q[2];
rz(-1.0609259) q[2];
sx q[2];
rz(0.0048268371) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9058075) q[1];
sx q[1];
rz(-0.90079325) q[1];
sx q[1];
rz(1.7520105) q[1];
rz(-pi) q[2];
rz(1.6082786) q[3];
sx q[3];
rz(-1.1082543) q[3];
sx q[3];
rz(-1.2087865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59914261) q[2];
sx q[2];
rz(-1.6833545) q[2];
sx q[2];
rz(-2.7899) q[2];
rz(1.1951949) q[3];
sx q[3];
rz(-1.1870793) q[3];
sx q[3];
rz(-0.024141969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.675932) q[0];
sx q[0];
rz(-0.52023482) q[0];
sx q[0];
rz(2.7886673) q[0];
rz(2.832761) q[1];
sx q[1];
rz(-1.0363204) q[1];
sx q[1];
rz(0.028506361) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7362979) q[0];
sx q[0];
rz(-2.5319244) q[0];
sx q[0];
rz(1.8788615) q[0];
x q[1];
rz(-1.9602922) q[2];
sx q[2];
rz(-1.2151698) q[2];
sx q[2];
rz(2.2782922) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5529768) q[1];
sx q[1];
rz(-1.7874831) q[1];
sx q[1];
rz(-0.79766794) q[1];
rz(-1.0794382) q[3];
sx q[3];
rz(-2.6565928) q[3];
sx q[3];
rz(2.5122143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21294022) q[2];
sx q[2];
rz(-3.0781367) q[2];
sx q[2];
rz(-2.1707936) q[2];
rz(-1.0468696) q[3];
sx q[3];
rz(-1.6887083) q[3];
sx q[3];
rz(0.48986062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30188072) q[0];
sx q[0];
rz(-1.3588926) q[0];
sx q[0];
rz(-3.0506328) q[0];
rz(1.2184527) q[1];
sx q[1];
rz(-2.8107042) q[1];
sx q[1];
rz(-2.5023696) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7122593) q[0];
sx q[0];
rz(-0.3540701) q[0];
sx q[0];
rz(0.39850111) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0193954) q[2];
sx q[2];
rz(-1.4336042) q[2];
sx q[2];
rz(0.73681632) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7057695) q[1];
sx q[1];
rz(-2.2652049) q[1];
sx q[1];
rz(-2.4468876) q[1];
rz(-pi) q[2];
rz(1.8973783) q[3];
sx q[3];
rz(-1.5726798) q[3];
sx q[3];
rz(-0.43549805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0794534) q[2];
sx q[2];
rz(-0.97525758) q[2];
sx q[2];
rz(-0.47387588) q[2];
rz(-1.8114932) q[3];
sx q[3];
rz(-2.2606943) q[3];
sx q[3];
rz(0.37548319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5652931) q[0];
sx q[0];
rz(-1.0813035) q[0];
sx q[0];
rz(-0.22258776) q[0];
rz(2.0780156) q[1];
sx q[1];
rz(-1.56366) q[1];
sx q[1];
rz(1.9482013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8603966) q[0];
sx q[0];
rz(-1.4544857) q[0];
sx q[0];
rz(0.75006811) q[0];
x q[1];
rz(1.3435828) q[2];
sx q[2];
rz(-1.6889204) q[2];
sx q[2];
rz(0.83641499) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30380043) q[1];
sx q[1];
rz(-1.1322316) q[1];
sx q[1];
rz(1.720251) q[1];
x q[2];
rz(2.9783747) q[3];
sx q[3];
rz(-1.8344518) q[3];
sx q[3];
rz(1.4782983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0003164) q[2];
sx q[2];
rz(-2.2681984) q[2];
sx q[2];
rz(-0.6401965) q[2];
rz(-0.021942465) q[3];
sx q[3];
rz(-1.7317737) q[3];
sx q[3];
rz(0.35023165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6496277) q[0];
sx q[0];
rz(-1.7438629) q[0];
sx q[0];
rz(2.6212027) q[0];
rz(-0.87977663) q[1];
sx q[1];
rz(-1.2326515) q[1];
sx q[1];
rz(-2.2487776) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8917903) q[0];
sx q[0];
rz(-3.043104) q[0];
sx q[0];
rz(-0.47827025) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73142902) q[2];
sx q[2];
rz(-0.9584223) q[2];
sx q[2];
rz(-2.2867672) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8930606) q[1];
sx q[1];
rz(-1.8527781) q[1];
sx q[1];
rz(-2.3470013) q[1];
rz(2.625301) q[3];
sx q[3];
rz(-1.9441351) q[3];
sx q[3];
rz(3.1331799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2616547) q[2];
sx q[2];
rz(-2.0377906) q[2];
sx q[2];
rz(3.0547764) q[2];
rz(-0.9451198) q[3];
sx q[3];
rz(-0.34543959) q[3];
sx q[3];
rz(-1.1300348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18411186) q[0];
sx q[0];
rz(-0.090395398) q[0];
sx q[0];
rz(-3.0775253) q[0];
rz(-2.0059026) q[1];
sx q[1];
rz(-1.7013197) q[1];
sx q[1];
rz(0.51220977) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50822608) q[0];
sx q[0];
rz(-2.6056943) q[0];
sx q[0];
rz(-0.6750489) q[0];
rz(-1.1697328) q[2];
sx q[2];
rz(-2.3980015) q[2];
sx q[2];
rz(2.0498073) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87742245) q[1];
sx q[1];
rz(-1.3119196) q[1];
sx q[1];
rz(1.4822465) q[1];
rz(-pi) q[2];
rz(-0.036831696) q[3];
sx q[3];
rz(-2.7324711) q[3];
sx q[3];
rz(-2.6927519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3247437) q[2];
sx q[2];
rz(-0.73255676) q[2];
sx q[2];
rz(2.0762439) q[2];
rz(-1.4893074) q[3];
sx q[3];
rz(-0.95732147) q[3];
sx q[3];
rz(-3.0491507) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5866933) q[0];
sx q[0];
rz(-0.41189343) q[0];
sx q[0];
rz(-0.33083415) q[0];
rz(1.2384442) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(0.27473658) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82652826) q[0];
sx q[0];
rz(-1.9103658) q[0];
sx q[0];
rz(1.1134321) q[0];
x q[1];
rz(3.0255453) q[2];
sx q[2];
rz(-0.2751285) q[2];
sx q[2];
rz(0.76062084) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3420136) q[1];
sx q[1];
rz(-1.430961) q[1];
sx q[1];
rz(2.0796156) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.072321691) q[3];
sx q[3];
rz(-0.9843217) q[3];
sx q[3];
rz(-2.9972635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1642509) q[2];
sx q[2];
rz(-2.3873316) q[2];
sx q[2];
rz(1.784262) q[2];
rz(0.50061289) q[3];
sx q[3];
rz(-1.5631792) q[3];
sx q[3];
rz(-1.64465) q[3];
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
rz(1.5464583) q[0];
sx q[0];
rz(-1.5820137) q[0];
sx q[0];
rz(-0.59376846) q[0];
rz(3.0733227) q[1];
sx q[1];
rz(-1.2208114) q[1];
sx q[1];
rz(0.023718871) q[1];
rz(-2.7958128) q[2];
sx q[2];
rz(-1.4152534) q[2];
sx q[2];
rz(-1.5243668) q[2];
rz(0.17791893) q[3];
sx q[3];
rz(-1.3861227) q[3];
sx q[3];
rz(1.9021195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
