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
rz(0.97419089) q[0];
sx q[0];
rz(-2.2926957) q[0];
sx q[0];
rz(1.726024) q[0];
rz(-0.8610791) q[1];
sx q[1];
rz(-1.6281035) q[1];
sx q[1];
rz(-2.3754062) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7216033) q[0];
sx q[0];
rz(-0.074062183) q[0];
sx q[0];
rz(-0.5038528) q[0];
rz(1.8433601) q[2];
sx q[2];
rz(-0.12530357) q[2];
sx q[2];
rz(-1.7150777) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76152245) q[1];
sx q[1];
rz(-2.9381611) q[1];
sx q[1];
rz(3.1043209) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9835408) q[3];
sx q[3];
rz(-1.6696602) q[3];
sx q[3];
rz(-0.017798558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2506931) q[2];
sx q[2];
rz(-0.74960274) q[2];
sx q[2];
rz(-3.0787943) q[2];
rz(-1.4248258) q[3];
sx q[3];
rz(-1.7664322) q[3];
sx q[3];
rz(-2.4019901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-3.0733136) q[0];
sx q[0];
rz(-2.8977019) q[0];
sx q[0];
rz(-2.0252315) q[0];
rz(-1.5694537) q[1];
sx q[1];
rz(-0.29466033) q[1];
sx q[1];
rz(2.0223298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9452764) q[0];
sx q[0];
rz(-2.2782562) q[0];
sx q[0];
rz(-1.5718801) q[0];
rz(-pi) q[1];
rz(0.92932911) q[2];
sx q[2];
rz(-0.57999883) q[2];
sx q[2];
rz(1.3652238) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83554448) q[1];
sx q[1];
rz(-1.5054963) q[1];
sx q[1];
rz(0.52564486) q[1];
rz(0.6898361) q[3];
sx q[3];
rz(-0.82131006) q[3];
sx q[3];
rz(0.46508712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.08673) q[2];
sx q[2];
rz(-2.3708998) q[2];
sx q[2];
rz(1.1718303) q[2];
rz(-2.150676) q[3];
sx q[3];
rz(-2.1358229) q[3];
sx q[3];
rz(-1.2539554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1427926) q[0];
sx q[0];
rz(-0.87738377) q[0];
sx q[0];
rz(-3.0679833) q[0];
rz(-0.40120801) q[1];
sx q[1];
rz(-0.3708655) q[1];
sx q[1];
rz(2.4295095) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129666) q[0];
sx q[0];
rz(-3.0168128) q[0];
sx q[0];
rz(-1.9931884) q[0];
rz(3.1396578) q[2];
sx q[2];
rz(-2.4855485) q[2];
sx q[2];
rz(-1.4081788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69171762) q[1];
sx q[1];
rz(-0.32539168) q[1];
sx q[1];
rz(-0.34145932) q[1];
rz(2.1041342) q[3];
sx q[3];
rz(-1.520562) q[3];
sx q[3];
rz(0.929757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3965343) q[2];
sx q[2];
rz(-1.859954) q[2];
sx q[2];
rz(-0.38953951) q[2];
rz(-0.39858308) q[3];
sx q[3];
rz(-1.7321209) q[3];
sx q[3];
rz(1.4220953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013537708) q[0];
sx q[0];
rz(-0.062352926) q[0];
sx q[0];
rz(-2.7259977) q[0];
rz(1.8889113) q[1];
sx q[1];
rz(-1.1377734) q[1];
sx q[1];
rz(-2.2223053) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.141826) q[0];
sx q[0];
rz(-1.5298889) q[0];
sx q[0];
rz(-1.7342939) q[0];
rz(-pi) q[1];
x q[1];
rz(0.083468632) q[2];
sx q[2];
rz(-2.8776904) q[2];
sx q[2];
rz(-0.47191274) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34621908) q[1];
sx q[1];
rz(-0.88118964) q[1];
sx q[1];
rz(-2.3350689) q[1];
rz(-pi) q[2];
rz(-2.0218798) q[3];
sx q[3];
rz(-2.2871784) q[3];
sx q[3];
rz(0.35600397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2933942) q[2];
sx q[2];
rz(-1.5092809) q[2];
sx q[2];
rz(2.1634114) q[2];
rz(0.50655043) q[3];
sx q[3];
rz(-1.6440441) q[3];
sx q[3];
rz(-1.0165366) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070234805) q[0];
sx q[0];
rz(-2.4429584) q[0];
sx q[0];
rz(-0.51682669) q[0];
rz(1.0454073) q[1];
sx q[1];
rz(-2.1357336) q[1];
sx q[1];
rz(2.7396835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088083) q[0];
sx q[0];
rz(-2.0044745) q[0];
sx q[0];
rz(0.78701235) q[0];
x q[1];
rz(1.044173) q[2];
sx q[2];
rz(-0.32794558) q[2];
sx q[2];
rz(1.8484268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0037108) q[1];
sx q[1];
rz(-1.5750843) q[1];
sx q[1];
rz(1.3715308) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7414209) q[3];
sx q[3];
rz(-0.62230357) q[3];
sx q[3];
rz(-2.4476652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4154633) q[2];
sx q[2];
rz(-1.9620506) q[2];
sx q[2];
rz(-2.9187875) q[2];
rz(-2.2807617) q[3];
sx q[3];
rz(-1.5102883) q[3];
sx q[3];
rz(1.2739325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.05805) q[0];
sx q[0];
rz(-1.2419751) q[0];
sx q[0];
rz(-0.85839957) q[0];
rz(0.79943132) q[1];
sx q[1];
rz(-2.1245427) q[1];
sx q[1];
rz(-1.3851091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7290447) q[0];
sx q[0];
rz(-2.8083889) q[0];
sx q[0];
rz(0.97107519) q[0];
x q[1];
rz(-3.0869003) q[2];
sx q[2];
rz(-1.2451951) q[2];
sx q[2];
rz(0.97132746) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4614758) q[1];
sx q[1];
rz(-1.6915186) q[1];
sx q[1];
rz(0.6013479) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0603745) q[3];
sx q[3];
rz(-1.7822232) q[3];
sx q[3];
rz(1.8544153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4110306) q[2];
sx q[2];
rz(-0.73341122) q[2];
sx q[2];
rz(-0.15156558) q[2];
rz(1.9568806) q[3];
sx q[3];
rz(-1.1368753) q[3];
sx q[3];
rz(2.241551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93722349) q[0];
sx q[0];
rz(-1.7974412) q[0];
sx q[0];
rz(-2.3580661) q[0];
rz(-0.20768684) q[1];
sx q[1];
rz(-2.2238104) q[1];
sx q[1];
rz(-1.2088561) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3850526) q[0];
sx q[0];
rz(-1.970743) q[0];
sx q[0];
rz(-0.51838309) q[0];
rz(-pi) q[1];
rz(-0.98548205) q[2];
sx q[2];
rz(-0.83067465) q[2];
sx q[2];
rz(-0.58500803) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5294339) q[1];
sx q[1];
rz(-2.4116257) q[1];
sx q[1];
rz(1.5117701) q[1];
rz(-pi) q[2];
rz(1.1999628) q[3];
sx q[3];
rz(-1.7177903) q[3];
sx q[3];
rz(-1.539409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8565159) q[2];
sx q[2];
rz(-0.74399844) q[2];
sx q[2];
rz(-0.12459717) q[2];
rz(-1.8018657) q[3];
sx q[3];
rz(-2.7281269) q[3];
sx q[3];
rz(-1.193803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1272005) q[0];
sx q[0];
rz(-1.7008282) q[0];
sx q[0];
rz(0.097323962) q[0];
rz(2.4864181) q[1];
sx q[1];
rz(-0.56071463) q[1];
sx q[1];
rz(-2.4526144) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5118588) q[0];
sx q[0];
rz(-1.0847158) q[0];
sx q[0];
rz(0.084045579) q[0];
x q[1];
rz(-1.9227486) q[2];
sx q[2];
rz(-1.5492288) q[2];
sx q[2];
rz(-1.5605614) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6072907) q[1];
sx q[1];
rz(-1.4431776) q[1];
sx q[1];
rz(-2.4625596) q[1];
x q[2];
rz(0.19676713) q[3];
sx q[3];
rz(-0.69878687) q[3];
sx q[3];
rz(-3.1270032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.24136209) q[2];
sx q[2];
rz(-2.1817544) q[2];
sx q[2];
rz(1.0660561) q[2];
rz(-2.7212932) q[3];
sx q[3];
rz(-1.2371141) q[3];
sx q[3];
rz(1.458781) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083753839) q[0];
sx q[0];
rz(-2.3729615) q[0];
sx q[0];
rz(0.88430697) q[0];
rz(3.0506483) q[1];
sx q[1];
rz(-0.93946409) q[1];
sx q[1];
rz(-1.9627862) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55806736) q[0];
sx q[0];
rz(-1.5060165) q[0];
sx q[0];
rz(-2.9811048) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5680204) q[2];
sx q[2];
rz(-0.66589543) q[2];
sx q[2];
rz(-1.8393054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2388666) q[1];
sx q[1];
rz(-1.7416493) q[1];
sx q[1];
rz(-2.577223) q[1];
rz(-pi) q[2];
rz(-1.6903773) q[3];
sx q[3];
rz(-1.371583) q[3];
sx q[3];
rz(1.8389033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74943298) q[2];
sx q[2];
rz(-1.1460816) q[2];
sx q[2];
rz(2.667099) q[2];
rz(0.042304603) q[3];
sx q[3];
rz(-1.6241122) q[3];
sx q[3];
rz(-1.8208985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45847356) q[0];
sx q[0];
rz(-2.0404158) q[0];
sx q[0];
rz(1.3018357) q[0];
rz(0.63109541) q[1];
sx q[1];
rz(-2.0986235) q[1];
sx q[1];
rz(1.016681) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4563718) q[0];
sx q[0];
rz(-2.4349762) q[0];
sx q[0];
rz(-2.0232692) q[0];
rz(-pi) q[1];
rz(-1.6696681) q[2];
sx q[2];
rz(-1.6244096) q[2];
sx q[2];
rz(0.19744548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8606372) q[1];
sx q[1];
rz(-1.5796575) q[1];
sx q[1];
rz(-0.85207664) q[1];
rz(-pi) q[2];
rz(-2.2325114) q[3];
sx q[3];
rz(-0.98041816) q[3];
sx q[3];
rz(0.31766858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8299228) q[2];
sx q[2];
rz(-0.72511464) q[2];
sx q[2];
rz(-2.1666918) q[2];
rz(1.4026862) q[3];
sx q[3];
rz(-0.97167531) q[3];
sx q[3];
rz(1.5506802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0065895157) q[0];
sx q[0];
rz(-2.4000744) q[0];
sx q[0];
rz(0.97493521) q[0];
rz(-1.4047752) q[1];
sx q[1];
rz(-1.732323) q[1];
sx q[1];
rz(2.7224532) q[1];
rz(-0.51437893) q[2];
sx q[2];
rz(-1.0399401) q[2];
sx q[2];
rz(1.8075338) q[2];
rz(2.7110161) q[3];
sx q[3];
rz(-1.4686573) q[3];
sx q[3];
rz(2.8051643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
