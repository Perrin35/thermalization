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
rz(-2.6211205) q[0];
sx q[0];
rz(-1.0360798) q[0];
sx q[0];
rz(0.64089027) q[0];
rz(-0.23735292) q[1];
sx q[1];
rz(3.7832398) q[1];
sx q[1];
rz(6.1982815) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9176506) q[0];
sx q[0];
rz(-0.71086796) q[0];
sx q[0];
rz(1.4213597) q[0];
rz(-pi) q[1];
rz(-0.48975421) q[2];
sx q[2];
rz(-0.066250525) q[2];
sx q[2];
rz(0.79023933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0247268) q[1];
sx q[1];
rz(-1.6849396) q[1];
sx q[1];
rz(0.29205871) q[1];
x q[2];
rz(-0.39438168) q[3];
sx q[3];
rz(-2.1254) q[3];
sx q[3];
rz(-1.7930052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25195965) q[2];
sx q[2];
rz(-1.9388988) q[2];
sx q[2];
rz(-1.6331875) q[2];
rz(2.3097322) q[3];
sx q[3];
rz(-0.32604495) q[3];
sx q[3];
rz(1.3445725) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93083301) q[0];
sx q[0];
rz(-2.5063214) q[0];
sx q[0];
rz(0.26500901) q[0];
rz(2.3459072) q[1];
sx q[1];
rz(-2.7964451) q[1];
sx q[1];
rz(-1.7204684) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1596589) q[0];
sx q[0];
rz(-2.2082885) q[0];
sx q[0];
rz(-0.33099126) q[0];
rz(-0.9515597) q[2];
sx q[2];
rz(-0.71647055) q[2];
sx q[2];
rz(1.2576511) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7184192) q[1];
sx q[1];
rz(-1.8224026) q[1];
sx q[1];
rz(-0.045028506) q[1];
rz(-pi) q[2];
rz(-1.4031409) q[3];
sx q[3];
rz(-1.5412313) q[3];
sx q[3];
rz(-0.61691446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77218324) q[2];
sx q[2];
rz(-0.99516827) q[2];
sx q[2];
rz(-2.4824202) q[2];
rz(2.2199421) q[3];
sx q[3];
rz(-2.0989336) q[3];
sx q[3];
rz(-1.5590182) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42940656) q[0];
sx q[0];
rz(-2.5262008) q[0];
sx q[0];
rz(1.8581101) q[0];
rz(2.7146924) q[1];
sx q[1];
rz(-0.59395298) q[1];
sx q[1];
rz(1.9740419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1438863) q[0];
sx q[0];
rz(-0.50181118) q[0];
sx q[0];
rz(-2.3104179) q[0];
rz(2.0757572) q[2];
sx q[2];
rz(-2.0450108) q[2];
sx q[2];
rz(1.9948639) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0223654) q[1];
sx q[1];
rz(-1.664838) q[1];
sx q[1];
rz(0.29081197) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.394672) q[3];
sx q[3];
rz(-1.307789) q[3];
sx q[3];
rz(2.5387678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.20649642) q[2];
sx q[2];
rz(-1.4859716) q[2];
sx q[2];
rz(3.0237954) q[2];
rz(-1.8360893) q[3];
sx q[3];
rz(-0.85956231) q[3];
sx q[3];
rz(0.99087244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369668) q[0];
sx q[0];
rz(-2.0898297) q[0];
sx q[0];
rz(3.1090609) q[0];
rz(2.1578728) q[1];
sx q[1];
rz(-0.81147057) q[1];
sx q[1];
rz(-2.5987015) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6257831) q[0];
sx q[0];
rz(-2.2383732) q[0];
sx q[0];
rz(1.9208917) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2603733) q[2];
sx q[2];
rz(-2.8146048) q[2];
sx q[2];
rz(0.56337683) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11527744) q[1];
sx q[1];
rz(-1.3168854) q[1];
sx q[1];
rz(0.8064958) q[1];
rz(2.9916688) q[3];
sx q[3];
rz(-2.3036851) q[3];
sx q[3];
rz(-0.91502193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5797609) q[2];
sx q[2];
rz(-1.1299955) q[2];
sx q[2];
rz(2.180991) q[2];
rz(2.1660755) q[3];
sx q[3];
rz(-2.2701008) q[3];
sx q[3];
rz(2.6196041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.3112711) q[0];
sx q[0];
rz(-0.68205849) q[0];
sx q[0];
rz(-0.3325381) q[0];
rz(2.5034816) q[1];
sx q[1];
rz(-2.5159914) q[1];
sx q[1];
rz(2.6854551) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35121954) q[0];
sx q[0];
rz(-2.995092) q[0];
sx q[0];
rz(-2.0149433) q[0];
x q[1];
rz(0.47686968) q[2];
sx q[2];
rz(-2.2945291) q[2];
sx q[2];
rz(-0.70666955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3553473) q[1];
sx q[1];
rz(-1.4817217) q[1];
sx q[1];
rz(-0.98615714) q[1];
rz(-pi) q[2];
rz(-2.0354603) q[3];
sx q[3];
rz(-0.70208462) q[3];
sx q[3];
rz(-0.15248577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0629603) q[2];
sx q[2];
rz(-2.3802064) q[2];
sx q[2];
rz(-2.6893993) q[2];
rz(2.8771628) q[3];
sx q[3];
rz(-2.9954438) q[3];
sx q[3];
rz(1.2368081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.0662769) q[0];
sx q[0];
rz(-2.3105268) q[0];
sx q[0];
rz(2.4237295) q[0];
rz(1.5526937) q[1];
sx q[1];
rz(-1.540686) q[1];
sx q[1];
rz(2.9343361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22482797) q[0];
sx q[0];
rz(-1.5981354) q[0];
sx q[0];
rz(0.44204373) q[0];
x q[1];
rz(-1.147993) q[2];
sx q[2];
rz(-1.6356653) q[2];
sx q[2];
rz(-2.1366304) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2237602) q[1];
sx q[1];
rz(-2.4759935) q[1];
sx q[1];
rz(0.83719801) q[1];
rz(-pi) q[2];
rz(3.0318694) q[3];
sx q[3];
rz(-1.3113113) q[3];
sx q[3];
rz(-2.5867227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76102105) q[2];
sx q[2];
rz(-2.0863159) q[2];
sx q[2];
rz(0.5611678) q[2];
rz(-1.0593972) q[3];
sx q[3];
rz(-2.0249517) q[3];
sx q[3];
rz(1.4580956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2137432) q[0];
sx q[0];
rz(-2.1708467) q[0];
sx q[0];
rz(-0.84681502) q[0];
rz(0.98833409) q[1];
sx q[1];
rz(-1.7653932) q[1];
sx q[1];
rz(0.83289897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578787) q[0];
sx q[0];
rz(-1.3728317) q[0];
sx q[0];
rz(-1.676487) q[0];
rz(-0.14212823) q[2];
sx q[2];
rz(-1.649039) q[2];
sx q[2];
rz(-2.4852666) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3551402) q[1];
sx q[1];
rz(-1.2017829) q[1];
sx q[1];
rz(0.92496101) q[1];
rz(3.0924721) q[3];
sx q[3];
rz(-1.6843154) q[3];
sx q[3];
rz(-3.1044131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46574584) q[2];
sx q[2];
rz(-2.4387359) q[2];
sx q[2];
rz(-2.3715026) q[2];
rz(3.0105524) q[3];
sx q[3];
rz(-1.659487) q[3];
sx q[3];
rz(-1.2469163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939664) q[0];
sx q[0];
rz(-0.93526953) q[0];
sx q[0];
rz(-2.6499709) q[0];
rz(-2.8288016) q[1];
sx q[1];
rz(-2.8520695) q[1];
sx q[1];
rz(0.082399592) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9461878) q[0];
sx q[0];
rz(-0.090443693) q[0];
sx q[0];
rz(-1.4294595) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5748207) q[2];
sx q[2];
rz(-0.9563891) q[2];
sx q[2];
rz(0.64726171) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0325378) q[1];
sx q[1];
rz(-1.3961642) q[1];
sx q[1];
rz(2.9777378) q[1];
rz(-pi) q[2];
rz(-1.2639323) q[3];
sx q[3];
rz(-2.8147449) q[3];
sx q[3];
rz(-2.2512448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6975434) q[2];
sx q[2];
rz(-1.2080344) q[2];
sx q[2];
rz(-0.018012878) q[2];
rz(0.87120122) q[3];
sx q[3];
rz(-2.8049991) q[3];
sx q[3];
rz(2.474031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95053259) q[0];
sx q[0];
rz(-2.58707) q[0];
sx q[0];
rz(0.46448034) q[0];
rz(-0.49098268) q[1];
sx q[1];
rz(-1.6698488) q[1];
sx q[1];
rz(-1.4674662) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6086836) q[0];
sx q[0];
rz(-0.94588806) q[0];
sx q[0];
rz(1.5708357) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2651132) q[2];
sx q[2];
rz(-2.0569498) q[2];
sx q[2];
rz(-2.7012239) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1388359) q[1];
sx q[1];
rz(-1.3228387) q[1];
sx q[1];
rz(3.0791705) q[1];
rz(0.76045658) q[3];
sx q[3];
rz(-2.1158614) q[3];
sx q[3];
rz(2.2732002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1269425) q[2];
sx q[2];
rz(-1.6501959) q[2];
sx q[2];
rz(-1.3313741) q[2];
rz(2.2018382) q[3];
sx q[3];
rz(-2.0846114) q[3];
sx q[3];
rz(1.9542998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2347655) q[0];
sx q[0];
rz(-0.48521388) q[0];
sx q[0];
rz(1.0377129) q[0];
rz(1.0543793) q[1];
sx q[1];
rz(-1.3155921) q[1];
sx q[1];
rz(1.0494999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1880555) q[0];
sx q[0];
rz(-0.049727289) q[0];
sx q[0];
rz(1.567107) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0590354) q[2];
sx q[2];
rz(-2.0768713) q[2];
sx q[2];
rz(-2.807694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8423374) q[1];
sx q[1];
rz(-1.6016869) q[1];
sx q[1];
rz(0.96986356) q[1];
rz(-pi) q[2];
x q[2];
rz(2.353999) q[3];
sx q[3];
rz(-2.3830722) q[3];
sx q[3];
rz(-2.1473928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1513169) q[2];
sx q[2];
rz(-1.3070725) q[2];
sx q[2];
rz(1.6952197) q[2];
rz(-2.9826048) q[3];
sx q[3];
rz(-1.9828601) q[3];
sx q[3];
rz(2.3915763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0010407) q[0];
sx q[0];
rz(-2.1646071) q[0];
sx q[0];
rz(0.80143308) q[0];
rz(1.4871545) q[1];
sx q[1];
rz(-2.9947037) q[1];
sx q[1];
rz(-0.53302232) q[1];
rz(-2.2764789) q[2];
sx q[2];
rz(-0.8556753) q[2];
sx q[2];
rz(-0.52134261) q[2];
rz(0.98900908) q[3];
sx q[3];
rz(-1.0186356) q[3];
sx q[3];
rz(2.9531425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
