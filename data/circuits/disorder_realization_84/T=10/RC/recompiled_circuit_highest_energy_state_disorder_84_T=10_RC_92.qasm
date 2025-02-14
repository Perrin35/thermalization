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
rz(-2.5007024) q[0];
rz(-0.23735292) q[1];
sx q[1];
rz(3.7832398) q[1];
sx q[1];
rz(6.1982815) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9083402) q[0];
sx q[0];
rz(-1.4734992) q[0];
sx q[0];
rz(0.86546524) q[0];
rz(-1.6019967) q[2];
sx q[2];
rz(-1.5123431) q[2];
sx q[2];
rz(1.8606869) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81604874) q[1];
sx q[1];
rz(-2.8286165) q[1];
sx q[1];
rz(-0.37892394) q[1];
rz(-pi) q[2];
rz(-2.1617684) q[3];
sx q[3];
rz(-1.9035467) q[3];
sx q[3];
rz(-2.7036288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25195965) q[2];
sx q[2];
rz(-1.2026938) q[2];
sx q[2];
rz(1.5084051) q[2];
rz(-2.3097322) q[3];
sx q[3];
rz(-0.32604495) q[3];
sx q[3];
rz(-1.3445725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2107596) q[0];
sx q[0];
rz(-0.63527125) q[0];
sx q[0];
rz(0.26500901) q[0];
rz(-0.79568544) q[1];
sx q[1];
rz(-0.34514752) q[1];
sx q[1];
rz(-1.4211242) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1596589) q[0];
sx q[0];
rz(-2.2082885) q[0];
sx q[0];
rz(0.33099126) q[0];
x q[1];
rz(-2.1876343) q[2];
sx q[2];
rz(-1.1797311) q[2];
sx q[2];
rz(0.80634889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15884033) q[1];
sx q[1];
rz(-1.6144061) q[1];
sx q[1];
rz(-1.8226472) q[1];
x q[2];
rz(1.4031409) q[3];
sx q[3];
rz(-1.5412313) q[3];
sx q[3];
rz(0.61691446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77218324) q[2];
sx q[2];
rz(-2.1464244) q[2];
sx q[2];
rz(0.65917242) q[2];
rz(-2.2199421) q[3];
sx q[3];
rz(-2.0989336) q[3];
sx q[3];
rz(1.5590182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7121861) q[0];
sx q[0];
rz(-0.61539188) q[0];
sx q[0];
rz(-1.8581101) q[0];
rz(0.42690024) q[1];
sx q[1];
rz(-2.5476397) q[1];
sx q[1];
rz(-1.1675507) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3385462) q[0];
sx q[0];
rz(-1.2075222) q[0];
sx q[0];
rz(2.7873894) q[0];
x q[1];
rz(-2.3857791) q[2];
sx q[2];
rz(-2.4632556) q[2];
sx q[2];
rz(-2.8752022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7181212) q[1];
sx q[1];
rz(-1.2813066) q[1];
sx q[1];
rz(-1.4726588) q[1];
rz(-pi) q[2];
rz(1.394672) q[3];
sx q[3];
rz(-1.307789) q[3];
sx q[3];
rz(0.60282487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.20649642) q[2];
sx q[2];
rz(-1.6556211) q[2];
sx q[2];
rz(3.0237954) q[2];
rz(1.8360893) q[3];
sx q[3];
rz(-2.2820303) q[3];
sx q[3];
rz(-2.1507202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8369668) q[0];
sx q[0];
rz(-1.0517629) q[0];
sx q[0];
rz(3.1090609) q[0];
rz(-2.1578728) q[1];
sx q[1];
rz(-0.81147057) q[1];
sx q[1];
rz(-0.54289114) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6257831) q[0];
sx q[0];
rz(-2.2383732) q[0];
sx q[0];
rz(-1.9208917) q[0];
rz(3.0383598) q[2];
sx q[2];
rz(-1.8816173) q[2];
sx q[2];
rz(2.9048267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9223843) q[1];
sx q[1];
rz(-2.3047949) q[1];
sx q[1];
rz(0.34511415) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.832295) q[3];
sx q[3];
rz(-1.4595539) q[3];
sx q[3];
rz(-0.55505607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5797609) q[2];
sx q[2];
rz(-2.0115972) q[2];
sx q[2];
rz(-0.96060166) q[2];
rz(-2.1660755) q[3];
sx q[3];
rz(-2.2701008) q[3];
sx q[3];
rz(0.52198854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.3112711) q[0];
sx q[0];
rz(-0.68205849) q[0];
sx q[0];
rz(-0.3325381) q[0];
rz(0.63811103) q[1];
sx q[1];
rz(-2.5159914) q[1];
sx q[1];
rz(-2.6854551) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79955937) q[0];
sx q[0];
rz(-1.7029951) q[0];
sx q[0];
rz(3.0782736) q[0];
x q[1];
rz(2.3535185) q[2];
sx q[2];
rz(-1.9219134) q[2];
sx q[2];
rz(1.9478363) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3553473) q[1];
sx q[1];
rz(-1.659871) q[1];
sx q[1];
rz(0.98615714) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92335574) q[3];
sx q[3];
rz(-1.8643987) q[3];
sx q[3];
rz(-1.7838316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0786324) q[2];
sx q[2];
rz(-0.76138622) q[2];
sx q[2];
rz(-2.6893993) q[2];
rz(-2.8771628) q[3];
sx q[3];
rz(-2.9954438) q[3];
sx q[3];
rz(1.9047846) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0662769) q[0];
sx q[0];
rz(-0.83106581) q[0];
sx q[0];
rz(-0.7178632) q[0];
rz(1.5526937) q[1];
sx q[1];
rz(-1.6009067) q[1];
sx q[1];
rz(-2.9343361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.358905) q[0];
sx q[0];
rz(-1.1289294) q[0];
sx q[0];
rz(-1.601041) q[0];
rz(1.7278094) q[2];
sx q[2];
rz(-0.42745379) q[2];
sx q[2];
rz(0.70892109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9178324) q[1];
sx q[1];
rz(-2.4759935) q[1];
sx q[1];
rz(-0.83719801) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9620216) q[3];
sx q[3];
rz(-0.28124725) q[3];
sx q[3];
rz(-2.1811671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3805716) q[2];
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
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2137432) q[0];
sx q[0];
rz(-2.1708467) q[0];
sx q[0];
rz(0.84681502) q[0];
rz(-0.98833409) q[1];
sx q[1];
rz(-1.7653932) q[1];
sx q[1];
rz(-0.83289897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083713934) q[0];
sx q[0];
rz(-1.3728317) q[0];
sx q[0];
rz(-1.676487) q[0];
rz(-pi) q[1];
rz(-1.6498327) q[2];
sx q[2];
rz(-1.7124868) q[2];
sx q[2];
rz(-2.2159383) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.78645241) q[1];
sx q[1];
rz(-1.2017829) q[1];
sx q[1];
rz(-2.2166316) q[1];
x q[2];
rz(-1.9774632) q[3];
sx q[3];
rz(-3.0179438) q[3];
sx q[3];
rz(0.37227889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6758468) q[2];
sx q[2];
rz(-2.4387359) q[2];
sx q[2];
rz(-0.77009002) q[2];
rz(-3.0105524) q[3];
sx q[3];
rz(-1.659487) q[3];
sx q[3];
rz(1.2469163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34762621) q[0];
sx q[0];
rz(-0.93526953) q[0];
sx q[0];
rz(0.4916218) q[0];
rz(2.8288016) q[1];
sx q[1];
rz(-0.28952315) q[1];
sx q[1];
rz(-3.0591931) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6254344) q[0];
sx q[0];
rz(-1.5835198) q[0];
sx q[0];
rz(-1.660343) q[0];
rz(-pi) q[1];
rz(-0.56677197) q[2];
sx q[2];
rz(-0.9563891) q[2];
sx q[2];
rz(0.64726171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1090549) q[1];
sx q[1];
rz(-1.7454285) q[1];
sx q[1];
rz(2.9777378) q[1];
rz(-0.1020482) q[3];
sx q[3];
rz(-1.2597435) q[3];
sx q[3];
rz(-1.2132258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6975434) q[2];
sx q[2];
rz(-1.9335582) q[2];
sx q[2];
rz(3.1235798) q[2];
rz(-0.87120122) q[3];
sx q[3];
rz(-0.33659354) q[3];
sx q[3];
rz(2.474031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
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
rz(-0.95053259) q[0];
sx q[0];
rz(-2.58707) q[0];
sx q[0];
rz(2.6771123) q[0];
rz(0.49098268) q[1];
sx q[1];
rz(-1.4717439) q[1];
sx q[1];
rz(1.6741265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1037285) q[0];
sx q[0];
rz(-1.5707644) q[0];
sx q[0];
rz(0.62490827) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6355714) q[2];
sx q[2];
rz(-1.3014761) q[2];
sx q[2];
rz(1.8647746) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3882093) q[1];
sx q[1];
rz(-2.8860554) q[1];
sx q[1];
rz(1.3292043) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2676226) q[3];
sx q[3];
rz(-2.2010816) q[3];
sx q[3];
rz(1.9809679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1269425) q[2];
sx q[2];
rz(-1.6501959) q[2];
sx q[2];
rz(-1.3313741) q[2];
rz(0.93975449) q[3];
sx q[3];
rz(-1.0569812) q[3];
sx q[3];
rz(-1.1872928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2347655) q[0];
sx q[0];
rz(-2.6563788) q[0];
sx q[0];
rz(1.0377129) q[0];
rz(1.0543793) q[1];
sx q[1];
rz(-1.3155921) q[1];
sx q[1];
rz(-2.0920928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9535372) q[0];
sx q[0];
rz(-0.049727289) q[0];
sx q[0];
rz(-1.5744857) q[0];
rz(-pi) q[1];
rz(0.56042258) q[2];
sx q[2];
rz(-1.9935521) q[2];
sx q[2];
rz(-2.1566856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8423374) q[1];
sx q[1];
rz(-1.6016869) q[1];
sx q[1];
rz(0.96986356) q[1];
rz(0.78759362) q[3];
sx q[3];
rz(-0.75852048) q[3];
sx q[3];
rz(0.99419981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1513169) q[2];
sx q[2];
rz(-1.8345202) q[2];
sx q[2];
rz(-1.6952197) q[2];
rz(0.15898786) q[3];
sx q[3];
rz(-1.9828601) q[3];
sx q[3];
rz(-0.75001636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1405519) q[0];
sx q[0];
rz(-0.97698553) q[0];
sx q[0];
rz(-2.3401596) q[0];
rz(1.6544381) q[1];
sx q[1];
rz(-0.14688891) q[1];
sx q[1];
rz(2.6085703) q[1];
rz(0.64143388) q[2];
sx q[2];
rz(-0.95856711) q[2];
sx q[2];
rz(-2.7480456) q[2];
rz(0.63538649) q[3];
sx q[3];
rz(-1.0839331) q[3];
sx q[3];
rz(-1.4270368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
