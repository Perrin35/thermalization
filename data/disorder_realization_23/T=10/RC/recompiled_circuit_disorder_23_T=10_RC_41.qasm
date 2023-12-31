OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(6.8847818) q[1];
sx q[1];
rz(9.8431982) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.943676) q[0];
sx q[0];
rz(-2.0352053) q[0];
sx q[0];
rz(0.29505131) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8848959) q[2];
sx q[2];
rz(-0.17527097) q[2];
sx q[2];
rz(2.0435317) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5814372) q[1];
sx q[1];
rz(-2.0247702) q[1];
sx q[1];
rz(0.27172383) q[1];
x q[2];
rz(2.0166287) q[3];
sx q[3];
rz(-0.48005193) q[3];
sx q[3];
rz(-2.0208841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3216386) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(2.3036172) q[2];
rz(0.49301246) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3943966) q[0];
sx q[0];
rz(-0.72421873) q[0];
sx q[0];
rz(-1.863742) q[0];
rz(0.17678075) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(2.7094254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1457739) q[0];
sx q[0];
rz(-0.60129014) q[0];
sx q[0];
rz(2.6390618) q[0];
rz(-pi) q[1];
rz(0.82168174) q[2];
sx q[2];
rz(-0.5558008) q[2];
sx q[2];
rz(2.1205714) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.00694) q[1];
sx q[1];
rz(-1.68774) q[1];
sx q[1];
rz(-1.1275396) q[1];
rz(-1.6099986) q[3];
sx q[3];
rz(-0.56534846) q[3];
sx q[3];
rz(0.70985868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(-0.48669997) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.8816201) q[3];
sx q[3];
rz(0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(-2.7242463) q[0];
rz(-1.6529282) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(-0.506385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7931472) q[0];
sx q[0];
rz(-1.2198824) q[0];
sx q[0];
rz(-0.1883513) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.330553) q[2];
sx q[2];
rz(-2.2908387) q[2];
sx q[2];
rz(0.70149295) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80458927) q[1];
sx q[1];
rz(-2.5924006) q[1];
sx q[1];
rz(0.6188436) q[1];
rz(-pi) q[2];
rz(-1.1261602) q[3];
sx q[3];
rz(-0.4053084) q[3];
sx q[3];
rz(-2.0733548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4425519) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(0.58602035) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9716924) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(2.3024094) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(1.5930088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65028134) q[0];
sx q[0];
rz(-0.45931739) q[0];
sx q[0];
rz(1.3643054) q[0];
rz(-pi) q[1];
rz(0.71009212) q[2];
sx q[2];
rz(-0.93821628) q[2];
sx q[2];
rz(-1.8284423) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6099445) q[1];
sx q[1];
rz(-2.2713486) q[1];
sx q[1];
rz(2.080426) q[1];
rz(-pi) q[2];
rz(-2.5769916) q[3];
sx q[3];
rz(-0.49800107) q[3];
sx q[3];
rz(2.6548487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8362391) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(3.0419066) q[2];
rz(0.95885197) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.3249741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56617671) q[0];
sx q[0];
rz(-1.7594936) q[0];
sx q[0];
rz(-2.8856522) q[0];
rz(-0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(0.76006132) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5998659) q[0];
sx q[0];
rz(-1.4540298) q[0];
sx q[0];
rz(1.6745425) q[0];
x q[1];
rz(1.7145833) q[2];
sx q[2];
rz(-2.723105) q[2];
sx q[2];
rz(-2.8029122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1000423) q[1];
sx q[1];
rz(-1.4428382) q[1];
sx q[1];
rz(-3.0138739) q[1];
rz(-pi) q[2];
rz(2.7759336) q[3];
sx q[3];
rz(-0.95698157) q[3];
sx q[3];
rz(-2.1706276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(2.467353) q[2];
rz(-2.9267866) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(0.017344346) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53133416) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.1791139) q[0];
rz(2.9367661) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(-1.0669605) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2039295) q[0];
sx q[0];
rz(-1.5852889) q[0];
sx q[0];
rz(3.1209164) q[0];
rz(-pi) q[1];
rz(-2.009379) q[2];
sx q[2];
rz(-1.9595993) q[2];
sx q[2];
rz(-2.4593381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8384335) q[1];
sx q[1];
rz(-1.618297) q[1];
sx q[1];
rz(0.82364239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1596998) q[3];
sx q[3];
rz(-2.3387863) q[3];
sx q[3];
rz(-1.6646977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71010464) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(2.5816494) q[2];
rz(2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(-2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.2840435) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(-0.85817671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9265384) q[0];
sx q[0];
rz(-2.9460322) q[0];
sx q[0];
rz(-2.2901448) q[0];
rz(-1.467534) q[2];
sx q[2];
rz(-0.41245663) q[2];
sx q[2];
rz(0.39706424) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3837636) q[1];
sx q[1];
rz(-1.541829) q[1];
sx q[1];
rz(-1.5555192) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9167561) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(-0.35017761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8537366) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(-0.89007968) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(-2.1530698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(-0.37471399) q[0];
rz(0.97887865) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6357248) q[0];
sx q[0];
rz(-0.74698193) q[0];
sx q[0];
rz(-0.58734679) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1153203) q[2];
sx q[2];
rz(-2.4233344) q[2];
sx q[2];
rz(-0.19933137) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.046557758) q[1];
sx q[1];
rz(-1.7525502) q[1];
sx q[1];
rz(-0.05038105) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54578636) q[3];
sx q[3];
rz(-1.7958399) q[3];
sx q[3];
rz(0.033566098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4902041) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.7780001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63012183) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.2619031) q[0];
rz(-0.17503861) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(-1.6040241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14311929) q[0];
sx q[0];
rz(-0.38030085) q[0];
sx q[0];
rz(-0.93912504) q[0];
rz(-1.6749009) q[2];
sx q[2];
rz(-1.1541919) q[2];
sx q[2];
rz(2.5660851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4496085) q[1];
sx q[1];
rz(-2.2712628) q[1];
sx q[1];
rz(-1.4222172) q[1];
rz(0.92026199) q[3];
sx q[3];
rz(-2.5839621) q[3];
sx q[3];
rz(-0.44616163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1016772) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(-0.76134479) q[2];
rz(2.2411761) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(0.049023978) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(0.21690579) q[0];
rz(2.5096109) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(-2.1868618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8666329) q[0];
sx q[0];
rz(-1.8909591) q[0];
sx q[0];
rz(0.76036705) q[0];
x q[1];
rz(-1.1502613) q[2];
sx q[2];
rz(-2.0805801) q[2];
sx q[2];
rz(1.2812986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.875647) q[1];
sx q[1];
rz(-1.6469643) q[1];
sx q[1];
rz(0.45653685) q[1];
rz(-pi) q[2];
rz(-0.44317742) q[3];
sx q[3];
rz(-0.9842397) q[3];
sx q[3];
rz(-0.12936684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0163991) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(-0.94474244) q[2];
rz(-2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(-2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64086296) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(0.89818556) q[2];
sx q[2];
rz(-2.2834416) q[2];
sx q[2];
rz(-1.7590547) q[2];
rz(2.0142043) q[3];
sx q[3];
rz(-1.3848806) q[3];
sx q[3];
rz(1.0564907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
