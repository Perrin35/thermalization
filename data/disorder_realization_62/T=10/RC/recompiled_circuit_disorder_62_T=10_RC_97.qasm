OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(4.1783279) q[1];
sx q[1];
rz(8.7611603) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9288869) q[0];
sx q[0];
rz(-1.1097254) q[0];
sx q[0];
rz(1.5972208) q[0];
rz(-pi) q[1];
rz(2.3637949) q[2];
sx q[2];
rz(-1.8282837) q[2];
sx q[2];
rz(-2.4553026) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9259778) q[1];
sx q[1];
rz(-1.8324592) q[1];
sx q[1];
rz(-1.91933) q[1];
x q[2];
rz(-3.033028) q[3];
sx q[3];
rz(-1.8130455) q[3];
sx q[3];
rz(3.0755175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87876451) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(-3.0905241) q[2];
rz(2.5845394) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(-0.59659514) q[0];
rz(-0.82582981) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.2260431) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37218371) q[0];
sx q[0];
rz(-0.49159494) q[0];
sx q[0];
rz(-0.43222897) q[0];
x q[1];
rz(0.82310279) q[2];
sx q[2];
rz(-0.96379333) q[2];
sx q[2];
rz(-2.5603103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.021124161) q[1];
sx q[1];
rz(-1.7906584) q[1];
sx q[1];
rz(3.1411509) q[1];
rz(3.1321208) q[3];
sx q[3];
rz(-0.99820271) q[3];
sx q[3];
rz(-0.026281683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(2.8105695) q[2];
rz(0.80667574) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9817292) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(1.249041) q[0];
rz(0.088009134) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(-2.1121315) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7824161) q[0];
sx q[0];
rz(-0.65890197) q[0];
sx q[0];
rz(1.2078148) q[0];
rz(-pi) q[1];
rz(1.8573895) q[2];
sx q[2];
rz(-2.2006052) q[2];
sx q[2];
rz(-2.0558002) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50476915) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(1.7117281) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73996468) q[3];
sx q[3];
rz(-1.4810586) q[3];
sx q[3];
rz(1.2147853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(2.0987299) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.625659) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0632616) q[0];
sx q[0];
rz(-0.69561361) q[0];
sx q[0];
rz(-2.9177833) q[0];
rz(-pi) q[1];
rz(-1.0117202) q[2];
sx q[2];
rz(-1.498073) q[2];
sx q[2];
rz(2.5474472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8856636) q[1];
sx q[1];
rz(-2.4293578) q[1];
sx q[1];
rz(-2.7182012) q[1];
rz(-pi) q[2];
rz(1.3514148) q[3];
sx q[3];
rz(-0.99321584) q[3];
sx q[3];
rz(-2.1882309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76379124) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.8466922) q[2];
rz(3.0002248) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-0.98658371) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(-1.6217344) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(2.246726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1631854) q[0];
sx q[0];
rz(-2.1911231) q[0];
sx q[0];
rz(2.083367) q[0];
rz(-0.72563719) q[2];
sx q[2];
rz(-2.0113532) q[2];
sx q[2];
rz(0.63647905) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4836854) q[1];
sx q[1];
rz(-0.92805082) q[1];
sx q[1];
rz(0.55773463) q[1];
rz(-3.1154247) q[3];
sx q[3];
rz(-2.8550365) q[3];
sx q[3];
rz(0.22931306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52508369) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(-1.9990702) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(-1.07553) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557945) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.7640132) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-0.46498743) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2842448) q[0];
sx q[0];
rz(-1.2372969) q[0];
sx q[0];
rz(-1.1112945) q[0];
rz(-2.8000185) q[2];
sx q[2];
rz(-1.5014868) q[2];
sx q[2];
rz(2.9937033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28593674) q[1];
sx q[1];
rz(-0.93614139) q[1];
sx q[1];
rz(-0.059307701) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24352169) q[3];
sx q[3];
rz(-0.72391073) q[3];
sx q[3];
rz(0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.650699) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-2.6965551) q[2];
rz(-0.93368357) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(1.4200462) q[0];
rz(-3.1177915) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(0.15596095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1056846) q[0];
sx q[0];
rz(-1.7419635) q[0];
sx q[0];
rz(0.25733421) q[0];
rz(2.4620352) q[2];
sx q[2];
rz(-2.5153) q[2];
sx q[2];
rz(1.9192413) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1289039) q[1];
sx q[1];
rz(-1.3959937) q[1];
sx q[1];
rz(3.0867982) q[1];
rz(-pi) q[2];
rz(-2.3994) q[3];
sx q[3];
rz(-0.67897292) q[3];
sx q[3];
rz(1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-2.0689266) q[2];
rz(-0.32564751) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496562) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-0.33777133) q[0];
rz(-1.0900963) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(2.8930194) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2598495) q[0];
sx q[0];
rz(-1.3215085) q[0];
sx q[0];
rz(2.6134406) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3209004) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(-1.5936268) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3358811) q[1];
sx q[1];
rz(-2.1599249) q[1];
sx q[1];
rz(2.5475403) q[1];
rz(1.7406086) q[3];
sx q[3];
rz(-1.6465934) q[3];
sx q[3];
rz(0.71789391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2934072) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(0.48197204) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329353) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(0.70156082) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(-1.3185906) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628172) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(-1.1286939) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4039688) q[2];
sx q[2];
rz(-0.80543033) q[2];
sx q[2];
rz(-1.9464303) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9642155) q[1];
sx q[1];
rz(-0.69987684) q[1];
sx q[1];
rz(-1.7978976) q[1];
rz(1.9570458) q[3];
sx q[3];
rz(-0.85856122) q[3];
sx q[3];
rz(1.6561179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7302154) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(2.7424157) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(-1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(1.8803966) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3189145) q[0];
sx q[0];
rz(-0.42352522) q[0];
sx q[0];
rz(-2.8828997) q[0];
rz(-pi) q[1];
rz(1.5930428) q[2];
sx q[2];
rz(-0.97264475) q[2];
sx q[2];
rz(2.1200137) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2640472) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(1.643814) q[1];
rz(-pi) q[2];
rz(-1.7301171) q[3];
sx q[3];
rz(-0.89309249) q[3];
sx q[3];
rz(-0.89980984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(-1.9231046) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8284843) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(3.1153395) q[2];
sx q[2];
rz(-0.99925169) q[2];
sx q[2];
rz(3.0796438) q[2];
rz(-0.20523397) q[3];
sx q[3];
rz(-0.60094613) q[3];
sx q[3];
rz(-0.080106674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
