OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4322296) q[0];
sx q[0];
rz(-0.95786434) q[0];
sx q[0];
rz(0.14444484) q[0];
rz(-2.5748409) q[1];
sx q[1];
rz(-2.6161939) q[1];
sx q[1];
rz(2.1638343) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.190783) q[0];
sx q[0];
rz(-1.8291744) q[0];
sx q[0];
rz(-1.4825312) q[0];
x q[1];
rz(2.0830886) q[2];
sx q[2];
rz(-2.5663178) q[2];
sx q[2];
rz(3.1303867) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49330083) q[1];
sx q[1];
rz(-0.19953218) q[1];
sx q[1];
rz(-0.7944016) q[1];
rz(-pi) q[2];
rz(2.9631859) q[3];
sx q[3];
rz(-0.79531407) q[3];
sx q[3];
rz(2.2291396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4686761) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(-2.9428234) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7070049) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(-0.36112753) q[0];
rz(-1.7065642) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(-0.8180058) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084413962) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(-1.8684698) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0188574) q[2];
sx q[2];
rz(-1.2018179) q[2];
sx q[2];
rz(-2.4462647) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6815731) q[1];
sx q[1];
rz(-2.1732554) q[1];
sx q[1];
rz(0.54089344) q[1];
rz(-pi) q[2];
rz(3.0844968) q[3];
sx q[3];
rz(-1.207587) q[3];
sx q[3];
rz(-0.2163987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.25257418) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(-1.4206295) q[2];
rz(-1.8255) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(-0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.4017568) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(-0.87093583) q[0];
rz(-0.31618205) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(-2.8443764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0270099) q[0];
sx q[0];
rz(-1.6342589) q[0];
sx q[0];
rz(-0.40513904) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88163968) q[2];
sx q[2];
rz(-2.214553) q[2];
sx q[2];
rz(1.8635441) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3618468) q[1];
sx q[1];
rz(-1.5231113) q[1];
sx q[1];
rz(-1.7576799) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66629569) q[3];
sx q[3];
rz(-2.1994281) q[3];
sx q[3];
rz(-0.48376885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7685984) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(2.7139943) q[2];
rz(1.1887431) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.7052085) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(0.77392459) q[0];
rz(0.71290839) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(-2.4598222) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49127689) q[0];
sx q[0];
rz(-0.33704764) q[0];
sx q[0];
rz(-2.2469673) q[0];
rz(1.2615471) q[2];
sx q[2];
rz(-2.5144858) q[2];
sx q[2];
rz(-2.4583465) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3649498) q[1];
sx q[1];
rz(-0.72307359) q[1];
sx q[1];
rz(-0.88100453) q[1];
rz(-pi) q[2];
rz(-1.2069615) q[3];
sx q[3];
rz(-1.2483276) q[3];
sx q[3];
rz(-2.913167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(-0.95820367) q[2];
rz(-2.0751674) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(-0.56030309) q[0];
rz(0.99984461) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(1.5195742) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27846562) q[0];
sx q[0];
rz(-1.0882049) q[0];
sx q[0];
rz(2.0762073) q[0];
x q[1];
rz(1.0137453) q[2];
sx q[2];
rz(-1.0336913) q[2];
sx q[2];
rz(-2.277166) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6474364) q[1];
sx q[1];
rz(-1.7687706) q[1];
sx q[1];
rz(-0.76311771) q[1];
x q[2];
rz(0.95789692) q[3];
sx q[3];
rz(-1.6727722) q[3];
sx q[3];
rz(1.0979872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6891629) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(-2.9186644) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(-0.071578659) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3361622) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(2.0715332) q[0];
rz(-1.3609715) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(-1.8575352) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43858466) q[0];
sx q[0];
rz(-2.2692338) q[0];
sx q[0];
rz(2.3914778) q[0];
x q[1];
rz(0.96096303) q[2];
sx q[2];
rz(-0.46513882) q[2];
sx q[2];
rz(-2.7679408) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2594086) q[1];
sx q[1];
rz(-0.55667294) q[1];
sx q[1];
rz(0.56682079) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1070126) q[3];
sx q[3];
rz(-1.8367193) q[3];
sx q[3];
rz(-1.769161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47508919) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(-0.63759032) q[2];
rz(0.83667886) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(-1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5752983) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(-2.7138846) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(2.2033851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8021916) q[0];
sx q[0];
rz(-1.2656478) q[0];
sx q[0];
rz(1.9992273) q[0];
rz(-pi) q[1];
rz(2.1392518) q[2];
sx q[2];
rz(-0.92407862) q[2];
sx q[2];
rz(-1.9288174) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2141014) q[1];
sx q[1];
rz(-2.2117105) q[1];
sx q[1];
rz(-0.44968857) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35258099) q[3];
sx q[3];
rz(-2.4028824) q[3];
sx q[3];
rz(1.1383566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2618835) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(1.7555457) q[2];
rz(1.322768) q[3];
sx q[3];
rz(-1.991792) q[3];
sx q[3];
rz(-0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8742074) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(1.7077131) q[0];
rz(-1.8677615) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(2.3103255) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87330504) q[0];
sx q[0];
rz(-1.1697066) q[0];
sx q[0];
rz(-1.7307161) q[0];
x q[1];
rz(0.12886329) q[2];
sx q[2];
rz(-2.2017751) q[2];
sx q[2];
rz(-1.250759) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23849328) q[1];
sx q[1];
rz(-2.0598754) q[1];
sx q[1];
rz(-2.3563983) q[1];
rz(-pi) q[2];
rz(2.9112178) q[3];
sx q[3];
rz(-1.3532234) q[3];
sx q[3];
rz(2.2090467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5635809) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-0.9643628) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(-2.4826629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37373856) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(-2.1642165) q[0];
rz(1.7550229) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(1.1057373) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39767299) q[0];
sx q[0];
rz(-2.0154675) q[0];
sx q[0];
rz(-1.7001274) q[0];
rz(-0.90412037) q[2];
sx q[2];
rz(-1.4881926) q[2];
sx q[2];
rz(-2.3711575) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3016262) q[1];
sx q[1];
rz(-1.7558388) q[1];
sx q[1];
rz(-1.3355096) q[1];
rz(-pi) q[2];
rz(-1.2826142) q[3];
sx q[3];
rz(-1.9087221) q[3];
sx q[3];
rz(-1.5087138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(0.14979714) q[2];
rz(1.7685361) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(-1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(0.1272442) q[0];
rz(1.4808902) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(-0.1677992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084899336) q[0];
sx q[0];
rz(-0.035824422) q[0];
sx q[0];
rz(0.56532677) q[0];
rz(-pi) q[1];
rz(0.22044214) q[2];
sx q[2];
rz(-0.61969212) q[2];
sx q[2];
rz(0.1534136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2044636) q[1];
sx q[1];
rz(-0.58848721) q[1];
sx q[1];
rz(-0.012339331) q[1];
x q[2];
rz(3.0182748) q[3];
sx q[3];
rz(-1.6880369) q[3];
sx q[3];
rz(0.58837147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1148791) q[2];
sx q[2];
rz(-2.202704) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(-3.051565) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(-0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5933843) q[0];
sx q[0];
rz(-1.1593288) q[0];
sx q[0];
rz(2.819084) q[0];
rz(-0.38800115) q[1];
sx q[1];
rz(-1.3996268) q[1];
sx q[1];
rz(-0.7849801) q[1];
rz(-2.3315196) q[2];
sx q[2];
rz(-1.0675061) q[2];
sx q[2];
rz(1.6455417) q[2];
rz(-2.9705863) q[3];
sx q[3];
rz(-1.0382367) q[3];
sx q[3];
rz(2.3259179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
