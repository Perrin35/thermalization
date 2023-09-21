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
rz(-1.707466) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(-2.5399962) q[1];
sx q[1];
rz(2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9040065) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(-2.0531274) q[0];
rz(-pi) q[1];
rz(1.7069874) q[2];
sx q[2];
rz(-1.681466) q[2];
sx q[2];
rz(1.1510804) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0093706) q[1];
sx q[1];
rz(-1.8144061) q[1];
sx q[1];
rz(1.1019215) q[1];
x q[2];
rz(2.9207346) q[3];
sx q[3];
rz(-2.0005895) q[3];
sx q[3];
rz(-0.62648279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3216386) q[2];
sx q[2];
rz(-1.7724089) q[2];
sx q[2];
rz(2.3036172) q[2];
rz(2.6485802) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74719602) q[0];
sx q[0];
rz(-0.72421873) q[0];
sx q[0];
rz(1.2778506) q[0];
rz(2.9648119) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(-2.7094254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85053274) q[0];
sx q[0];
rz(-1.2948372) q[0];
sx q[0];
rz(0.54131298) q[0];
rz(1.1439267) q[2];
sx q[2];
rz(-1.2032713) q[2];
sx q[2];
rz(1.2183684) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.00694) q[1];
sx q[1];
rz(-1.68774) q[1];
sx q[1];
rz(1.1275396) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1357972) q[3];
sx q[3];
rz(-1.5917935) q[3];
sx q[3];
rz(-0.82783031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54923487) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(0.48669997) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(1.6529282) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(-2.6352077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15784141) q[0];
sx q[0];
rz(-2.7451773) q[0];
sx q[0];
rz(1.0979963) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26489139) q[2];
sx q[2];
rz(-2.3893917) q[2];
sx q[2];
rz(1.0571935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8295146) q[1];
sx q[1];
rz(-1.2631589) q[1];
sx q[1];
rz(0.46244669) q[1];
x q[2];
rz(0.18249986) q[3];
sx q[3];
rz(-1.9347408) q[3];
sx q[3];
rz(-2.5516627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4425519) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(-0.59147269) q[2];
rz(-2.5555723) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(-1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(-2.3024094) q[0];
rz(3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(1.5930088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4913113) q[0];
sx q[0];
rz(-2.6822753) q[0];
sx q[0];
rz(-1.7772872) q[0];
x q[1];
rz(0.71009212) q[2];
sx q[2];
rz(-0.93821628) q[2];
sx q[2];
rz(1.3131504) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3241987) q[1];
sx q[1];
rz(-0.84016582) q[1];
sx q[1];
rz(2.6170931) q[1];
rz(-pi) q[2];
rz(2.5769916) q[3];
sx q[3];
rz(-0.49800107) q[3];
sx q[3];
rz(0.48674395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30535355) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(0.099686064) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(-1.3249741) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56617671) q[0];
sx q[0];
rz(-1.7594936) q[0];
sx q[0];
rz(2.8856522) q[0];
rz(-2.6804965) q[1];
sx q[1];
rz(-2.0979116) q[1];
sx q[1];
rz(0.76006132) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041199112) q[0];
sx q[0];
rz(-1.6738335) q[0];
sx q[0];
rz(-0.11739199) q[0];
rz(-1.1561469) q[2];
sx q[2];
rz(-1.6290602) q[2];
sx q[2];
rz(2.0410048) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6872245) q[1];
sx q[1];
rz(-1.6974653) q[1];
sx q[1];
rz(1.441799) q[1];
x q[2];
rz(2.2171668) q[3];
sx q[3];
rz(-1.8674208) q[3];
sx q[3];
rz(0.38277205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(2.467353) q[2];
rz(-0.21480602) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(3.1242483) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53133416) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(1.1791139) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(-2.0746322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93766312) q[0];
sx q[0];
rz(-1.5852889) q[0];
sx q[0];
rz(3.1209164) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42491575) q[2];
sx q[2];
rz(-1.1669461) q[2];
sx q[2];
rz(2.0770819) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3031591) q[1];
sx q[1];
rz(-1.5232956) q[1];
sx q[1];
rz(2.3179503) q[1];
rz(-pi) q[2];
rz(-0.79640572) q[3];
sx q[3];
rz(-1.4561597) q[3];
sx q[3];
rz(-2.93626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.431488) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(-0.55994326) q[2];
rz(-2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85754919) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(1.42111) q[0];
rz(0.20206085) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(2.2834159) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35408033) q[0];
sx q[0];
rz(-1.6991827) q[0];
sx q[0];
rz(1.7187198) q[0];
x q[1];
rz(-1.6740587) q[2];
sx q[2];
rz(-2.729136) q[2];
sx q[2];
rz(-2.7445284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1865902) q[1];
sx q[1];
rz(-1.5860671) q[1];
sx q[1];
rz(3.112622) q[1];
x q[2];
rz(2.5704727) q[3];
sx q[3];
rz(-1.8652328) q[3];
sx q[3];
rz(-1.7366228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28785607) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(1.3254335) q[2];
rz(2.251513) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(-0.37471399) q[0];
rz(-2.162714) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50586787) q[0];
sx q[0];
rz(-2.3946107) q[0];
sx q[0];
rz(0.58734679) q[0];
x q[1];
rz(-0.71808727) q[2];
sx q[2];
rz(-1.5880843) q[2];
sx q[2];
rz(-1.3912488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6264682) q[1];
sx q[1];
rz(-1.6203468) q[1];
sx q[1];
rz(-1.7527761) q[1];
x q[2];
rz(-0.54578636) q[3];
sx q[3];
rz(-1.3457527) q[3];
sx q[3];
rz(0.033566098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4902041) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(2.6728969) q[2];
rz(1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.7780001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.8796896) q[0];
rz(2.966554) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(-1.5375686) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8309098) q[0];
sx q[0];
rz(-1.791782) q[0];
sx q[0];
rz(1.2587147) q[0];
rz(-1.6749009) q[2];
sx q[2];
rz(-1.9874007) q[2];
sx q[2];
rz(-2.5660851) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9242212) q[1];
sx q[1];
rz(-1.6842168) q[1];
sx q[1];
rz(2.4356615) q[1];
x q[2];
rz(-2.0314625) q[3];
sx q[3];
rz(-1.8970282) q[3];
sx q[3];
rz(0.5512475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1016772) q[2];
sx q[2];
rz(-0.95255178) q[2];
sx q[2];
rz(0.76134479) q[2];
rz(-0.90041655) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(-0.049023978) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3392357) q[0];
sx q[0];
rz(-0.46827066) q[0];
sx q[0];
rz(0.21690579) q[0];
rz(2.5096109) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(-2.1868618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58669421) q[0];
sx q[0];
rz(-2.2838755) q[0];
sx q[0];
rz(-1.9999534) q[0];
rz(-pi) q[1];
rz(0.63072272) q[2];
sx q[2];
rz(-2.4927757) q[2];
sx q[2];
rz(1.1184675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54143822) q[1];
sx q[1];
rz(-0.46240515) q[1];
sx q[1];
rz(2.9701783) q[1];
rz(-pi) q[2];
rz(-2.6984152) q[3];
sx q[3];
rz(-0.9842397) q[3];
sx q[3];
rz(-3.0122258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1251936) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(-0.94474244) q[2];
rz(2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5007297) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(0.8846994) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(-2.3064638) q[2];
sx q[2];
rz(-2.0616812) q[2];
sx q[2];
rz(2.4731935) q[2];
rz(1.9839722) q[3];
sx q[3];
rz(-0.4784085) q[3];
sx q[3];
rz(2.9984409) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
