OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(4.099457) q[0];
sx q[0];
rz(9.2803331) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(-2.1638343) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7841871) q[0];
sx q[0];
rz(-1.4854684) q[0];
sx q[0];
rz(-2.8822495) q[0];
rz(-pi) q[1];
rz(2.833843) q[2];
sx q[2];
rz(-1.0767184) q[2];
sx q[2];
rz(-0.57927629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6482918) q[1];
sx q[1];
rz(-0.19953218) q[1];
sx q[1];
rz(-2.3471911) q[1];
rz(-0.78731491) q[3];
sx q[3];
rz(-1.4437321) q[3];
sx q[3];
rz(0.78391677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4686761) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(-0.93227512) q[2];
rz(0.19876924) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(-0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4345877) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(-0.36112753) q[0];
rz(1.4350285) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(2.3235869) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60048238) q[0];
sx q[0];
rz(-0.4821018) q[0];
sx q[0];
rz(0.62645285) q[0];
rz(-pi) q[1];
rz(-1.1992707) q[2];
sx q[2];
rz(-1.4563592) q[2];
sx q[2];
rz(2.3105846) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4600196) q[1];
sx q[1];
rz(-2.1732554) q[1];
sx q[1];
rz(2.6006992) q[1];
rz(1.2070451) q[3];
sx q[3];
rz(-1.6241637) q[3];
sx q[3];
rz(-1.766891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8890185) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(1.4206295) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(-2.7533598) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(-2.2706568) q[0];
rz(2.8254106) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(0.2972163) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5706144) q[0];
sx q[0];
rz(-1.9750722) q[0];
sx q[0];
rz(-1.6398318) q[0];
x q[1];
rz(2.3702413) q[2];
sx q[2];
rz(-1.0372247) q[2];
sx q[2];
rz(2.9750864) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7797459) q[1];
sx q[1];
rz(-1.6184813) q[1];
sx q[1];
rz(1.7576799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3171595) q[3];
sx q[3];
rz(-1.0472877) q[3];
sx q[3];
rz(1.5200966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3729942) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4363842) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(0.77392459) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(0.68177044) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9457152) q[0];
sx q[0];
rz(-1.8316852) q[0];
sx q[0];
rz(2.9257141) q[0];
rz(-1.2615471) q[2];
sx q[2];
rz(-2.5144858) q[2];
sx q[2];
rz(-0.68324616) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7766429) q[1];
sx q[1];
rz(-0.72307359) q[1];
sx q[1];
rz(0.88100453) q[1];
x q[2];
rz(2.3247129) q[3];
sx q[3];
rz(-0.48135346) q[3];
sx q[3];
rz(-0.64827418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1327847) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(-0.95820367) q[2];
rz(-2.0751674) q[3];
sx q[3];
rz(-1.8582148) q[3];
sx q[3];
rz(-0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(-0.56030309) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(1.6220185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1513838) q[0];
sx q[0];
rz(-0.68400331) q[0];
sx q[0];
rz(0.74599501) q[0];
rz(-pi) q[1];
rz(-1.0137453) q[2];
sx q[2];
rz(-2.1079014) q[2];
sx q[2];
rz(0.86442664) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.8738185) q[1];
sx q[1];
rz(-0.78332892) q[1];
sx q[1];
rz(-2.8591213) q[1];
rz(-pi) q[2];
rz(-2.1836957) q[3];
sx q[3];
rz(-1.6727722) q[3];
sx q[3];
rz(1.0979872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4524298) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(2.9186644) q[2];
rz(-0.034742268) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8054304) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(2.0715332) q[0];
rz(1.7806212) q[1];
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
rz(-2.5491935) q[0];
sx q[0];
rz(-2.1201029) q[0];
sx q[0];
rz(2.4248289) q[0];
rz(-pi) q[1];
rz(0.27988866) q[2];
sx q[2];
rz(-1.9473238) q[2];
sx q[2];
rz(2.1044452) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9575427) q[1];
sx q[1];
rz(-1.8584538) q[1];
sx q[1];
rz(-0.48344739) q[1];
x q[2];
rz(-2.1070126) q[3];
sx q[3];
rz(-1.3048733) q[3];
sx q[3];
rz(1.3724316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.47508919) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(0.63759032) q[2];
rz(0.83667886) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5752983) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(-0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(-0.93820757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095031247) q[0];
sx q[0];
rz(-1.9782269) q[0];
sx q[0];
rz(0.33336063) q[0];
rz(1.0023408) q[2];
sx q[2];
rz(-2.217514) q[2];
sx q[2];
rz(1.2127753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36233703) q[1];
sx q[1];
rz(-1.926683) q[1];
sx q[1];
rz(-0.87902714) q[1];
rz(-pi) q[2];
rz(2.7890117) q[3];
sx q[3];
rz(-0.73871021) q[3];
sx q[3];
rz(1.1383566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2618835) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(1.7555457) q[2];
rz(1.322768) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8742074) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(1.7077131) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(0.83126718) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48152637) q[0];
sx q[0];
rz(-0.43018451) q[0];
sx q[0];
rz(-0.35920401) q[0];
rz(-pi) q[1];
rz(-1.3966884) q[2];
sx q[2];
rz(-0.6422407) q[2];
sx q[2];
rz(-1.0345392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2483406) q[1];
sx q[1];
rz(-2.2447531) q[1];
sx q[1];
rz(-2.2158951) q[1];
rz(-pi) q[2];
rz(1.3475111) q[3];
sx q[3];
rz(-1.7956453) q[3];
sx q[3];
rz(0.68883483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5780118) q[2];
sx q[2];
rz(-1.7687904) q[2];
sx q[2];
rz(0.9643628) q[2];
rz(1.9780805) q[3];
sx q[3];
rz(-0.61146277) q[3];
sx q[3];
rz(2.4826629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37373856) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(-0.97737616) q[0];
rz(1.3865698) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(-2.0358553) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2290105) q[0];
sx q[0];
rz(-1.454103) q[0];
sx q[0];
rz(0.44793655) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90412037) q[2];
sx q[2];
rz(-1.6534001) q[2];
sx q[2];
rz(2.3711575) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.92377201) q[1];
sx q[1];
rz(-0.2982699) q[1];
sx q[1];
rz(-2.24733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2826142) q[3];
sx q[3];
rz(-1.2328706) q[3];
sx q[3];
rz(1.5087138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5332807) q[2];
sx q[2];
rz(-0.38485843) q[2];
sx q[2];
rz(0.14979714) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(-1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6479284) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(-3.0143484) q[0];
rz(-1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(-0.1677992) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084899336) q[0];
sx q[0];
rz(-0.035824422) q[0];
sx q[0];
rz(0.56532677) q[0];
x q[1];
rz(-1.4160412) q[2];
sx q[2];
rz(-2.1733279) q[2];
sx q[2];
rz(-0.11520152) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64393109) q[1];
sx q[1];
rz(-1.5776458) q[1];
sx q[1];
rz(-0.58845206) q[1];
x q[2];
rz(3.0182748) q[3];
sx q[3];
rz(-1.6880369) q[3];
sx q[3];
rz(0.58837147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1148791) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(0.76114571) q[2];
rz(-3.051565) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(-2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5933843) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(-2.3315196) q[2];
sx q[2];
rz(-1.0675061) q[2];
sx q[2];
rz(1.6455417) q[2];
rz(1.0317867) q[3];
sx q[3];
rz(-1.4236593) q[3];
sx q[3];
rz(0.66766213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
