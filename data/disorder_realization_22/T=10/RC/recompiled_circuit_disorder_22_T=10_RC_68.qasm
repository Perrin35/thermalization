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
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(0.97775835) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.190783) q[0];
sx q[0];
rz(-1.3124183) q[0];
sx q[0];
rz(-1.4825312) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0852282) q[2];
sx q[2];
rz(-1.8407514) q[2];
sx q[2];
rz(1.1411238) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.49330083) q[1];
sx q[1];
rz(-0.19953218) q[1];
sx q[1];
rz(2.3471911) q[1];
x q[2];
rz(1.3917189) q[3];
sx q[3];
rz(-0.79154166) q[3];
sx q[3];
rz(0.66034987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4686761) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(2.1762302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4345877) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(-0.36112753) q[0];
rz(-1.4350285) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(2.3235869) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60048238) q[0];
sx q[0];
rz(-2.6594909) q[0];
sx q[0];
rz(-0.62645285) q[0];
rz(-pi) q[1];
rz(-1.9423219) q[2];
sx q[2];
rz(-1.6852334) q[2];
sx q[2];
rz(2.3105846) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4600196) q[1];
sx q[1];
rz(-2.1732554) q[1];
sx q[1];
rz(-0.54089344) q[1];
rz(-pi) q[2];
rz(1.4217671) q[3];
sx q[3];
rz(-2.7741198) q[3];
sx q[3];
rz(3.0847103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8890185) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(-1.7209631) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(-0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7398359) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(2.2706568) q[0];
rz(-2.8254106) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(-0.2972163) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5706144) q[0];
sx q[0];
rz(-1.1665205) q[0];
sx q[0];
rz(-1.5017609) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3702413) q[2];
sx q[2];
rz(-1.0372247) q[2];
sx q[2];
rz(-2.9750864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3618468) q[1];
sx q[1];
rz(-1.6184813) q[1];
sx q[1];
rz(1.7576799) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66629569) q[3];
sx q[3];
rz(-0.94216457) q[3];
sx q[3];
rz(-2.6578238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7685984) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(-0.42759839) q[2];
rz(-1.1887431) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(0.52743131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4363842) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(2.3676681) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(-0.68177044) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6503158) q[0];
sx q[0];
rz(-0.33704764) q[0];
sx q[0];
rz(-2.2469673) q[0];
rz(-pi) q[1];
rz(1.2615471) q[2];
sx q[2];
rz(-0.62710688) q[2];
sx q[2];
rz(2.4583465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75979739) q[1];
sx q[1];
rz(-2.0054381) q[1];
sx q[1];
rz(-2.1684907) q[1];
rz(0.34337266) q[3];
sx q[3];
rz(-1.2265172) q[3];
sx q[3];
rz(-1.6791277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(-0.95820367) q[2];
rz(2.0751674) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(2.7619894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.20345774) q[1];
sx q[1];
rz(1.6220185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1006267) q[0];
sx q[0];
rz(-2.0140411) q[0];
sx q[0];
rz(2.6020781) q[0];
rz(1.0137453) q[2];
sx q[2];
rz(-1.0336913) q[2];
sx q[2];
rz(0.86442664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.49415627) q[1];
sx q[1];
rz(-1.7687706) q[1];
sx q[1];
rz(2.3784749) q[1];
x q[2];
rz(-1.394746) q[3];
sx q[3];
rz(-0.62024833) q[3];
sx q[3];
rz(2.8125417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6891629) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(-2.9186644) q[2];
rz(-0.034742268) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(-0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.3361622) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(1.0700595) q[0];
rz(-1.3609715) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(1.8575352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.703008) q[0];
sx q[0];
rz(-0.87235886) q[0];
sx q[0];
rz(0.75011487) q[0];
x q[1];
rz(2.861704) q[2];
sx q[2];
rz(-1.9473238) q[2];
sx q[2];
rz(1.0371475) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9575427) q[1];
sx q[1];
rz(-1.8584538) q[1];
sx q[1];
rz(-2.6581453) q[1];
x q[2];
rz(2.0606023) q[3];
sx q[3];
rz(-0.59270699) q[3];
sx q[3];
rz(0.61471516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6665035) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(0.63759032) q[2];
rz(-0.83667886) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(2.2033851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0465614) q[0];
sx q[0];
rz(-1.1633658) q[0];
sx q[0];
rz(2.808232) q[0];
rz(-pi) q[1];
rz(0.73056716) q[2];
sx q[2];
rz(-1.1267203) q[2];
sx q[2];
rz(0.72545746) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5349622) q[1];
sx q[1];
rz(-2.3772847) q[1];
sx q[1];
rz(1.0431837) q[1];
rz(-1.8754962) q[3];
sx q[3];
rz(-0.88677553) q[3];
sx q[3];
rz(1.5414433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87970916) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(1.7555457) q[2];
rz(1.8188247) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26738527) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(2.3103255) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2682876) q[0];
sx q[0];
rz(-1.1697066) q[0];
sx q[0];
rz(-1.7307161) q[0];
rz(1.3966884) q[2];
sx q[2];
rz(-2.499352) q[2];
sx q[2];
rz(-1.0345392) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89325209) q[1];
sx q[1];
rz(-0.89683956) q[1];
sx q[1];
rz(0.92569758) q[1];
rz(1.3475111) q[3];
sx q[3];
rz(-1.3459473) q[3];
sx q[3];
rz(-0.68883483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5635809) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-2.1772299) q[2];
rz(1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(-0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.3061378) q[1];
sx q[1];
rz(2.0358553) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7439197) q[0];
sx q[0];
rz(-2.0154675) q[0];
sx q[0];
rz(-1.4414653) q[0];
rz(0.90412037) q[2];
sx q[2];
rz(-1.6534001) q[2];
sx q[2];
rz(-2.3711575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8399664) q[1];
sx q[1];
rz(-1.3857538) q[1];
sx q[1];
rz(-1.3355096) q[1];
x q[2];
rz(0.68007277) q[3];
sx q[3];
rz(-0.44049965) q[3];
sx q[3];
rz(-0.77914731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(2.9917955) q[2];
rz(-1.3730565) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(-1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(3.0143484) q[0];
rz(-1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(-0.1677992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0566933) q[0];
sx q[0];
rz(-0.035824422) q[0];
sx q[0];
rz(2.5762659) q[0];
x q[1];
rz(2.5334353) q[2];
sx q[2];
rz(-1.4434575) q[2];
sx q[2];
rz(-1.5437768) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2192981) q[1];
sx q[1];
rz(-0.98235992) q[1];
sx q[1];
rz(-1.5625619) q[1];
x q[2];
rz(-2.3778902) q[3];
sx q[3];
rz(-2.9716431) q[3];
sx q[3];
rz(-2.9156239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1148791) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(0.76114571) q[2];
rz(-0.090027697) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.54820838) q[0];
sx q[0];
rz(-1.1593288) q[0];
sx q[0];
rz(2.819084) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(2.3315196) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(-0.17100632) q[3];
sx q[3];
rz(-2.1033559) q[3];
sx q[3];
rz(-0.81567473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
