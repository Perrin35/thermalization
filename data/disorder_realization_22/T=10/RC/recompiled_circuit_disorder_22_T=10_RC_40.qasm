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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.30774967) q[2];
sx q[2];
rz(-1.0767184) q[2];
sx q[2];
rz(2.5623164) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2797151) q[1];
sx q[1];
rz(-1.7126843) q[1];
sx q[1];
rz(-3.0008297) q[1];
rz(2.3542777) q[3];
sx q[3];
rz(-1.6978605) q[3];
sx q[3];
rz(-0.78391677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4686761) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(-0.19876924) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4345877) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(-2.7804651) q[0];
rz(1.7065642) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(2.3235869) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6011071) q[0];
sx q[0];
rz(-1.2955106) q[0];
sx q[0];
rz(2.7406373) q[0];
rz(-pi) q[1];
rz(-1.9423219) q[2];
sx q[2];
rz(-1.6852334) q[2];
sx q[2];
rz(2.3105846) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9243014) q[1];
sx q[1];
rz(-1.1326619) q[1];
sx q[1];
rz(0.89466722) q[1];
rz(-pi) q[2];
rz(-1.2070451) q[3];
sx q[3];
rz(-1.6241637) q[3];
sx q[3];
rz(-1.3747017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8890185) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(1.7209631) q[2];
rz(-1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(-0.38823286) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7398359) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(0.87093583) q[0];
rz(2.8254106) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(2.8443764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0270099) q[0];
sx q[0];
rz(-1.5073338) q[0];
sx q[0];
rz(0.40513904) q[0];
rz(-pi) q[1];
rz(2.3702413) q[2];
sx q[2];
rz(-1.0372247) q[2];
sx q[2];
rz(-0.16650621) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45589046) q[1];
sx q[1];
rz(-2.9487902) q[1];
sx q[1];
rz(-1.3193858) q[1];
x q[2];
rz(0.8244332) q[3];
sx q[3];
rz(-2.094305) q[3];
sx q[3];
rz(1.5200966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3729942) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(2.7139943) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(-0.52743131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7052085) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(2.3676681) q[0];
rz(-0.71290839) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(2.4598222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43142372) q[0];
sx q[0];
rz(-1.7792601) q[0];
sx q[0];
rz(1.8375977) q[0];
x q[1];
rz(-2.1749928) q[2];
sx q[2];
rz(-1.3912429) q[2];
sx q[2];
rz(-2.0008848) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75979739) q[1];
sx q[1];
rz(-1.1361546) q[1];
sx q[1];
rz(0.97310193) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34337266) q[3];
sx q[3];
rz(-1.2265172) q[3];
sx q[3];
rz(1.6791277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(-2.183389) q[2];
rz(1.0664252) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-1.6220185) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9902089) q[0];
sx q[0];
rz(-0.68400331) q[0];
sx q[0];
rz(-0.74599501) q[0];
rz(-2.5298169) q[2];
sx q[2];
rz(-2.0423186) q[2];
sx q[2];
rz(0.39786354) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8738185) q[1];
sx q[1];
rz(-0.78332892) q[1];
sx q[1];
rz(-2.8591213) q[1];
rz(-3.0171379) q[3];
sx q[3];
rz(-2.180047) q[3];
sx q[3];
rz(2.5973158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6891629) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(-0.2229283) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(-0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3361622) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(1.0700595) q[0];
rz(-1.7806212) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(-1.8575352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7362471) q[0];
sx q[0];
rz(-0.9760455) q[0];
sx q[0];
rz(2.2527762) q[0];
rz(-0.27988866) q[2];
sx q[2];
rz(-1.1942689) q[2];
sx q[2];
rz(2.1044452) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2594086) q[1];
sx q[1];
rz(-2.5849197) q[1];
sx q[1];
rz(2.5747719) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8347557) q[3];
sx q[3];
rz(-1.0553428) q[3];
sx q[3];
rz(-0.043434867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47508919) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(-2.5040023) q[2];
rz(-2.3049138) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(-1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56629431) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(2.7138846) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(0.93820757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81323775) q[0];
sx q[0];
rz(-0.52044808) q[0];
sx q[0];
rz(0.92207272) q[0];
rz(2.5222048) q[2];
sx q[2];
rz(-2.3084547) q[2];
sx q[2];
rz(2.7433861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7792556) q[1];
sx q[1];
rz(-1.926683) q[1];
sx q[1];
rz(-0.87902714) q[1];
rz(-pi) q[2];
rz(-2.7890117) q[3];
sx q[3];
rz(-0.73871021) q[3];
sx q[3];
rz(-1.1383566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.87970916) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(1.3860469) q[2];
rz(1.8188247) q[3];
sx q[3];
rz(-1.991792) q[3];
sx q[3];
rz(0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8742074) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(-1.7077131) q[0];
rz(-1.8677615) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(-2.3103255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76038218) q[0];
sx q[0];
rz(-1.7179278) q[0];
sx q[0];
rz(-0.40572625) q[0];
x q[1];
rz(2.2057461) q[2];
sx q[2];
rz(-1.6747464) q[2];
sx q[2];
rz(2.7452591) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89325209) q[1];
sx q[1];
rz(-2.2447531) q[1];
sx q[1];
rz(-2.2158951) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7940815) q[3];
sx q[3];
rz(-1.7956453) q[3];
sx q[3];
rz(-2.4527578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5635809) q[2];
sx q[2];
rz(-1.7687904) q[2];
sx q[2];
rz(0.9643628) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(0.65892974) q[3];
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
rz(-0.37373856) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(-2.1642165) q[0];
rz(-1.7550229) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(-1.1057373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9125821) q[0];
sx q[0];
rz(-1.454103) q[0];
sx q[0];
rz(0.44793655) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90412037) q[2];
sx q[2];
rz(-1.4881926) q[2];
sx q[2];
rz(0.77043515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9165009) q[1];
sx q[1];
rz(-1.3396001) q[1];
sx q[1];
rz(-2.9514312) q[1];
rz(-2.4615199) q[3];
sx q[3];
rz(-2.701093) q[3];
sx q[3];
rz(0.77914731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5332807) q[2];
sx q[2];
rz(-0.38485843) q[2];
sx q[2];
rz(-2.9917955) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6479284) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(-0.1272442) q[0];
rz(-1.6607025) q[1];
sx q[1];
rz(-0.36247411) q[1];
sx q[1];
rz(-2.9737934) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6608749) q[0];
sx q[0];
rz(-1.6010451) q[0];
sx q[0];
rz(-1.5515996) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9211505) q[2];
sx q[2];
rz(-0.61969212) q[2];
sx q[2];
rz(-2.9881791) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4976616) q[1];
sx q[1];
rz(-1.5639468) q[1];
sx q[1];
rz(2.5531406) q[1];
x q[2];
rz(-0.12331788) q[3];
sx q[3];
rz(-1.6880369) q[3];
sx q[3];
rz(-2.5532212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1148791) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(-2.3804469) q[2];
rz(-3.051565) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
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
rz(0.38800115) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(-0.81007304) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(0.17100632) q[3];
sx q[3];
rz(-1.0382367) q[3];
sx q[3];
rz(2.3259179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
