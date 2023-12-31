OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(-3.1414519) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(-1.1934086) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407335) q[0];
sx q[0];
rz(-1.8574323) q[0];
sx q[0];
rz(-0.1851693) q[0];
rz(-pi) q[1];
rz(-0.46618669) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(0.28238645) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57230091) q[1];
sx q[1];
rz(-2.3094059) q[1];
sx q[1];
rz(0.65170793) q[1];
x q[2];
rz(1.3532577) q[3];
sx q[3];
rz(-1.6780403) q[3];
sx q[3];
rz(1.4835964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45941916) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.9127282) q[2];
rz(-1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(1.0128101) q[0];
rz(-3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(2.0181296) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7682122) q[0];
sx q[0];
rz(-2.0313615) q[0];
sx q[0];
rz(0.064231355) q[0];
rz(1.5078817) q[2];
sx q[2];
rz(-2.3496369) q[2];
sx q[2];
rz(-2.6346249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9065735) q[1];
sx q[1];
rz(-2.0779013) q[1];
sx q[1];
rz(2.1706235) q[1];
x q[2];
rz(-0.65621891) q[3];
sx q[3];
rz(-2.0341633) q[3];
sx q[3];
rz(0.046063395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(0.91903764) q[2];
rz(-0.67409003) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27750257) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(0.69349849) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(1.1330053) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61921652) q[0];
sx q[0];
rz(-1.4285061) q[0];
sx q[0];
rz(-0.0033358047) q[0];
rz(2.0793545) q[2];
sx q[2];
rz(-0.8478176) q[2];
sx q[2];
rz(1.6194956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6269835) q[1];
sx q[1];
rz(-1.8585464) q[1];
sx q[1];
rz(2.1327553) q[1];
rz(-pi) q[2];
rz(-0.62882256) q[3];
sx q[3];
rz(-1.0477133) q[3];
sx q[3];
rz(-0.76997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(-1.2934925) q[2];
rz(-0.039316468) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(-1.1608634) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(0.13555759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40977749) q[0];
sx q[0];
rz(-1.8827794) q[0];
sx q[0];
rz(-0.72809763) q[0];
rz(1.3877669) q[2];
sx q[2];
rz(-2.7515538) q[2];
sx q[2];
rz(-0.70034617) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23952661) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(0.86218254) q[1];
rz(-pi) q[2];
rz(-2.0849864) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(1.7900975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23665145) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(0.044163477) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(2.1283545) q[0];
rz(-3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-1.0838881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.097682) q[0];
sx q[0];
rz(-1.327276) q[0];
sx q[0];
rz(2.9527412) q[0];
x q[1];
rz(2.8903557) q[2];
sx q[2];
rz(-1.305797) q[2];
sx q[2];
rz(0.83515177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.52884) q[1];
sx q[1];
rz(-1.4626979) q[1];
sx q[1];
rz(1.0428863) q[1];
rz(-0.089133457) q[3];
sx q[3];
rz(-0.51557487) q[3];
sx q[3];
rz(2.6775132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(0.24442913) q[2];
rz(-0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(-3.0474512) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(2.24618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0953513) q[0];
sx q[0];
rz(-1.6086676) q[0];
sx q[0];
rz(-2.7992159) q[0];
x q[1];
rz(-2.5080639) q[2];
sx q[2];
rz(-1.3760566) q[2];
sx q[2];
rz(0.63846904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0370334) q[1];
sx q[1];
rz(-1.7393141) q[1];
sx q[1];
rz(0.044635459) q[1];
rz(-pi) q[2];
x q[2];
rz(1.496109) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(-2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.133693) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(0.80319476) q[2];
rz(-1.9512272) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-2.7289594) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(-2.6224526) q[0];
rz(0.58147645) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(1.8849467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6709267) q[0];
sx q[0];
rz(-1.5304655) q[0];
sx q[0];
rz(1.8223902) q[0];
rz(-pi) q[1];
rz(-2.0357382) q[2];
sx q[2];
rz(-1.3824116) q[2];
sx q[2];
rz(0.18907324) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7278442) q[1];
sx q[1];
rz(-1.8520978) q[1];
sx q[1];
rz(-1.6378228) q[1];
rz(-pi) q[2];
rz(-1.2506966) q[3];
sx q[3];
rz(-0.88722908) q[3];
sx q[3];
rz(-0.86910955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1533623) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(-2.2310232) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(1.6114657) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751188) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(0.30817729) q[0];
rz(-3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(0.38696188) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68650866) q[0];
sx q[0];
rz(-2.5065656) q[0];
sx q[0];
rz(1.1839266) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1808124) q[2];
sx q[2];
rz(-2.3443601) q[2];
sx q[2];
rz(-2.0955992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4657198) q[1];
sx q[1];
rz(-0.67589251) q[1];
sx q[1];
rz(-0.11996108) q[1];
rz(-pi) q[2];
x q[2];
rz(1.79027) q[3];
sx q[3];
rz(-2.05728) q[3];
sx q[3];
rz(-2.6284077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(0.44000885) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(2.0429042) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(0.049302014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44137529) q[0];
sx q[0];
rz(-0.8964296) q[0];
sx q[0];
rz(-2.3250439) q[0];
rz(-pi) q[1];
rz(2.2531613) q[2];
sx q[2];
rz(-1.9677791) q[2];
sx q[2];
rz(1.9156485) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.8763435) q[1];
sx q[1];
rz(-2.0298404) q[1];
sx q[1];
rz(2.6190119) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4008972) q[3];
sx q[3];
rz(-1.9441838) q[3];
sx q[3];
rz(-1.5965366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0631642) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(-2.4196529) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(-2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9938875) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(2.06185) q[0];
rz(1.059277) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(1.4019029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7584383) q[0];
sx q[0];
rz(-2.8935195) q[0];
sx q[0];
rz(1.0323348) q[0];
rz(-pi) q[1];
rz(-0.14342587) q[2];
sx q[2];
rz(-0.18866814) q[2];
sx q[2];
rz(0.27486899) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55587308) q[1];
sx q[1];
rz(-1.266529) q[1];
sx q[1];
rz(0.85822206) q[1];
rz(-pi) q[2];
rz(-0.28016443) q[3];
sx q[3];
rz(-0.6503085) q[3];
sx q[3];
rz(-3.0349817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(-2.5907717) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.14810066) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(0.91611721) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(-0.79402906) q[2];
sx q[2];
rz(-0.90894884) q[2];
sx q[2];
rz(-0.85469645) q[2];
rz(-0.017833088) q[3];
sx q[3];
rz(-2.087839) q[3];
sx q[3];
rz(2.9057333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
