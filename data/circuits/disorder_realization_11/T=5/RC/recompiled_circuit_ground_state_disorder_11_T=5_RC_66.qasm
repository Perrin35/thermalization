OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.017555822) q[0];
sx q[0];
rz(-2.9889844) q[0];
sx q[0];
rz(0.18327644) q[0];
rz(-2.2493989) q[1];
sx q[1];
rz(-2.0018938) q[1];
sx q[1];
rz(0.5782063) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86872549) q[0];
sx q[0];
rz(-1.5807165) q[0];
sx q[0];
rz(0.091477576) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.033138795) q[2];
sx q[2];
rz(-1.6553179) q[2];
sx q[2];
rz(3.0042841) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7664585) q[1];
sx q[1];
rz(-0.69978461) q[1];
sx q[1];
rz(-2.3699939) q[1];
x q[2];
rz(3.105025) q[3];
sx q[3];
rz(-2.0362377) q[3];
sx q[3];
rz(-2.5643666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9649428) q[2];
sx q[2];
rz(-2.2735333) q[2];
sx q[2];
rz(-3.0435666) q[2];
rz(0.81884223) q[3];
sx q[3];
rz(-3.0374073) q[3];
sx q[3];
rz(0.63481832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1408511) q[0];
sx q[0];
rz(-0.31743693) q[0];
sx q[0];
rz(2.4733518) q[0];
rz(2.5375598) q[1];
sx q[1];
rz(-3.1112473) q[1];
sx q[1];
rz(0.86413962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0111496) q[0];
sx q[0];
rz(-0.98447039) q[0];
sx q[0];
rz(-1.4545813) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7145442) q[2];
sx q[2];
rz(-0.93736726) q[2];
sx q[2];
rz(-2.065393) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.171794) q[1];
sx q[1];
rz(-1.8117597) q[1];
sx q[1];
rz(-1.765224) q[1];
rz(-0.92738028) q[3];
sx q[3];
rz(-2.815554) q[3];
sx q[3];
rz(3.0969574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16022564) q[2];
sx q[2];
rz(-1.9619433) q[2];
sx q[2];
rz(0.78440624) q[2];
rz(-1.2330327) q[3];
sx q[3];
rz(-0.27850702) q[3];
sx q[3];
rz(-1.2109141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8104372) q[0];
sx q[0];
rz(-1.5758608) q[0];
sx q[0];
rz(1.021215) q[0];
rz(-1.5039697) q[1];
sx q[1];
rz(-0.31871381) q[1];
sx q[1];
rz(2.5655897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9224145) q[0];
sx q[0];
rz(-1.6873345) q[0];
sx q[0];
rz(-0.4187421) q[0];
x q[1];
rz(1.1812649) q[2];
sx q[2];
rz(-1.3329055) q[2];
sx q[2];
rz(-2.212461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29167029) q[1];
sx q[1];
rz(-1.2322767) q[1];
sx q[1];
rz(0.68847076) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3542451) q[3];
sx q[3];
rz(-1.8365905) q[3];
sx q[3];
rz(0.98255052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7716498) q[2];
sx q[2];
rz(-0.45054951) q[2];
sx q[2];
rz(-0.096262781) q[2];
rz(-0.43109584) q[3];
sx q[3];
rz(-2.1754706) q[3];
sx q[3];
rz(-0.14804429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8226606) q[0];
sx q[0];
rz(-3.0156101) q[0];
sx q[0];
rz(-2.8462963) q[0];
rz(0.88116208) q[1];
sx q[1];
rz(-2.9144139) q[1];
sx q[1];
rz(-2.6491902) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32689708) q[0];
sx q[0];
rz(-1.7806223) q[0];
sx q[0];
rz(-2.9710415) q[0];
x q[1];
rz(3.0922331) q[2];
sx q[2];
rz(-1.4774895) q[2];
sx q[2];
rz(0.24053426) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4223692) q[1];
sx q[1];
rz(-1.2417267) q[1];
sx q[1];
rz(0.13767875) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2232699) q[3];
sx q[3];
rz(-1.1545968) q[3];
sx q[3];
rz(2.2496441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1262576) q[2];
sx q[2];
rz(-2.8110562) q[2];
sx q[2];
rz(0.34273657) q[2];
rz(-2.1526509) q[3];
sx q[3];
rz(-1.1611232) q[3];
sx q[3];
rz(-0.69800085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2547176) q[0];
sx q[0];
rz(-2.8509792) q[0];
sx q[0];
rz(-1.2683723) q[0];
rz(0.36240029) q[1];
sx q[1];
rz(-2.2762894) q[1];
sx q[1];
rz(1.5579582) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1311296) q[0];
sx q[0];
rz(-2.4005425) q[0];
sx q[0];
rz(2.4416832) q[0];
rz(1.8475632) q[2];
sx q[2];
rz(-1.74969) q[2];
sx q[2];
rz(-1.8362294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4851361) q[1];
sx q[1];
rz(-0.97378661) q[1];
sx q[1];
rz(-1.5890501) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5823808) q[3];
sx q[3];
rz(-2.5413168) q[3];
sx q[3];
rz(1.6299316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4871939) q[2];
sx q[2];
rz(-2.0687658) q[2];
sx q[2];
rz(0.85713345) q[2];
rz(-2.1060139) q[3];
sx q[3];
rz(-2.3643957) q[3];
sx q[3];
rz(-0.63151675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.421748) q[0];
sx q[0];
rz(-2.6578465) q[0];
sx q[0];
rz(-0.33472043) q[0];
rz(2.1597629) q[1];
sx q[1];
rz(-1.5900759) q[1];
sx q[1];
rz(2.2614711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28583254) q[0];
sx q[0];
rz(-1.8632392) q[0];
sx q[0];
rz(-2.9797735) q[0];
rz(0.99079972) q[2];
sx q[2];
rz(-2.8282165) q[2];
sx q[2];
rz(-2.6277781) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0047061027) q[1];
sx q[1];
rz(-2.0396173) q[1];
sx q[1];
rz(-0.15699082) q[1];
rz(-pi) q[2];
rz(0.37161785) q[3];
sx q[3];
rz(-1.1303969) q[3];
sx q[3];
rz(-0.49345106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6191787) q[2];
sx q[2];
rz(-0.48888561) q[2];
sx q[2];
rz(1.4948814) q[2];
rz(-1.8374247) q[3];
sx q[3];
rz(-2.2924278) q[3];
sx q[3];
rz(-2.9554534) q[3];
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
rz(-2.3592767) q[0];
sx q[0];
rz(-2.9501811) q[0];
sx q[0];
rz(3.0707448) q[0];
rz(0.62295667) q[1];
sx q[1];
rz(-2.7406335) q[1];
sx q[1];
rz(-2.7080022) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1529652) q[0];
sx q[0];
rz(-0.30634964) q[0];
sx q[0];
rz(-3.0393837) q[0];
x q[1];
rz(-0.18871267) q[2];
sx q[2];
rz(-2.8450539) q[2];
sx q[2];
rz(1.485757) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5035491) q[1];
sx q[1];
rz(-1.2473599) q[1];
sx q[1];
rz(-1.4694609) q[1];
rz(-pi) q[2];
rz(2.6979792) q[3];
sx q[3];
rz(-1.6742236) q[3];
sx q[3];
rz(-0.87956968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5101461) q[2];
sx q[2];
rz(-0.53312174) q[2];
sx q[2];
rz(2.5041637) q[2];
rz(1.140945) q[3];
sx q[3];
rz(-2.7009522) q[3];
sx q[3];
rz(-0.91658896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3797853) q[0];
sx q[0];
rz(-2.8768235) q[0];
sx q[0];
rz(2.8763212) q[0];
rz(3.1130262) q[1];
sx q[1];
rz(-0.18762372) q[1];
sx q[1];
rz(0.92715895) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4521479) q[0];
sx q[0];
rz(-1.5981705) q[0];
sx q[0];
rz(0.22449799) q[0];
rz(-pi) q[1];
rz(0.36225944) q[2];
sx q[2];
rz(-0.39316828) q[2];
sx q[2];
rz(2.071381) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3520364) q[1];
sx q[1];
rz(-1.6733783) q[1];
sx q[1];
rz(-0.7898777) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42334564) q[3];
sx q[3];
rz(-1.4676499) q[3];
sx q[3];
rz(-2.1684017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.022543) q[2];
sx q[2];
rz(-1.0485342) q[2];
sx q[2];
rz(1.0724462) q[2];
rz(-0.33506814) q[3];
sx q[3];
rz(-3.0266893) q[3];
sx q[3];
rz(-2.9160299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-0.91113126) q[0];
sx q[0];
rz(-1.9879531) q[0];
sx q[0];
rz(-2.352584) q[0];
rz(0.59793961) q[1];
sx q[1];
rz(-1.1006678) q[1];
sx q[1];
rz(-0.71802968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5563468) q[0];
sx q[0];
rz(-1.6119154) q[0];
sx q[0];
rz(-1.4151735) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4903551) q[2];
sx q[2];
rz(-0.1165963) q[2];
sx q[2];
rz(1.9060022) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78290547) q[1];
sx q[1];
rz(-2.1456239) q[1];
sx q[1];
rz(2.3692822) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6634462) q[3];
sx q[3];
rz(-1.7481553) q[3];
sx q[3];
rz(-0.58045829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52704063) q[2];
sx q[2];
rz(-2.8895832) q[2];
sx q[2];
rz(1.4177812) q[2];
rz(-1.1082015) q[3];
sx q[3];
rz(-0.36644822) q[3];
sx q[3];
rz(-0.38780701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80477667) q[0];
sx q[0];
rz(-2.5122061) q[0];
sx q[0];
rz(2.8861073) q[0];
rz(0.77087036) q[1];
sx q[1];
rz(-1.5893385) q[1];
sx q[1];
rz(-1.593387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6219538) q[0];
sx q[0];
rz(-2.0072091) q[0];
sx q[0];
rz(-3.0163517) q[0];
rz(0.22996567) q[2];
sx q[2];
rz(-0.88092025) q[2];
sx q[2];
rz(1.3892737) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.080049971) q[1];
sx q[1];
rz(-1.1876186) q[1];
sx q[1];
rz(1.8661168) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53867619) q[3];
sx q[3];
rz(-2.4823501) q[3];
sx q[3];
rz(1.074312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8678681) q[2];
sx q[2];
rz(-0.74213433) q[2];
sx q[2];
rz(1.254427) q[2];
rz(-2.6015688) q[3];
sx q[3];
rz(-3.0906782) q[3];
sx q[3];
rz(2.36256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.222432) q[0];
sx q[0];
rz(-1.5800911) q[0];
sx q[0];
rz(2.0240361) q[0];
rz(0.22668214) q[1];
sx q[1];
rz(-3.0381028) q[1];
sx q[1];
rz(1.4729952) q[1];
rz(0.99074715) q[2];
sx q[2];
rz(-1.2540091) q[2];
sx q[2];
rz(-2.2768496) q[2];
rz(-2.4189714) q[3];
sx q[3];
rz(-1.327805) q[3];
sx q[3];
rz(-1.1205825) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
