OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(5.1685652) q[0];
sx q[0];
rz(9.4246372) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(4.1058022) q[1];
sx q[1];
rz(10.618187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50085917) q[0];
sx q[0];
rz(-1.2841604) q[0];
sx q[0];
rz(0.1851693) q[0];
rz(-pi) q[1];
rz(-2.5932181) q[2];
sx q[2];
rz(-1.3142685) q[2];
sx q[2];
rz(-2.246849) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6686033) q[1];
sx q[1];
rz(-2.0358634) q[1];
sx q[1];
rz(-2.4238062) q[1];
rz(-3.0317806) q[3];
sx q[3];
rz(-1.7870652) q[3];
sx q[3];
rz(3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6821735) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(1.4131644) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(-1.0128101) q[0];
rz(0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-1.123463) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2259953) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(-2.0321839) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6337109) q[2];
sx q[2];
rz(-0.79195576) q[2];
sx q[2];
rz(2.6346249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1264798) q[1];
sx q[1];
rz(-2.0868595) q[1];
sx q[1];
rz(2.5491787) q[1];
x q[2];
rz(-0.68617679) q[3];
sx q[3];
rz(-0.78305972) q[3];
sx q[3];
rz(0.99883294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79364395) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(-0.91903764) q[2];
rz(-2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(1.526171) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640901) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(0.69349849) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(-2.0085874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61921652) q[0];
sx q[0];
rz(-1.7130865) q[0];
sx q[0];
rz(-3.1382568) q[0];
rz(-pi) q[1];
rz(-1.0622382) q[2];
sx q[2];
rz(-0.8478176) q[2];
sx q[2];
rz(1.6194956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.36724597) q[1];
sx q[1];
rz(-0.62421747) q[1];
sx q[1];
rz(1.063785) q[1];
rz(-pi) q[2];
rz(2.3660925) q[3];
sx q[3];
rz(-0.79458487) q[3];
sx q[3];
rz(0.19898181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(3.1022762) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.8815276) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.9807293) q[0];
rz(2.2456031) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(0.13555759) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4281222) q[0];
sx q[0];
rz(-2.2566416) q[0];
sx q[0];
rz(1.1629521) q[0];
x q[1];
rz(1.9549471) q[2];
sx q[2];
rz(-1.5015366) q[2];
sx q[2];
rz(-2.4406976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23952661) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(-0.86218254) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8454148) q[3];
sx q[3];
rz(-1.7227968) q[3];
sx q[3];
rz(2.8677468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9049412) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(-0.87990749) q[2];
rz(-3.0974292) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039466) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(0.049731072) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(-2.0577046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5729382) q[0];
sx q[0];
rz(-1.7540115) q[0];
sx q[0];
rz(1.3230447) q[0];
x q[1];
rz(1.8439699) q[2];
sx q[2];
rz(-1.3285085) q[2];
sx q[2];
rz(-0.66852409) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6127527) q[1];
sx q[1];
rz(-1.6788947) q[1];
sx q[1];
rz(-2.0987064) q[1];
rz(-pi) q[2];
rz(1.5203939) q[3];
sx q[3];
rz(-2.084123) q[3];
sx q[3];
rz(-0.36171519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9087387) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(0.50306815) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(-2.24618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0462414) q[0];
sx q[0];
rz(-1.5329251) q[0];
sx q[0];
rz(-0.34237679) q[0];
rz(-pi) q[1];
rz(1.8108098) q[2];
sx q[2];
rz(-0.95108205) q[2];
sx q[2];
rz(1.0735219) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0370334) q[1];
sx q[1];
rz(-1.7393141) q[1];
sx q[1];
rz(-0.044635459) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5395245) q[3];
sx q[3];
rz(-3.0668689) q[3];
sx q[3];
rz(0.63262313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(-1.1903654) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(-0.58147645) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(1.8849467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4706659) q[0];
sx q[0];
rz(-1.5304655) q[0];
sx q[0];
rz(1.3192024) q[0];
rz(-pi) q[1];
rz(1.9728327) q[2];
sx q[2];
rz(-0.49905825) q[2];
sx q[2];
rz(1.7390342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1756806) q[1];
sx q[1];
rz(-1.5064081) q[1];
sx q[1];
rz(-2.8596911) q[1];
rz(-pi) q[2];
rz(-1.2506966) q[3];
sx q[3];
rz(-0.88722908) q[3];
sx q[3];
rz(-0.86910955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1533623) q[2];
sx q[2];
rz(-1.0299269) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0664739) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(-0.30817729) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(-0.38696188) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9403801) q[0];
sx q[0];
rz(-1.3450755) q[0];
sx q[0];
rz(2.1696027) q[0];
rz(2.6110821) q[2];
sx q[2];
rz(-2.1973917) q[2];
sx q[2];
rz(-0.26041398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-2.4657001) q[1];
sx q[1];
rz(-0.11996108) q[1];
rz(-pi) q[2];
x q[2];
rz(1.79027) q[3];
sx q[3];
rz(-1.0843127) q[3];
sx q[3];
rz(2.6284077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(2.7015838) q[2];
rz(-0.7157588) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(2.0429042) q[0];
rz(2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54287275) q[0];
sx q[0];
rz(-2.1763986) q[0];
sx q[0];
rz(-0.70830317) q[0];
rz(-pi) q[1];
rz(0.88843139) q[2];
sx q[2];
rz(-1.9677791) q[2];
sx q[2];
rz(1.2259442) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1028565) q[1];
sx q[1];
rz(-0.68118459) q[1];
sx q[1];
rz(-0.78050651) q[1];
rz(2.6155869) q[3];
sx q[3];
rz(-2.3283539) q[3];
sx q[3];
rz(0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(-0.72193974) q[2];
rz(-2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(2.06185) q[0];
rz(-2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(1.7396897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7125177) q[0];
sx q[0];
rz(-1.4445462) q[0];
sx q[0];
rz(1.35668) q[0];
rz(-pi) q[1];
rz(-0.18677588) q[2];
sx q[2];
rz(-1.5439856) q[2];
sx q[2];
rz(-1.1550127) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8733752) q[1];
sx q[1];
rz(-0.89726071) q[1];
sx q[1];
rz(2.7482277) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8614282) q[3];
sx q[3];
rz(-0.6503085) q[3];
sx q[3];
rz(-0.10661099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(-1.6213017) q[2];
rz(2.5907717) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(-2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14810066) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(-0.91611721) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(-0.73268391) q[2];
sx q[2];
rz(-0.97326836) q[2];
sx q[2];
rz(1.2748981) q[2];
rz(-3.1237596) q[3];
sx q[3];
rz(-1.0537536) q[3];
sx q[3];
rz(-0.23585933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
