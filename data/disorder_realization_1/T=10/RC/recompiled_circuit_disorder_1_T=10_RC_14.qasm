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
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50085917) q[0];
sx q[0];
rz(-1.8574323) q[0];
sx q[0];
rz(0.1851693) q[0];
x q[1];
rz(-2.675406) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(-0.28238645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8661583) q[1];
sx q[1];
rz(-2.1992116) q[1];
sx q[1];
rz(0.98316146) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0317806) q[3];
sx q[3];
rz(-1.3545274) q[3];
sx q[3];
rz(-3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45941916) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(1.2288644) q[2];
rz(1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(-1.123463) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9119308) q[0];
sx q[0];
rz(-0.46470416) q[0];
sx q[0];
rz(1.4421411) q[0];
rz(-2.3617619) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(1.1080527) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23501913) q[1];
sx q[1];
rz(-1.0636914) q[1];
sx q[1];
rz(2.1706235) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65621891) q[3];
sx q[3];
rz(-2.0341633) q[3];
sx q[3];
rz(-0.046063395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79364395) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(2.222555) q[2];
rz(-2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640901) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(2.0085874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1904859) q[0];
sx q[0];
rz(-1.5740984) q[0];
sx q[0];
rz(-1.4285054) q[0];
rz(2.6373367) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(-0.91932636) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7743467) q[1];
sx q[1];
rz(-2.5173752) q[1];
sx q[1];
rz(2.0778076) q[1];
rz(2.5127701) q[3];
sx q[3];
rz(-1.0477133) q[3];
sx q[3];
rz(-0.76997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2514078) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.2599729) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(-1.1608634) q[0];
rz(2.2456031) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(3.0060351) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7134705) q[0];
sx q[0];
rz(-2.2566416) q[0];
sx q[0];
rz(1.1629521) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1866456) q[2];
sx q[2];
rz(-1.5015366) q[2];
sx q[2];
rz(-0.70089507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.23952661) q[1];
sx q[1];
rz(-2.8956928) q[1];
sx q[1];
rz(0.86218254) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0849864) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(1.3514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23665145) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(2.2616852) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-2.0577046) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7679374) q[0];
sx q[0];
rz(-0.3070139) q[0];
sx q[0];
rz(-0.92371773) q[0];
rz(-pi) q[1];
rz(-1.8439699) q[2];
sx q[2];
rz(-1.8130842) q[2];
sx q[2];
rz(-0.66852409) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6127527) q[1];
sx q[1];
rz(-1.6788947) q[1];
sx q[1];
rz(-1.0428863) q[1];
rz(-pi) q[2];
rz(-0.089133457) q[3];
sx q[3];
rz(-0.51557487) q[3];
sx q[3];
rz(-0.46407947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9087387) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(0.24442913) q[2];
rz(2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(0.89541268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46195128) q[0];
sx q[0];
rz(-1.9129176) q[0];
sx q[0];
rz(1.530594) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8199811) q[2];
sx q[2];
rz(-2.4827637) q[2];
sx q[2];
rz(-2.4668601) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9858866) q[1];
sx q[1];
rz(-2.967318) q[1];
sx q[1];
rz(1.3143015) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6454837) q[3];
sx q[3];
rz(-1.5731305) q[3];
sx q[3];
rz(2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0078997) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(0.80319476) q[2];
rz(-1.1903654) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(-2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(-2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(-1.8849467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6709267) q[0];
sx q[0];
rz(-1.5304655) q[0];
sx q[0];
rz(1.8223902) q[0];
rz(-pi) q[1];
rz(2.9314552) q[2];
sx q[2];
rz(-1.1147095) q[2];
sx q[2];
rz(1.2880585) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9650967) q[1];
sx q[1];
rz(-0.28897044) q[1];
sx q[1];
rz(-0.22775905) q[1];
x q[2];
rz(-2.4323746) q[3];
sx q[3];
rz(-1.3243444) q[3];
sx q[3];
rz(-0.49530464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1533623) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.3640277) q[2];
rz(-0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664739) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(0.30817729) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(2.7546308) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2012126) q[0];
sx q[0];
rz(-1.7965172) q[0];
sx q[0];
rz(2.1696027) q[0];
x q[1];
rz(0.53051051) q[2];
sx q[2];
rz(-0.94420099) q[2];
sx q[2];
rz(-0.26041398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3124233) q[1];
sx q[1];
rz(-2.2409391) q[1];
sx q[1];
rz(-1.6664684) q[1];
rz(-pi) q[2];
rz(-0.49658637) q[3];
sx q[3];
rz(-1.3771309) q[3];
sx q[3];
rz(0.95369875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62362921) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(2.7015838) q[2];
rz(-0.7157588) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14116645) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(-1.0986885) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-0.049302014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7002174) q[0];
sx q[0];
rz(-0.8964296) q[0];
sx q[0];
rz(-0.81654878) q[0];
x q[1];
rz(-0.98408913) q[2];
sx q[2];
rz(-0.77312914) q[2];
sx q[2];
rz(-3.0423321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2652492) q[1];
sx q[1];
rz(-1.1117522) q[1];
sx q[1];
rz(2.6190119) q[1];
x q[2];
rz(-2.4008972) q[3];
sx q[3];
rz(-1.9441838) q[3];
sx q[3];
rz(-1.5965366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14770517) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(2.06185) q[0];
rz(2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.7396897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3106874) q[0];
sx q[0];
rz(-1.3584104) q[0];
sx q[0];
rz(3.0124245) q[0];
x q[1];
rz(-0.14342587) q[2];
sx q[2];
rz(-2.9529245) q[2];
sx q[2];
rz(2.8667237) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.55587308) q[1];
sx q[1];
rz(-1.8750637) q[1];
sx q[1];
rz(2.2833706) q[1];
rz(-0.28016443) q[3];
sx q[3];
rz(-2.4912842) q[3];
sx q[3];
rz(3.0349817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4828651) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(-0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.993492) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(0.91611721) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(2.312071) q[2];
sx q[2];
rz(-0.98486949) q[2];
sx q[2];
rz(-2.9688901) q[2];
rz(1.5394474) q[3];
sx q[3];
rz(-0.51732224) q[3];
sx q[3];
rz(2.8696685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
