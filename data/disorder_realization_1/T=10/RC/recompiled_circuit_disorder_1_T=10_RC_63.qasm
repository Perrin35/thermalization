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
rz(-1.1146201) q[0];
sx q[0];
rz(-0.00014076509) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1228468) q[0];
sx q[0];
rz(-1.393264) q[0];
sx q[0];
rz(1.8621423) q[0];
rz(-pi) q[1];
rz(-1.2725856) q[2];
sx q[2];
rz(-2.0993005) q[2];
sx q[2];
rz(-0.82982153) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5692917) q[1];
sx q[1];
rz(-2.3094059) q[1];
sx q[1];
rz(2.4898847) q[1];
rz(1.1080997) q[3];
sx q[3];
rz(-2.8994312) q[3];
sx q[3];
rz(0.36377454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6821735) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.9127282) q[2];
rz(1.4131644) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(-1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(-2.1287825) q[0];
rz(-3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(1.123463) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2259953) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(-2.0321839) q[0];
x q[1];
rz(2.3617619) q[2];
sx q[2];
rz(-1.6155598) q[2];
sx q[2];
rz(-2.0335399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.23501913) q[1];
sx q[1];
rz(-2.0779013) q[1];
sx q[1];
rz(0.97096918) q[1];
rz(-2.4554159) q[3];
sx q[3];
rz(-0.78305972) q[3];
sx q[3];
rz(-0.99883294) q[3];
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
rz(2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(-1.526171) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(1.1330053) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1904859) q[0];
sx q[0];
rz(-1.5740984) q[0];
sx q[0];
rz(1.4285054) q[0];
x q[1];
rz(-2.6373367) q[2];
sx q[2];
rz(-0.85668737) q[2];
sx q[2];
rz(2.2222663) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36724597) q[1];
sx q[1];
rz(-0.62421747) q[1];
sx q[1];
rz(-2.0778076) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77550019) q[3];
sx q[3];
rz(-0.79458487) q[3];
sx q[3];
rz(0.19898181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8901849) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816198) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(1.1608634) q[0];
rz(2.2456031) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(-0.13555759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318152) q[0];
sx q[0];
rz(-1.8827794) q[0];
sx q[0];
rz(-0.72809763) q[0];
rz(-1.1866456) q[2];
sx q[2];
rz(-1.5015366) q[2];
sx q[2];
rz(-2.4406976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1783274) q[1];
sx q[1];
rz(-1.3849003) q[1];
sx q[1];
rz(-0.16190298) q[1];
x q[2];
rz(-1.8454148) q[3];
sx q[3];
rz(-1.7227968) q[3];
sx q[3];
rz(-0.27384588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9049412) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(-2.2616852) q[2];
rz(0.044163477) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(-0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1039466) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(-0.049731072) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(1.0838881) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7679374) q[0];
sx q[0];
rz(-2.8345788) q[0];
sx q[0];
rz(0.92371773) q[0];
x q[1];
rz(-1.8439699) q[2];
sx q[2];
rz(-1.3285085) q[2];
sx q[2];
rz(0.66852409) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6127527) q[1];
sx q[1];
rz(-1.4626979) q[1];
sx q[1];
rz(-2.0987064) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51387042) q[3];
sx q[3];
rz(-1.5268945) q[3];
sx q[3];
rz(-1.9572789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9087387) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(-0.24442913) q[2];
rz(2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844834) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(0.094141468) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(-2.24618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6796414) q[0];
sx q[0];
rz(-1.9129176) q[0];
sx q[0];
rz(-1.530594) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32161153) q[2];
sx q[2];
rz(-0.65882896) q[2];
sx q[2];
rz(2.4668601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6678644) q[1];
sx q[1];
rz(-1.5267936) q[1];
sx q[1];
rz(1.4021137) q[1];
x q[2];
rz(1.6454837) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(0.80319476) q[2];
rz(1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(-2.7289594) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068709277) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(-2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(1.2566459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4706659) q[0];
sx q[0];
rz(-1.6111272) q[0];
sx q[0];
rz(-1.3192024) q[0];
x q[1];
rz(-0.2101375) q[2];
sx q[2];
rz(-1.1147095) q[2];
sx q[2];
rz(1.2880585) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41374846) q[1];
sx q[1];
rz(-1.2894948) q[1];
sx q[1];
rz(-1.6378228) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4323746) q[3];
sx q[3];
rz(-1.8172482) q[3];
sx q[3];
rz(2.646288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98823035) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(1.3640277) q[2];
rz(-0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751188) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(-2.8334154) q[0];
rz(3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(-0.38696188) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.455084) q[0];
sx q[0];
rz(-2.5065656) q[0];
sx q[0];
rz(1.9576661) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53051051) q[2];
sx q[2];
rz(-0.94420099) q[2];
sx q[2];
rz(-0.26041398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9427529) q[1];
sx q[1];
rz(-1.4958592) q[1];
sx q[1];
rz(2.4692175) q[1];
x q[2];
rz(-2.6450063) q[3];
sx q[3];
rz(-1.7644617) q[3];
sx q[3];
rz(0.95369875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62362921) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(-2.0429042) q[0];
rz(2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44137529) q[0];
sx q[0];
rz(-2.2451631) q[0];
sx q[0];
rz(-2.3250439) q[0];
rz(0.49528893) q[2];
sx q[2];
rz(-2.1914748) q[2];
sx q[2];
rz(-0.6492614) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.038736131) q[1];
sx q[1];
rz(-0.68118459) q[1];
sx q[1];
rz(-0.78050651) q[1];
x q[2];
rz(-2.4008972) q[3];
sx q[3];
rz(-1.9441838) q[3];
sx q[3];
rz(1.545056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.07842841) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(-0.72193974) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(-0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14770517) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(-2.06185) q[0];
rz(2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.7396897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3106874) q[0];
sx q[0];
rz(-1.7831823) q[0];
sx q[0];
rz(3.0124245) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14342587) q[2];
sx q[2];
rz(-2.9529245) q[2];
sx q[2];
rz(0.27486899) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55587308) q[1];
sx q[1];
rz(-1.8750637) q[1];
sx q[1];
rz(0.85822206) q[1];
x q[2];
rz(-0.63125061) q[3];
sx q[3];
rz(-1.7389986) q[3];
sx q[3];
rz(-1.9025308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(-1.520291) q[2];
rz(2.5907717) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-2.3475636) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(3.1237596) q[3];
sx q[3];
rz(-2.087839) q[3];
sx q[3];
rz(2.9057333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
