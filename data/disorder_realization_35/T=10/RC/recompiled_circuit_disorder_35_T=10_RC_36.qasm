OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(4.5067956) q[0];
sx q[0];
rz(11.542008) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(1.1693118) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071872358) q[0];
sx q[0];
rz(-0.9099996) q[0];
sx q[0];
rz(1.1138492) q[0];
rz(-pi) q[1];
rz(0.92584445) q[2];
sx q[2];
rz(-1.0174417) q[2];
sx q[2];
rz(2.4990621) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.484326) q[1];
sx q[1];
rz(-2.6063759) q[1];
sx q[1];
rz(-0.80429299) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4952881) q[3];
sx q[3];
rz(-1.8823874) q[3];
sx q[3];
rz(0.15234767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26596507) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(-0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009636119) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(-1.0377201) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71887165) q[0];
sx q[0];
rz(-1.0178716) q[0];
sx q[0];
rz(1.8296815) q[0];
rz(-pi) q[1];
rz(-2.3677164) q[2];
sx q[2];
rz(-2.4485588) q[2];
sx q[2];
rz(-2.3351923) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4251551) q[1];
sx q[1];
rz(-1.5967224) q[1];
sx q[1];
rz(-1.9772711) q[1];
rz(-pi) q[2];
rz(0.81289566) q[3];
sx q[3];
rz(-1.1162236) q[3];
sx q[3];
rz(1.9302492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4804046) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(-1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5304853) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(-2.5684165) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1646106) q[0];
sx q[0];
rz(-1.6342388) q[0];
sx q[0];
rz(2.001686) q[0];
rz(2.3689752) q[2];
sx q[2];
rz(-1.801991) q[2];
sx q[2];
rz(0.85341838) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1307615) q[1];
sx q[1];
rz(-1.2847932) q[1];
sx q[1];
rz(-2.7421013) q[1];
rz(-0.8543386) q[3];
sx q[3];
rz(-2.5594098) q[3];
sx q[3];
rz(-2.4454988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(-2.9476681) q[2];
rz(3.0443232) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28213421) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(2.5909246) q[0];
rz(1.1286873) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(2.7788924) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.761128) q[0];
sx q[0];
rz(-1.8109545) q[0];
sx q[0];
rz(-2.2855177) q[0];
x q[1];
rz(-1.0417468) q[2];
sx q[2];
rz(-0.036465557) q[2];
sx q[2];
rz(2.1349825) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4421694) q[1];
sx q[1];
rz(-1.6740834) q[1];
sx q[1];
rz(2.5073754) q[1];
rz(-pi) q[2];
rz(0.94621559) q[3];
sx q[3];
rz(-1.5120301) q[3];
sx q[3];
rz(-1.1138686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-0.80250424) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(3.127393) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.682122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97809726) q[0];
sx q[0];
rz(-0.28143829) q[0];
sx q[0];
rz(-1.0174169) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2541788) q[2];
sx q[2];
rz(-1.1537342) q[2];
sx q[2];
rz(2.8487157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1294714) q[1];
sx q[1];
rz(-2.8399889) q[1];
sx q[1];
rz(-2.4232037) q[1];
x q[2];
rz(-1.0829809) q[3];
sx q[3];
rz(-1.2851464) q[3];
sx q[3];
rz(2.3398427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(2.2996976) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275948) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(-2.5573964) q[0];
rz(-1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-3.0029283) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5748782) q[0];
sx q[0];
rz(-3.1209271) q[0];
sx q[0];
rz(-2.979435) q[0];
rz(-pi) q[1];
rz(1.4939098) q[2];
sx q[2];
rz(-1.7711519) q[2];
sx q[2];
rz(0.69918699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2648592) q[1];
sx q[1];
rz(-0.85163341) q[1];
sx q[1];
rz(-2.369146) q[1];
rz(-pi) q[2];
rz(2.3994282) q[3];
sx q[3];
rz(-1.8932749) q[3];
sx q[3];
rz(0.025067586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(1.3343875) q[2];
rz(-1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(2.6050674) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(0.62430635) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1034575) q[0];
sx q[0];
rz(-2.3893444) q[0];
sx q[0];
rz(0.30970807) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1937218) q[2];
sx q[2];
rz(-1.1827173) q[2];
sx q[2];
rz(0.58448234) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0956456) q[1];
sx q[1];
rz(-2.1150757) q[1];
sx q[1];
rz(0.73927684) q[1];
rz(-pi) q[2];
rz(2.1948366) q[3];
sx q[3];
rz(-1.1723926) q[3];
sx q[3];
rz(-2.1053673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6162993) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(2.7590511) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(3.0338045) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.265825) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(0.13993046) q[0];
rz(1.6775999) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(-3.0126742) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873868) q[0];
sx q[0];
rz(-1.4652068) q[0];
sx q[0];
rz(1.1734074) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56477408) q[2];
sx q[2];
rz(-0.97988765) q[2];
sx q[2];
rz(1.5608982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.376437) q[1];
sx q[1];
rz(-0.51424485) q[1];
sx q[1];
rz(0.57904412) q[1];
rz(0.017459004) q[3];
sx q[3];
rz(-2.401712) q[3];
sx q[3];
rz(0.12714126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.8918442) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58383656) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(1.7519959) q[0];
x q[1];
rz(-1.1662912) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(-1.5664958) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0235325) q[1];
sx q[1];
rz(-0.91731056) q[1];
sx q[1];
rz(-2.1937624) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3343072) q[3];
sx q[3];
rz(-0.83179501) q[3];
sx q[3];
rz(0.47154271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(-0.11432153) q[2];
rz(-1.8814686) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(1.1589706) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(2.7217216) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(2.5949123) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4237758) q[0];
sx q[0];
rz(-0.86200889) q[0];
sx q[0];
rz(2.8481759) q[0];
rz(-1.7231862) q[2];
sx q[2];
rz(-1.1001462) q[2];
sx q[2];
rz(1.5002804) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5200561) q[1];
sx q[1];
rz(-1.6076419) q[1];
sx q[1];
rz(2.8523793) q[1];
rz(2.5952026) q[3];
sx q[3];
rz(-0.26248172) q[3];
sx q[3];
rz(2.4266092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(-3.1372916) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(2.6314541) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99988408) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-1.6179686) q[2];
sx q[2];
rz(-1.5559559) q[2];
sx q[2];
rz(2.422239) q[2];
rz(2.1605282) q[3];
sx q[3];
rz(-2.5909501) q[3];
sx q[3];
rz(0.29028374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
