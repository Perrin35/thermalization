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
rz(-1.7763897) q[0];
sx q[0];
rz(2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(2.6095698) q[1];
sx q[1];
rz(11.397059) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3496075) q[0];
sx q[0];
rz(-1.926593) q[0];
sx q[0];
rz(-2.4277359) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92584445) q[2];
sx q[2];
rz(-2.1241509) q[2];
sx q[2];
rz(0.64253053) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3633903) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(-1.9744639) q[1];
x q[2];
rz(1.4952881) q[3];
sx q[3];
rz(-1.2592053) q[3];
sx q[3];
rz(-2.989245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.26596507) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(-0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009636119) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(1.0377201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1514725) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(-0.56818509) q[0];
rz(-pi) q[1];
rz(2.6056387) q[2];
sx q[2];
rz(-1.1079271) q[2];
sx q[2];
rz(3.0218389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2983919) q[1];
sx q[1];
rz(-1.9771264) q[1];
sx q[1];
rz(3.1133679) q[1];
rz(0.81289566) q[3];
sx q[3];
rz(-1.1162236) q[3];
sx q[3];
rz(1.9302492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4804046) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(3.0569055) q[2];
rz(0.37880138) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(-1.9975196) q[3];
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
rz(1.5304853) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(-0.95570046) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5983551) q[0];
sx q[0];
rz(-2.7063473) q[0];
sx q[0];
rz(-1.4198562) q[0];
x q[1];
rz(2.3689752) q[2];
sx q[2];
rz(-1.801991) q[2];
sx q[2];
rz(-2.2881743) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5830071) q[1];
sx q[1];
rz(-1.188394) q[1];
sx q[1];
rz(1.8797727) q[1];
x q[2];
rz(-0.8543386) q[3];
sx q[3];
rz(-0.58218282) q[3];
sx q[3];
rz(2.4454988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(0.19392459) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28213421) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-0.36270025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013214839) q[0];
sx q[0];
rz(-2.2608739) q[0];
sx q[0];
rz(-0.31353686) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0417468) q[2];
sx q[2];
rz(-0.036465557) q[2];
sx q[2];
rz(2.1349825) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1945222) q[1];
sx q[1];
rz(-2.2010989) q[1];
sx q[1];
rz(1.6987726) q[1];
rz(-pi) q[2];
rz(1.6710715) q[3];
sx q[3];
rz(-2.5146211) q[3];
sx q[3];
rz(0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(3.1385699) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18810774) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(1.4594706) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.592011) q[0];
sx q[0];
rz(-1.3322543) q[0];
sx q[0];
rz(0.15079389) q[0];
rz(-pi) q[1];
rz(-2.1826477) q[2];
sx q[2];
rz(-0.78275567) q[2];
sx q[2];
rz(-2.3252955) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0121213) q[1];
sx q[1];
rz(-2.8399889) q[1];
sx q[1];
rz(2.4232037) q[1];
x q[2];
rz(-0.32096433) q[3];
sx q[3];
rz(-2.0372314) q[3];
sx q[3];
rz(2.5209559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(2.1250336) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(1.564933) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1275948) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(-1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-3.0029283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4045227) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(-1.5674595) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20093341) q[2];
sx q[2];
rz(-1.6461419) q[2];
sx q[2];
rz(-2.2546525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8530635) q[1];
sx q[1];
rz(-2.1398586) q[1];
sx q[1];
rz(-0.89785518) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3994282) q[3];
sx q[3];
rz(-1.2483178) q[3];
sx q[3];
rz(-0.025067586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.3343875) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(-3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(0.53652525) q[0];
rz(0.58553186) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(-2.5172863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37492232) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(-1.8486345) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.675266) q[2];
sx q[2];
rz(-2.1412686) q[2];
sx q[2];
rz(1.8898659) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0956456) q[1];
sx q[1];
rz(-2.1150757) q[1];
sx q[1];
rz(0.73927684) q[1];
rz(2.1950486) q[3];
sx q[3];
rz(-0.72580273) q[3];
sx q[3];
rz(3.1012227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5252934) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(2.7590511) q[2];
rz(0.031490695) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.265825) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-3.0016622) q[0];
rz(1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(-0.12891842) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3710204) q[0];
sx q[0];
rz(-2.7311374) q[0];
sx q[0];
rz(-1.8380941) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56477408) q[2];
sx q[2];
rz(-0.97988765) q[2];
sx q[2];
rz(-1.5608982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8184549) q[1];
sx q[1];
rz(-1.2982681) q[1];
sx q[1];
rz(-2.6998991) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5548607) q[3];
sx q[3];
rz(-0.83105479) q[3];
sx q[3];
rz(0.15077886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63697469) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(2.9028153) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.8918442) q[0];
rz(0.016013913) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(1.790766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58383656) q[0];
sx q[0];
rz(-1.6769874) q[0];
sx q[0];
rz(1.7519959) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5596041) q[2];
sx q[2];
rz(-0.66170035) q[2];
sx q[2];
rz(-0.87994196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0061134) q[1];
sx q[1];
rz(-2.0524426) q[1];
sx q[1];
rz(-2.3856132) q[1];
x q[2];
rz(-1.8072855) q[3];
sx q[3];
rz(-2.3097976) q[3];
sx q[3];
rz(-0.47154271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(-0.11432153) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.1589706) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(0.54668033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0944259) q[0];
sx q[0];
rz(-1.7921653) q[0];
sx q[0];
rz(0.84035994) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6662153) q[2];
sx q[2];
rz(-1.7065085) q[2];
sx q[2];
rz(3.140608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1018131) q[1];
sx q[1];
rz(-1.8598078) q[1];
sx q[1];
rz(-1.5323557) q[1];
rz(-pi) q[2];
rz(-0.22565266) q[3];
sx q[3];
rz(-1.7060346) q[3];
sx q[3];
rz(1.7547363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(0.004301087) q[2];
rz(0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(2.6314541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99988408) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-1.5236241) q[2];
sx q[2];
rz(-1.5856367) q[2];
sx q[2];
rz(-0.71935364) q[2];
rz(1.0989582) q[3];
sx q[3];
rz(-1.2755339) q[3];
sx q[3];
rz(-0.76225029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];