OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(3.3942437) q[0];
sx q[0];
rz(10.630339) q[0];
rz(2.0134917) q[1];
sx q[1];
rz(-1.3265346) q[1];
sx q[1];
rz(0.74572745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0416314) q[0];
sx q[0];
rz(-1.3486762) q[0];
sx q[0];
rz(3.0801386) q[0];
rz(-1.5775852) q[2];
sx q[2];
rz(-1.7951843) q[2];
sx q[2];
rz(2.9780088) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89167833) q[1];
sx q[1];
rz(-0.39977705) q[1];
sx q[1];
rz(-1.0067902) q[1];
x q[2];
rz(-2.0449355) q[3];
sx q[3];
rz(-1.7563987) q[3];
sx q[3];
rz(1.7222278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2241609) q[2];
sx q[2];
rz(-1.5459583) q[2];
sx q[2];
rz(-1.6708299) q[2];
rz(1.6338232) q[3];
sx q[3];
rz(-3.1247415) q[3];
sx q[3];
rz(-0.9233709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725752) q[0];
sx q[0];
rz(-1.1976396) q[0];
sx q[0];
rz(-1.5684599) q[0];
rz(2.9720427) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(3.0068908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41551521) q[0];
sx q[0];
rz(-1.2252136) q[0];
sx q[0];
rz(-0.4536566) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7472166) q[2];
sx q[2];
rz(-0.069709502) q[2];
sx q[2];
rz(1.4331872) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1253648) q[1];
sx q[1];
rz(-1.1332336) q[1];
sx q[1];
rz(-1.9047649) q[1];
x q[2];
rz(1.1504796) q[3];
sx q[3];
rz(-2.947071) q[3];
sx q[3];
rz(-1.9357301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0160825) q[2];
sx q[2];
rz(-1.6566015) q[2];
sx q[2];
rz(-0.15277319) q[2];
rz(-1.7759391) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(0.16807817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.633054) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(-0.48164865) q[0];
rz(-0.18474361) q[1];
sx q[1];
rz(-1.3667204) q[1];
sx q[1];
rz(-0.99536037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1448875) q[0];
sx q[0];
rz(-1.2249682) q[0];
sx q[0];
rz(2.8976827) q[0];
rz(-pi) q[1];
rz(3.0929933) q[2];
sx q[2];
rz(-1.6311797) q[2];
sx q[2];
rz(2.8655665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.17653325) q[1];
sx q[1];
rz(-0.73672026) q[1];
sx q[1];
rz(0.45335575) q[1];
x q[2];
rz(-1.5257201) q[3];
sx q[3];
rz(-2.5451676) q[3];
sx q[3];
rz(3.0395122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2153726) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(-3.0769297) q[2];
rz(-2.0926545) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(1.8612727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30508405) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(-0.82103658) q[0];
rz(-0.07846421) q[1];
sx q[1];
rz(-1.4469701) q[1];
sx q[1];
rz(-2.1441114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0745463) q[0];
sx q[0];
rz(-0.95854488) q[0];
sx q[0];
rz(1.9590098) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5115963) q[2];
sx q[2];
rz(-1.534605) q[2];
sx q[2];
rz(2.9179887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6523903) q[1];
sx q[1];
rz(-1.308131) q[1];
sx q[1];
rz(-1.0242382) q[1];
rz(-pi) q[2];
rz(3.0313052) q[3];
sx q[3];
rz(-2.6607041) q[3];
sx q[3];
rz(2.4824924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0492101) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(1.9878261) q[2];
rz(2.9554101) q[3];
sx q[3];
rz(-0.081705339) q[3];
sx q[3];
rz(0.33128273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.004772923) q[0];
sx q[0];
rz(-2.3240219) q[0];
sx q[0];
rz(1.1432884) q[0];
rz(1.1413057) q[1];
sx q[1];
rz(-0.80991304) q[1];
sx q[1];
rz(-2.6236261) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.98949) q[0];
sx q[0];
rz(-0.94606646) q[0];
sx q[0];
rz(0.89390198) q[0];
rz(-pi) q[1];
rz(1.5557655) q[2];
sx q[2];
rz(-1.5698729) q[2];
sx q[2];
rz(0.24713384) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6828047) q[1];
sx q[1];
rz(-2.7066351) q[1];
sx q[1];
rz(0.87587237) q[1];
x q[2];
rz(0.2980026) q[3];
sx q[3];
rz(-1.7116705) q[3];
sx q[3];
rz(2.6948787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(-2.8138568) q[2];
rz(1.1945126) q[3];
sx q[3];
rz(-3.0142398) q[3];
sx q[3];
rz(0.85668844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1389393) q[0];
sx q[0];
rz(-0.26873538) q[0];
sx q[0];
rz(-0.56185454) q[0];
rz(-1.4909164) q[1];
sx q[1];
rz(-1.6249388) q[1];
sx q[1];
rz(0.098310016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94842978) q[0];
sx q[0];
rz(-1.5353445) q[0];
sx q[0];
rz(1.9928785) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19616429) q[2];
sx q[2];
rz(-3.1409396) q[2];
sx q[2];
rz(-3.0214423) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5206388) q[1];
sx q[1];
rz(-0.47157447) q[1];
sx q[1];
rz(-3.0388799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0344072) q[3];
sx q[3];
rz(-2.0258198) q[3];
sx q[3];
rz(-2.1912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.085122846) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(1.0657715) q[2];
rz(0.39984518) q[3];
sx q[3];
rz(-2.6078434) q[3];
sx q[3];
rz(1.2679509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9999076) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.6837233) q[0];
rz(1.1204002) q[1];
sx q[1];
rz(-0.1381865) q[1];
sx q[1];
rz(-0.33946005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0908112) q[0];
sx q[0];
rz(-3.0635186) q[0];
sx q[0];
rz(2.9705621) q[0];
rz(2.5869114) q[2];
sx q[2];
rz(-1.5825869) q[2];
sx q[2];
rz(1.5761216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5684699) q[1];
sx q[1];
rz(-0.085612729) q[1];
sx q[1];
rz(1.5517615) q[1];
x q[2];
rz(-1.6704329) q[3];
sx q[3];
rz(-1.908034) q[3];
sx q[3];
rz(-1.8189614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3642984) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(2.778229) q[2];
rz(-2.0511138) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(1.1183848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.61984396) q[0];
sx q[0];
rz(-2.813297) q[0];
sx q[0];
rz(-2.2186665) q[0];
rz(-1.6587616) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(-3.1373851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9516966) q[0];
sx q[0];
rz(-1.0605264) q[0];
sx q[0];
rz(-2.4416591) q[0];
rz(-pi) q[1];
rz(-1.9816859) q[2];
sx q[2];
rz(-1.5653603) q[2];
sx q[2];
rz(-1.5856575) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5397268) q[1];
sx q[1];
rz(-1.6352354) q[1];
sx q[1];
rz(-1.0947202) q[1];
x q[2];
rz(-1.6170623) q[3];
sx q[3];
rz(-2.1827201) q[3];
sx q[3];
rz(2.9718252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5620455) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(1.2008249) q[2];
rz(1.7480525) q[3];
sx q[3];
rz(-0.0037007185) q[3];
sx q[3];
rz(0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414108) q[0];
sx q[0];
rz(-0.44898471) q[0];
sx q[0];
rz(-1.9218943) q[0];
rz(1.8118743) q[1];
sx q[1];
rz(-2.0025573) q[1];
sx q[1];
rz(-0.16389287) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9326646) q[0];
sx q[0];
rz(-2.744356) q[0];
sx q[0];
rz(1.1139289) q[0];
rz(-pi) q[1];
rz(1.8437949) q[2];
sx q[2];
rz(-2.1381209) q[2];
sx q[2];
rz(-1.8253872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6299705) q[1];
sx q[1];
rz(-3.0093806) q[1];
sx q[1];
rz(-2.0396784) q[1];
rz(-pi) q[2];
rz(0.71299841) q[3];
sx q[3];
rz(-0.87942356) q[3];
sx q[3];
rz(-2.612243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4280052) q[2];
sx q[2];
rz(-3.0516477) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(-3.0981787) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6503705) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(2.6553335) q[0];
rz(-0.69475118) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(2.6659226) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2025288) q[0];
sx q[0];
rz(-0.049204218) q[0];
sx q[0];
rz(-2.3243107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86464244) q[2];
sx q[2];
rz(-1.5448031) q[2];
sx q[2];
rz(-0.038250462) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.025534677) q[1];
sx q[1];
rz(-0.60760159) q[1];
sx q[1];
rz(0.011354253) q[1];
x q[2];
rz(-0.52880295) q[3];
sx q[3];
rz(-2.7861054) q[3];
sx q[3];
rz(-2.5606009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4645369) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(-1.2976868) q[2];
rz(1.248598) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(0.36778522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1260592) q[0];
sx q[0];
rz(-1.561469) q[0];
sx q[0];
rz(-1.4373686) q[0];
rz(0.88232782) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(-3.0493564) q[2];
sx q[2];
rz(-1.606745) q[2];
sx q[2];
rz(1.8435115) q[2];
rz(-2.9024057) q[3];
sx q[3];
rz(-1.5755972) q[3];
sx q[3];
rz(1.5920873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
