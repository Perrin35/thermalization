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
rz(2.6095698) q[1];
sx q[1];
rz(11.397059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7919851) q[0];
sx q[0];
rz(-1.926593) q[0];
sx q[0];
rz(-2.4277359) q[0];
rz(-pi) q[1];
rz(0.77180441) q[2];
sx q[2];
rz(-0.82320854) q[2];
sx q[2];
rz(-2.8231205) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77820233) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(-1.9744639) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4952881) q[3];
sx q[3];
rz(-1.8823874) q[3];
sx q[3];
rz(2.989245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8756276) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(1.8189836) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-0.009636119) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(-2.1038726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.422721) q[0];
sx q[0];
rz(-1.0178716) q[0];
sx q[0];
rz(-1.3119112) q[0];
rz(-pi) q[1];
rz(-0.77387626) q[2];
sx q[2];
rz(-0.69303382) q[2];
sx q[2];
rz(-2.3351923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4251551) q[1];
sx q[1];
rz(-1.5448703) q[1];
sx q[1];
rz(1.1643216) q[1];
rz(-pi) q[2];
rz(-2.5492937) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(-2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(-1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5304853) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(-0.95570046) q[0];
rz(0.39069191) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(-2.5684165) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5983551) q[0];
sx q[0];
rz(-2.7063473) q[0];
sx q[0];
rz(1.7217365) q[0];
rz(-pi) q[1];
rz(-0.77261749) q[2];
sx q[2];
rz(-1.801991) q[2];
sx q[2];
rz(-2.2881743) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9925476) q[1];
sx q[1];
rz(-2.6547975) q[1];
sx q[1];
rz(2.4942314) q[1];
rz(-pi) q[2];
rz(-1.1099986) q[3];
sx q[3];
rz(-1.9402383) q[3];
sx q[3];
rz(2.8957469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(0.19392459) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3804647) q[0];
sx q[0];
rz(-1.3306381) q[0];
sx q[0];
rz(-0.856075) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0998459) q[2];
sx q[2];
rz(-3.1051271) q[2];
sx q[2];
rz(2.1349825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94707045) q[1];
sx q[1];
rz(-2.2010989) q[1];
sx q[1];
rz(1.6987726) q[1];
rz(-pi) q[2];
rz(-1.6710715) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(-2.7659622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(2.4827042) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-0.014199646) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(1.682122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5495816) q[0];
sx q[0];
rz(-1.8093384) q[0];
sx q[0];
rz(-0.15079389) q[0];
rz(-2.622501) q[2];
sx q[2];
rz(-2.1862098) q[2];
sx q[2];
rz(-1.5965243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1294714) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(-2.4232037) q[1];
rz(-pi) q[2];
rz(-1.0110537) q[3];
sx q[3];
rz(-2.5821745) q[3];
sx q[3];
rz(1.884348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013997812) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(-3.0029283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16620557) q[0];
sx q[0];
rz(-1.5741325) q[0];
sx q[0];
rz(0.020394527) q[0];
rz(-pi) q[1];
rz(-1.4939098) q[2];
sx q[2];
rz(-1.7711519) q[2];
sx q[2];
rz(-0.69918699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2648592) q[1];
sx q[1];
rz(-2.2899592) q[1];
sx q[1];
rz(-2.369146) q[1];
rz(-0.45913978) q[3];
sx q[3];
rz(-0.79677478) q[3];
sx q[3];
rz(1.2129984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.5360864) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(-pi/2) q[2];
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
x q[1];
rz(2.675266) q[2];
sx q[2];
rz(-2.1412686) q[2];
sx q[2];
rz(1.2517267) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1075322) q[1];
sx q[1];
rz(-0.95648396) q[1];
sx q[1];
rz(0.88453102) q[1];
x q[2];
rz(-0.94675605) q[3];
sx q[3];
rz(-1.1723926) q[3];
sx q[3];
rz(1.0362253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5252934) q[2];
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
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87576762) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(0.13993046) q[0];
rz(1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(3.0126742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0807954) q[0];
sx q[0];
rz(-1.9658488) q[0];
sx q[0];
rz(-3.0271544) q[0];
rz(0.56477408) q[2];
sx q[2];
rz(-2.161705) q[2];
sx q[2];
rz(-1.5608982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.376437) q[1];
sx q[1];
rz(-2.6273478) q[1];
sx q[1];
rz(-2.5625485) q[1];
rz(-2.4017879) q[3];
sx q[3];
rz(-1.5590258) q[3];
sx q[3];
rz(1.4307601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(-2.9028153) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.8918442) q[0];
rz(-0.016013913) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(-1.790766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6791145) q[0];
sx q[0];
rz(-2.9318641) q[0];
sx q[0];
rz(2.1049343) q[0];
rz(-2.5596041) q[2];
sx q[2];
rz(-2.4798923) q[2];
sx q[2];
rz(0.87994196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0061134) q[1];
sx q[1];
rz(-1.0891501) q[1];
sx q[1];
rz(0.7559795) q[1];
rz(0.25165598) q[3];
sx q[3];
rz(-0.76905426) q[3];
sx q[3];
rz(-2.3264399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(0.54668033) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28323805) q[0];
sx q[0];
rz(-0.75728098) q[0];
sx q[0];
rz(1.2454633) q[0];
x q[1];
rz(2.6662153) q[2];
sx q[2];
rz(-1.7065085) q[2];
sx q[2];
rz(-0.00098468653) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.039779546) q[1];
sx q[1];
rz(-1.2817849) q[1];
sx q[1];
rz(-1.6092369) q[1];
rz(-pi) q[2];
rz(-0.22565266) q[3];
sx q[3];
rz(-1.7060346) q[3];
sx q[3];
rz(-1.3868563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(-2.1440078) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99988408) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(-0.44395631) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(1.8757204) q[2];
sx q[2];
rz(-3.0921428) q[2];
sx q[2];
rz(-2.5947239) q[2];
rz(-2.0426345) q[3];
sx q[3];
rz(-1.2755339) q[3];
sx q[3];
rz(-0.76225029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
