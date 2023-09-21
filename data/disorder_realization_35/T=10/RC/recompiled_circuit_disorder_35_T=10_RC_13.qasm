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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7919851) q[0];
sx q[0];
rz(-1.2149997) q[0];
sx q[0];
rz(0.7138568) q[0];
rz(2.4835303) q[2];
sx q[2];
rz(-1.0339289) q[2];
sx q[2];
rz(1.8368349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4989657) q[1];
sx q[1];
rz(-1.9470012) q[1];
sx q[1];
rz(0.39019231) q[1];
x q[2];
rz(-0.3124247) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(-1.7463328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(1.3226091) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.1319565) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(-2.6665376) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(-2.1038726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8897734) q[0];
sx q[0];
rz(-2.5368241) q[0];
sx q[0];
rz(-0.39321995) q[0];
rz(-pi) q[1];
rz(0.53595397) q[2];
sx q[2];
rz(-1.1079271) q[2];
sx q[2];
rz(0.11975372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91451007) q[1];
sx q[1];
rz(-0.40725476) q[1];
sx q[1];
rz(1.6362908) q[1];
rz(-pi) q[2];
rz(-2.188835) q[3];
sx q[3];
rz(-2.2817094) q[3];
sx q[3];
rz(0.074912138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4804046) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(-0.37880138) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(-2.5684165) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5983551) q[0];
sx q[0];
rz(-2.7063473) q[0];
sx q[0];
rz(-1.7217365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2531883) q[2];
sx q[2];
rz(-2.3177958) q[2];
sx q[2];
rz(2.2044646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1307615) q[1];
sx q[1];
rz(-1.2847932) q[1];
sx q[1];
rz(2.7421013) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8543386) q[3];
sx q[3];
rz(-0.58218282) q[3];
sx q[3];
rz(-2.4454988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8594584) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(-2.5909246) q[0];
rz(-2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(0.36270025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6837316) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(-1.213221) q[0];
x q[1];
rz(0.018410725) q[2];
sx q[2];
rz(-1.6022748) q[2];
sx q[2];
rz(1.6056431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4094761) q[1];
sx q[1];
rz(-2.5001657) q[1];
sx q[1];
rz(0.17318053) q[1];
rz(-pi) q[2];
rz(-0.072399541) q[3];
sx q[3];
rz(-0.94745938) q[3];
sx q[3];
rz(2.6423531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4576733) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18810774) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(-0.017379934) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(1.4594706) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.592011) q[0];
sx q[0];
rz(-1.8093384) q[0];
sx q[0];
rz(2.9907988) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2541788) q[2];
sx q[2];
rz(-1.1537342) q[2];
sx q[2];
rz(-2.8487157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1294714) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(0.71838897) q[1];
rz(-1.0829809) q[3];
sx q[3];
rz(-1.2851464) q[3];
sx q[3];
rz(-0.80174996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.83546272) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(3.1275948) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(-1.9110514) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(3.0029283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9753871) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(-3.1211981) q[0];
rz(-pi) q[1];
rz(2.7799941) q[2];
sx q[2];
rz(-2.9271759) q[2];
sx q[2];
rz(2.8117361) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2648592) q[1];
sx q[1];
rz(-2.2899592) q[1];
sx q[1];
rz(-0.77244669) q[1];
rz(-pi) q[2];
rz(-2.3994282) q[3];
sx q[3];
rz(-1.2483178) q[3];
sx q[3];
rz(0.025067586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77506322) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.3343875) q[2];
rz(1.1602317) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(3.0464723) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5360864) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(-0.62430635) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3793959) q[0];
sx q[0];
rz(-1.3610098) q[0];
sx q[0];
rz(-2.4136153) q[0];
rz(-0.94787089) q[2];
sx q[2];
rz(-1.1827173) q[2];
sx q[2];
rz(0.58448234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.149951) q[1];
sx q[1];
rz(-0.88643622) q[1];
sx q[1];
rz(-0.73189862) q[1];
x q[2];
rz(-0.47847139) q[3];
sx q[3];
rz(-1.0020743) q[3];
sx q[3];
rz(2.3346321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6162993) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-2.7590511) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0607973) q[0];
sx q[0];
rz(-1.1757438) q[0];
sx q[0];
rz(3.0271544) q[0];
rz(-pi) q[1];
rz(2.5768186) q[2];
sx q[2];
rz(-2.161705) q[2];
sx q[2];
rz(1.5608982) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7651556) q[1];
sx q[1];
rz(-0.51424485) q[1];
sx q[1];
rz(0.57904412) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4017879) q[3];
sx q[3];
rz(-1.5825669) q[3];
sx q[3];
rz(-1.4307601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(-0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58383656) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(1.3895967) q[0];
rz(-2.5647854) q[2];
sx q[2];
rz(-1.9153321) q[2];
sx q[2];
rz(2.92958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2495888) q[1];
sx q[1];
rz(-0.87000404) q[1];
sx q[1];
rz(2.4904817) q[1];
rz(-pi) q[2];
rz(1.8072855) q[3];
sx q[3];
rz(-2.3097976) q[3];
sx q[3];
rz(0.47154271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(-1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(0.21324883) q[0];
rz(-2.7217216) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(0.54668033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(0.84035994) q[0];
rz(-pi) q[1];
rz(1.7231862) q[2];
sx q[2];
rz(-2.0414464) q[2];
sx q[2];
rz(-1.6413123) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.039779546) q[1];
sx q[1];
rz(-1.2817849) q[1];
sx q[1];
rz(1.6092369) q[1];
rz(-0.54639001) q[3];
sx q[3];
rz(-2.8791109) q[3];
sx q[3];
rz(0.7149834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1417086) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(2.6976363) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(-1.8757204) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
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