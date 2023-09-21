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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5380733) q[0];
sx q[0];
rz(-0.78342122) q[0];
sx q[0];
rz(2.6253683) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3697882) q[2];
sx q[2];
rz(-0.82320854) q[2];
sx q[2];
rz(-2.8231205) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77820233) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(-1.1671288) q[1];
x q[2];
rz(1.6463046) q[3];
sx q[3];
rz(-1.8823874) q[3];
sx q[3];
rz(0.15234767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26596507) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(1.3226091) q[2];
rz(2.8406075) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(-1.0377201) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25181928) q[0];
sx q[0];
rz(-2.5368241) q[0];
sx q[0];
rz(0.39321995) q[0];
x q[1];
rz(-2.6056387) q[2];
sx q[2];
rz(-1.1079271) q[2];
sx q[2];
rz(-3.0218389) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2270826) q[1];
sx q[1];
rz(-2.7343379) q[1];
sx q[1];
rz(-1.6362908) q[1];
rz(0.95275767) q[3];
sx q[3];
rz(-0.85988322) q[3];
sx q[3];
rz(3.0666805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(3.0569055) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(0.95570046) q[0];
rz(0.39069191) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(0.57317615) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7645435) q[0];
sx q[0];
rz(-2.0007613) q[0];
sx q[0];
rz(-3.0717875) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8163221) q[2];
sx q[2];
rz(-0.79954445) q[2];
sx q[2];
rz(0.48649597) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0108311) q[1];
sx q[1];
rz(-1.8567994) q[1];
sx q[1];
rz(2.7421013) q[1];
x q[2];
rz(-2.0315941) q[3];
sx q[3];
rz(-1.9402383) q[3];
sx q[3];
rz(0.24584578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
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
rz(-0.36270025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1283778) q[0];
sx q[0];
rz(-0.88071874) q[0];
sx q[0];
rz(-0.31353686) q[0];
x q[1];
rz(-0.018410725) q[2];
sx q[2];
rz(-1.6022748) q[2];
sx q[2];
rz(1.5359495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1945222) q[1];
sx q[1];
rz(-2.2010989) q[1];
sx q[1];
rz(-1.6987726) q[1];
rz(2.1953771) q[3];
sx q[3];
rz(-1.5120301) q[3];
sx q[3];
rz(1.1138686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(-2.4827042) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.1936369) q[1];
sx q[1];
rz(-1.682122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97809726) q[0];
sx q[0];
rz(-2.8601544) q[0];
sx q[0];
rz(-2.1241758) q[0];
rz(-0.95894496) q[2];
sx q[2];
rz(-2.358837) q[2];
sx q[2];
rz(-2.3252955) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0046878) q[1];
sx q[1];
rz(-1.7675752) q[1];
sx q[1];
rz(2.911527) q[1];
rz(-0.32096433) q[3];
sx q[3];
rz(-1.1043613) q[3];
sx q[3];
rz(-2.5209559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(2.2996976) q[2];
rz(-2.1250336) q[3];
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
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-0.13866436) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4045227) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(1.5674595) q[0];
x q[1];
rz(2.9406592) q[2];
sx q[2];
rz(-1.4954508) q[2];
sx q[2];
rz(-0.88694015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8530635) q[1];
sx q[1];
rz(-2.1398586) q[1];
sx q[1];
rz(2.2437375) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1451374) q[3];
sx q[3];
rz(-2.2666551) q[3];
sx q[3];
rz(-1.3130207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.8072051) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(-0.1383055) q[1];
sx q[1];
rz(-0.62430635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1034575) q[0];
sx q[0];
rz(-0.75224829) q[0];
sx q[0];
rz(0.30970807) q[0];
rz(0.46632669) q[2];
sx q[2];
rz(-1.000324) q[2];
sx q[2];
rz(-1.8898659) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.034060409) q[1];
sx q[1];
rz(-0.95648396) q[1];
sx q[1];
rz(0.88453102) q[1];
rz(-pi) q[2];
rz(-2.6631213) q[3];
sx q[3];
rz(-2.1395184) q[3];
sx q[3];
rz(2.3346321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6162993) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-2.7590511) q[2];
rz(0.031490695) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-3.0338045) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(3.0016622) q[0];
rz(1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(3.0126742) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55420586) q[0];
sx q[0];
rz(-1.6763858) q[0];
sx q[0];
rz(1.9681853) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56477408) q[2];
sx q[2];
rz(-2.161705) q[2];
sx q[2];
rz(1.5608982) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0205295) q[1];
sx q[1];
rz(-1.1464835) q[1];
sx q[1];
rz(-1.2709649) q[1];
rz(-pi) q[2];
rz(1.5548607) q[3];
sx q[3];
rz(-2.3105379) q[3];
sx q[3];
rz(0.15077886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(-1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.36214608) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(-1.8918442) q[0];
rz(-0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58383656) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(-1.7519959) q[0];
rz(-1.1662912) q[2];
sx q[2];
rz(-2.1098237) q[2];
sx q[2];
rz(1.5664958) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2495888) q[1];
sx q[1];
rz(-2.2715886) q[1];
sx q[1];
rz(-0.65111098) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25165598) q[3];
sx q[3];
rz(-2.3725384) q[3];
sx q[3];
rz(2.3264399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(-3.0272711) q[2];
rz(-1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.982622) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262912) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(2.9283438) q[0];
rz(2.7217216) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-0.54668033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0471668) q[0];
sx q[0];
rz(-1.7921653) q[0];
sx q[0];
rz(-0.84035994) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8516413) q[2];
sx q[2];
rz(-2.6486514) q[2];
sx q[2];
rz(-1.314756) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.039779546) q[1];
sx q[1];
rz(-1.2817849) q[1];
sx q[1];
rz(-1.5323557) q[1];
rz(1.7095079) q[3];
sx q[3];
rz(-1.3472392) q[3];
sx q[3];
rz(-0.15299882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(-3.1372916) q[2];
rz(2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-0.014856902) q[2];
sx q[2];
rz(-1.6179634) q[2];
sx q[2];
rz(0.85214324) q[2];
rz(-0.98106445) q[3];
sx q[3];
rz(-2.5909501) q[3];
sx q[3];
rz(0.29028374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
