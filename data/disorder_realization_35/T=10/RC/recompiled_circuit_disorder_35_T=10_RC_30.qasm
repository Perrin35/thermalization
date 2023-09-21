OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(1.024363) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5380733) q[0];
sx q[0];
rz(-2.3581714) q[0];
sx q[0];
rz(0.51622434) q[0];
rz(-pi) q[1];
rz(2.3697882) q[2];
sx q[2];
rz(-0.82320854) q[2];
sx q[2];
rz(2.8231205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4989657) q[1];
sx q[1];
rz(-1.9470012) q[1];
sx q[1];
rz(-2.7514003) q[1];
rz(0.3124247) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(-1.3952599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(1.8189836) q[2];
rz(2.8406075) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1865632) q[1];
sx q[1];
rz(1.0377201) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99012016) q[0];
sx q[0];
rz(-1.7904141) q[0];
sx q[0];
rz(-0.56818509) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53595397) q[2];
sx q[2];
rz(-2.0336656) q[2];
sx q[2];
rz(-3.0218389) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.91451007) q[1];
sx q[1];
rz(-0.40725476) q[1];
sx q[1];
rz(1.5053019) q[1];
rz(-pi) q[2];
x q[2];
rz(2.328697) q[3];
sx q[3];
rz(-1.1162236) q[3];
sx q[3];
rz(-1.9302492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(3.0569055) q[2];
rz(2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5304853) q[0];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5983551) q[0];
sx q[0];
rz(-2.7063473) q[0];
sx q[0];
rz(1.4198562) q[0];
rz(0.32527058) q[2];
sx q[2];
rz(-0.79954445) q[2];
sx q[2];
rz(2.6550967) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14904505) q[1];
sx q[1];
rz(-2.6547975) q[1];
sx q[1];
rz(2.4942314) q[1];
rz(-0.40804789) q[3];
sx q[3];
rz(-1.9983665) q[3];
sx q[3];
rz(-1.6392631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0009784) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(3.0443232) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(-2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28213421) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(-0.55066806) q[0];
rz(1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-2.7788924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013214839) q[0];
sx q[0];
rz(-2.2608739) q[0];
sx q[0];
rz(-0.31353686) q[0];
rz(-pi) q[1];
rz(-0.018410725) q[2];
sx q[2];
rz(-1.6022748) q[2];
sx q[2];
rz(1.5359495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73211654) q[1];
sx q[1];
rz(-0.64142694) q[1];
sx q[1];
rz(-0.17318053) q[1];
rz(-2.1953771) q[3];
sx q[3];
rz(-1.6295625) q[3];
sx q[3];
rz(-2.0277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.68391934) q[2];
sx q[2];
rz(-1.0689015) q[2];
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
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18810774) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(-1.4594706) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0844903) q[0];
sx q[0];
rz(-1.717289) q[0];
sx q[0];
rz(-1.811972) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2541788) q[2];
sx q[2];
rz(-1.9878584) q[2];
sx q[2];
rz(2.8487157) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1294714) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(-2.4232037) q[1];
rz(-2.0586117) q[3];
sx q[3];
rz(-1.8564463) q[3];
sx q[3];
rz(-0.80174996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.5766597) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(1.9110514) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(-0.13866436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16620557) q[0];
sx q[0];
rz(-1.5741325) q[0];
sx q[0];
rz(-0.020394527) q[0];
rz(-1.4939098) q[2];
sx q[2];
rz(-1.3704408) q[2];
sx q[2];
rz(0.69918699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0181959) q[1];
sx q[1];
rz(-1.0180078) q[1];
sx q[1];
rz(-0.68560302) q[1];
rz(-pi) q[2];
rz(-2.6824529) q[3];
sx q[3];
rz(-2.3448179) q[3];
sx q[3];
rz(1.2129984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.3343875) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3793959) q[0];
sx q[0];
rz(-1.3610098) q[0];
sx q[0];
rz(-0.72797738) q[0];
rz(0.95958556) q[2];
sx q[2];
rz(-0.72003905) q[2];
sx q[2];
rz(-2.6401273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1075322) q[1];
sx q[1];
rz(-2.1851087) q[1];
sx q[1];
rz(0.88453102) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1950486) q[3];
sx q[3];
rz(-0.72580273) q[3];
sx q[3];
rz(-3.1012227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(3.110102) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-3.0016622) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(0.12891842) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77057225) q[0];
sx q[0];
rz(-0.41045529) q[0];
sx q[0];
rz(-1.8380941) q[0];
rz(-2.2419937) q[2];
sx q[2];
rz(-1.1102144) q[2];
sx q[2];
rz(-0.34924289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0205295) q[1];
sx q[1];
rz(-1.1464835) q[1];
sx q[1];
rz(-1.2709649) q[1];
rz(-pi) q[2];
rz(2.4017879) q[3];
sx q[3];
rz(-1.5825669) q[3];
sx q[3];
rz(-1.7108325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(0.62409419) q[2];
rz(-2.9028153) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6791145) q[0];
sx q[0];
rz(-2.9318641) q[0];
sx q[0];
rz(1.0366584) q[0];
rz(-pi) q[1];
rz(0.58198858) q[2];
sx q[2];
rz(-2.4798923) q[2];
sx q[2];
rz(0.87994196) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13547922) q[1];
sx q[1];
rz(-2.0524426) q[1];
sx q[1];
rz(2.3856132) q[1];
rz(-0.25165598) q[3];
sx q[3];
rz(-0.76905426) q[3];
sx q[3];
rz(-0.8151527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0525557) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(0.11432153) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.982622) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91530144) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(2.9283438) q[0];
rz(-2.7217216) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-2.5949123) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0471668) q[0];
sx q[0];
rz(-1.7921653) q[0];
sx q[0];
rz(0.84035994) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7231862) q[2];
sx q[2];
rz(-1.1001462) q[2];
sx q[2];
rz(1.5002804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.17391275) q[1];
sx q[1];
rz(-2.850107) q[1];
sx q[1];
rz(3.013054) q[1];
x q[2];
rz(1.4320847) q[3];
sx q[3];
rz(-1.7943534) q[3];
sx q[3];
rz(2.9885938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13835779) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(-0.004301087) q[2];
rz(2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99988408) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(-0.44395631) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(-0.014856902) q[2];
sx q[2];
rz(-1.6179634) q[2];
sx q[2];
rz(0.85214324) q[2];
rz(-2.1605282) q[3];
sx q[3];
rz(-0.55064252) q[3];
sx q[3];
rz(-2.8513089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
