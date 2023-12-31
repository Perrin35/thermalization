OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(0.00014076509) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(-2.1773832) q[1];
sx q[1];
rz(1.1934086) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084273987) q[0];
sx q[0];
rz(-2.8017375) q[0];
sx q[0];
rz(1.0124595) q[0];
x q[1];
rz(2.675406) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(-2.8592062) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4729893) q[1];
sx q[1];
rz(-2.0358634) q[1];
sx q[1];
rz(-2.4238062) q[1];
rz(-pi) q[2];
rz(-3.0317806) q[3];
sx q[3];
rz(-1.3545274) q[3];
sx q[3];
rz(-3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.45941916) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.9127282) q[2];
rz(1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(-1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(-1.0128101) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(2.0181296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7682122) q[0];
sx q[0];
rz(-2.0313615) q[0];
sx q[0];
rz(3.0773613) q[0];
rz(-2.3617619) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(-2.0335399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9065735) q[1];
sx q[1];
rz(-1.0636914) q[1];
sx q[1];
rz(2.1706235) q[1];
rz(-pi) q[2];
rz(1.0081604) q[3];
sx q[3];
rz(-2.1481272) q[3];
sx q[3];
rz(-1.2853704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(-2.222555) q[2];
rz(0.67409003) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(-1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(0.69349849) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(-1.1330053) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5223761) q[0];
sx q[0];
rz(-1.7130865) q[0];
sx q[0];
rz(-0.0033358047) q[0];
rz(-pi) q[1];
rz(2.3511231) q[2];
sx q[2];
rz(-1.1970453) q[2];
sx q[2];
rz(2.8369396) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6269835) q[1];
sx q[1];
rz(-1.2830462) q[1];
sx q[1];
rz(-1.0088374) q[1];
rz(-0.95136178) q[3];
sx q[3];
rz(-2.1054483) q[3];
sx q[3];
rz(1.9922647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8901849) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816198) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(-1.9807293) q[0];
rz(2.2456031) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(0.13555759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40977749) q[0];
sx q[0];
rz(-1.2588132) q[0];
sx q[0];
rz(-0.72809763) q[0];
rz(-pi) q[1];
rz(0.074684871) q[2];
sx q[2];
rz(-1.9539781) q[2];
sx q[2];
rz(2.243724) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96326522) q[1];
sx q[1];
rz(-1.3849003) q[1];
sx q[1];
rz(-2.9796897) q[1];
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
rz(0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039466) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(1.0132382) q[0];
rz(0.049731072) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-1.0838881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5686544) q[0];
sx q[0];
rz(-1.3875811) q[0];
sx q[0];
rz(-1.818548) q[0];
rz(1.8439699) q[2];
sx q[2];
rz(-1.8130842) q[2];
sx q[2];
rz(-2.4730686) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6127527) q[1];
sx q[1];
rz(-1.4626979) q[1];
sx q[1];
rz(-2.0987064) q[1];
x q[2];
rz(-3.0524592) q[3];
sx q[3];
rz(-0.51557487) q[3];
sx q[3];
rz(-2.6775132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.23285398) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(2.8971635) q[2];
rz(0.43236732) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2844834) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(0.094141468) q[0];
rz(0.17177467) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(0.89541268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6796414) q[0];
sx q[0];
rz(-1.2286751) q[0];
sx q[0];
rz(1.530594) q[0];
rz(-pi) q[1];
rz(-0.32161153) q[2];
sx q[2];
rz(-0.65882896) q[2];
sx q[2];
rz(-2.4668601) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10455924) q[1];
sx q[1];
rz(-1.4022786) q[1];
sx q[1];
rz(0.044635459) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6454837) q[3];
sx q[3];
rz(-1.5731305) q[3];
sx q[3];
rz(-2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0078997) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(-1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(-1.8849467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.885868) q[0];
sx q[0];
rz(-0.2547383) q[0];
sx q[0];
rz(1.4101009) q[0];
rz(-1.1058544) q[2];
sx q[2];
rz(-1.3824116) q[2];
sx q[2];
rz(2.9525194) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9650967) q[1];
sx q[1];
rz(-2.8526222) q[1];
sx q[1];
rz(-2.9138336) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36862936) q[3];
sx q[3];
rz(-2.3978524) q[3];
sx q[3];
rz(-1.7891235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98823035) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751188) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(2.7546308) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.455084) q[0];
sx q[0];
rz(-2.5065656) q[0];
sx q[0];
rz(-1.9576661) q[0];
rz(0.96078028) q[2];
sx q[2];
rz(-0.79723251) q[2];
sx q[2];
rz(2.0955992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-2.4657001) q[1];
sx q[1];
rz(3.0216316) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49658637) q[3];
sx q[3];
rz(-1.3771309) q[3];
sx q[3];
rz(-2.1878939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62362921) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(-0.44000885) q[2];
rz(2.4258339) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(-2.4138342) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(0.049302014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5987199) q[0];
sx q[0];
rz(-2.1763986) q[0];
sx q[0];
rz(-2.4332895) q[0];
x q[1];
rz(-2.2531613) q[2];
sx q[2];
rz(-1.9677791) q[2];
sx q[2];
rz(-1.9156485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1028565) q[1];
sx q[1];
rz(-0.68118459) q[1];
sx q[1];
rz(0.78050651) q[1];
rz(-2.0588166) q[3];
sx q[3];
rz(-0.89142311) q[3];
sx q[3];
rz(2.8454012) q[3];
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
rz(0.94349629) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9938875) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(-1.0797427) q[0];
rz(2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(1.4019029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.429075) q[0];
sx q[0];
rz(-1.6970465) q[0];
sx q[0];
rz(-1.7849126) q[0];
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
rz(-2.5857196) q[1];
sx q[1];
rz(-1.266529) q[1];
sx q[1];
rz(2.2833706) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8614282) q[3];
sx q[3];
rz(-2.4912842) q[3];
sx q[3];
rz(0.10661099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(-1.6213017) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(-0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.993492) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(-0.91611721) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(0.79402906) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(-0.017833088) q[3];
sx q[3];
rz(-2.087839) q[3];
sx q[3];
rz(2.9057333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
