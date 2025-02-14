OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1462829) q[0];
sx q[0];
rz(0.028385552) q[0];
sx q[0];
rz(10.382494) q[0];
rz(-2.7081642) q[1];
sx q[1];
rz(-1.6038916) q[1];
sx q[1];
rz(-2.5133207) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4074898) q[0];
sx q[0];
rz(-0.39034319) q[0];
sx q[0];
rz(-0.32456188) q[0];
rz(-0.45189721) q[2];
sx q[2];
rz(-1.1152281) q[2];
sx q[2];
rz(1.4939552) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1294259) q[1];
sx q[1];
rz(-0.91246997) q[1];
sx q[1];
rz(-2.7321187) q[1];
rz(1.3468993) q[3];
sx q[3];
rz(-2.7198987) q[3];
sx q[3];
rz(-1.2661915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4665224) q[2];
sx q[2];
rz(-1.0202946) q[2];
sx q[2];
rz(2.7543219) q[2];
rz(-1.1747423) q[3];
sx q[3];
rz(-1.6956804) q[3];
sx q[3];
rz(0.89850473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77884498) q[0];
sx q[0];
rz(-1.7411106) q[0];
sx q[0];
rz(-1.8722906) q[0];
rz(1.6954039) q[1];
sx q[1];
rz(-1.8480999) q[1];
sx q[1];
rz(-0.92099014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5808153) q[0];
sx q[0];
rz(-1.2755503) q[0];
sx q[0];
rz(-0.92714374) q[0];
rz(-pi) q[1];
rz(0.66514449) q[2];
sx q[2];
rz(-1.4262305) q[2];
sx q[2];
rz(-1.7139561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.16352902) q[1];
sx q[1];
rz(-2.7394751) q[1];
sx q[1];
rz(1.0671713) q[1];
rz(-pi) q[2];
rz(-1.5637573) q[3];
sx q[3];
rz(-1.3448997) q[3];
sx q[3];
rz(1.5510538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15128073) q[2];
sx q[2];
rz(-1.9419779) q[2];
sx q[2];
rz(-2.3178103) q[2];
rz(1.6173877) q[3];
sx q[3];
rz(-1.5533605) q[3];
sx q[3];
rz(-1.2113781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23564944) q[0];
sx q[0];
rz(-0.88053572) q[0];
sx q[0];
rz(2.549951) q[0];
rz(-0.94332424) q[1];
sx q[1];
rz(-2.0843518) q[1];
sx q[1];
rz(-1.0234157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9236301) q[0];
sx q[0];
rz(-1.4537088) q[0];
sx q[0];
rz(0.10770672) q[0];
x q[1];
rz(2.1322971) q[2];
sx q[2];
rz(-0.63234419) q[2];
sx q[2];
rz(2.655766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4893787) q[1];
sx q[1];
rz(-1.4783597) q[1];
sx q[1];
rz(-3.1104065) q[1];
x q[2];
rz(-2.4571506) q[3];
sx q[3];
rz(-1.3108284) q[3];
sx q[3];
rz(0.79567676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8817875) q[2];
sx q[2];
rz(-1.2299579) q[2];
sx q[2];
rz(-1.7477431) q[2];
rz(1.7143837) q[3];
sx q[3];
rz(-2.2673159) q[3];
sx q[3];
rz(-0.29903308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-2.3706936) q[0];
sx q[0];
rz(-0.40632668) q[0];
sx q[0];
rz(0.64737493) q[0];
rz(0.24230832) q[1];
sx q[1];
rz(-1.4360042) q[1];
sx q[1];
rz(-1.4586331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74170106) q[0];
sx q[0];
rz(-1.8938602) q[0];
sx q[0];
rz(1.7123651) q[0];
rz(-2.2307441) q[2];
sx q[2];
rz(-2.0767412) q[2];
sx q[2];
rz(-2.8252223) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43964115) q[1];
sx q[1];
rz(-2.6547013) q[1];
sx q[1];
rz(1.6488579) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3620831) q[3];
sx q[3];
rz(-1.0406896) q[3];
sx q[3];
rz(0.63941832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.709939) q[2];
sx q[2];
rz(-0.84448758) q[2];
sx q[2];
rz(-1.2246152) q[2];
rz(2.8905408) q[3];
sx q[3];
rz(-0.71728388) q[3];
sx q[3];
rz(-2.203598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5421903) q[0];
sx q[0];
rz(-2.0090065) q[0];
sx q[0];
rz(-2.9436924) q[0];
rz(1.9056162) q[1];
sx q[1];
rz(-0.37207741) q[1];
sx q[1];
rz(3.1370251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76659996) q[0];
sx q[0];
rz(-0.99762929) q[0];
sx q[0];
rz(-1.8660924) q[0];
rz(-1.7103689) q[2];
sx q[2];
rz(-2.6032631) q[2];
sx q[2];
rz(-2.0135422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1612322) q[1];
sx q[1];
rz(-1.8440108) q[1];
sx q[1];
rz(1.4153764) q[1];
rz(-0.35085939) q[3];
sx q[3];
rz(-1.7688732) q[3];
sx q[3];
rz(-0.7367062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4602451) q[2];
sx q[2];
rz(-2.6474417) q[2];
sx q[2];
rz(0.34026185) q[2];
rz(1.5334689) q[3];
sx q[3];
rz(-1.3932649) q[3];
sx q[3];
rz(1.6218012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.915864) q[0];
sx q[0];
rz(-0.69750834) q[0];
sx q[0];
rz(-1.1874143) q[0];
rz(-1.7200708) q[1];
sx q[1];
rz(-2.7233796) q[1];
sx q[1];
rz(0.13030599) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73459816) q[0];
sx q[0];
rz(-0.13665527) q[0];
sx q[0];
rz(-2.0137441) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58660686) q[2];
sx q[2];
rz(-0.97525037) q[2];
sx q[2];
rz(-1.4968703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9355311) q[1];
sx q[1];
rz(-1.5714679) q[1];
sx q[1];
rz(2.9737581) q[1];
x q[2];
rz(3.1316787) q[3];
sx q[3];
rz(-1.8335206) q[3];
sx q[3];
rz(-0.81068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80798951) q[2];
sx q[2];
rz(-1.5895867) q[2];
sx q[2];
rz(2.1017334) q[2];
rz(2.9676159) q[3];
sx q[3];
rz(-1.1287374) q[3];
sx q[3];
rz(-2.4719293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60798764) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(-0.84130353) q[0];
rz(-1.816642) q[1];
sx q[1];
rz(-2.1957896) q[1];
sx q[1];
rz(-1.673505) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45638915) q[0];
sx q[0];
rz(-0.51420553) q[0];
sx q[0];
rz(0.77001621) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5910452) q[2];
sx q[2];
rz(-0.51873461) q[2];
sx q[2];
rz(0.27257365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7699112) q[1];
sx q[1];
rz(-1.103319) q[1];
sx q[1];
rz(-0.63422306) q[1];
rz(-pi) q[2];
rz(-2.2696247) q[3];
sx q[3];
rz(-2.6879394) q[3];
sx q[3];
rz(3.0986971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2738721) q[2];
sx q[2];
rz(-2.2404859) q[2];
sx q[2];
rz(-0.64485288) q[2];
rz(-1.0817179) q[3];
sx q[3];
rz(-2.7223301) q[3];
sx q[3];
rz(2.4206415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69407392) q[0];
sx q[0];
rz(-2.9677291) q[0];
sx q[0];
rz(-0.41282594) q[0];
rz(1.6089571) q[1];
sx q[1];
rz(-0.91333476) q[1];
sx q[1];
rz(0.70721165) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.083867) q[0];
sx q[0];
rz(-2.3207156) q[0];
sx q[0];
rz(2.7146401) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17036713) q[2];
sx q[2];
rz(-2.3922386) q[2];
sx q[2];
rz(-3.0339277) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3601348) q[1];
sx q[1];
rz(-2.4943922) q[1];
sx q[1];
rz(2.5639046) q[1];
rz(-pi) q[2];
rz(1.8052638) q[3];
sx q[3];
rz(-0.55113652) q[3];
sx q[3];
rz(-2.7296327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1327208) q[2];
sx q[2];
rz(-2.8421695) q[2];
sx q[2];
rz(0.51187619) q[2];
rz(0.68073186) q[3];
sx q[3];
rz(-1.778089) q[3];
sx q[3];
rz(2.9752922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39636382) q[0];
sx q[0];
rz(-1.5928716) q[0];
sx q[0];
rz(2.6863099) q[0];
rz(2.0782754) q[1];
sx q[1];
rz(-0.89439193) q[1];
sx q[1];
rz(-0.071970073) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4819099) q[0];
sx q[0];
rz(-1.558466) q[0];
sx q[0];
rz(3.1083049) q[0];
rz(0.44092245) q[2];
sx q[2];
rz(-1.5498232) q[2];
sx q[2];
rz(-1.4628589) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1942206) q[1];
sx q[1];
rz(-0.76057077) q[1];
sx q[1];
rz(2.0921699) q[1];
rz(0.91225635) q[3];
sx q[3];
rz(-2.7317762) q[3];
sx q[3];
rz(1.8456822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7887743) q[2];
sx q[2];
rz(-0.41663751) q[2];
sx q[2];
rz(-2.519506) q[2];
rz(-0.82428733) q[3];
sx q[3];
rz(-1.0930748) q[3];
sx q[3];
rz(-0.85339671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1085827) q[0];
sx q[0];
rz(-2.0350631) q[0];
sx q[0];
rz(-2.1906817) q[0];
rz(1.7225601) q[1];
sx q[1];
rz(-2.6585237) q[1];
sx q[1];
rz(-2.5249265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86980235) q[0];
sx q[0];
rz(-1.5464725) q[0];
sx q[0];
rz(0.64016827) q[0];
rz(-pi) q[1];
rz(0.51787106) q[2];
sx q[2];
rz(-2.4642337) q[2];
sx q[2];
rz(0.022938722) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93393713) q[1];
sx q[1];
rz(-0.21156921) q[1];
sx q[1];
rz(-3.1283698) q[1];
rz(-pi) q[2];
rz(-2.4973013) q[3];
sx q[3];
rz(-2.3228495) q[3];
sx q[3];
rz(-0.14547023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0509433) q[2];
sx q[2];
rz(-1.1530777) q[2];
sx q[2];
rz(-2.0743745) q[2];
rz(-1.0956592) q[3];
sx q[3];
rz(-1.9039543) q[3];
sx q[3];
rz(0.9637951) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6296366) q[0];
sx q[0];
rz(-0.98573276) q[0];
sx q[0];
rz(-1.0712256) q[0];
rz(-1.4818954) q[1];
sx q[1];
rz(-2.0493458) q[1];
sx q[1];
rz(0.58072166) q[1];
rz(0.26255519) q[2];
sx q[2];
rz(-1.5947744) q[2];
sx q[2];
rz(-2.6626432) q[2];
rz(-0.4023424) q[3];
sx q[3];
rz(-2.8035223) q[3];
sx q[3];
rz(0.30717472) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
