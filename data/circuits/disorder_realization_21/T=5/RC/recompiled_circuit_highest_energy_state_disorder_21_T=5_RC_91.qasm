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
rz(-3.1132071) q[0];
sx q[0];
rz(-0.95771587) q[0];
rz(-2.7081642) q[1];
sx q[1];
rz(1.537701) q[1];
sx q[1];
rz(8.796506) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0032438) q[0];
sx q[0];
rz(-1.4491557) q[0];
sx q[0];
rz(-0.37183372) q[0];
rz(0.84282173) q[2];
sx q[2];
rz(-0.6302689) q[2];
sx q[2];
rz(-2.3284019) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0121668) q[1];
sx q[1];
rz(-0.91246997) q[1];
sx q[1];
rz(0.40947394) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1584627) q[3];
sx q[3];
rz(-1.6618007) q[3];
sx q[3];
rz(-3.0418398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6750703) q[2];
sx q[2];
rz(-2.1212981) q[2];
sx q[2];
rz(0.38727078) q[2];
rz(1.9668503) q[3];
sx q[3];
rz(-1.4459123) q[3];
sx q[3];
rz(-0.89850473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627477) q[0];
sx q[0];
rz(-1.400482) q[0];
sx q[0];
rz(-1.269302) q[0];
rz(-1.4461888) q[1];
sx q[1];
rz(-1.8480999) q[1];
sx q[1];
rz(-0.92099014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816278) q[0];
sx q[0];
rz(-2.4422997) q[0];
sx q[0];
rz(-1.101732) q[0];
x q[1];
rz(2.4764482) q[2];
sx q[2];
rz(-1.7153622) q[2];
sx q[2];
rz(1.4276366) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.16352902) q[1];
sx q[1];
rz(-0.40211758) q[1];
sx q[1];
rz(-2.0744214) q[1];
x q[2];
rz(-1.5778353) q[3];
sx q[3];
rz(-1.796693) q[3];
sx q[3];
rz(1.5510538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15128073) q[2];
sx q[2];
rz(-1.9419779) q[2];
sx q[2];
rz(-0.82378236) q[2];
rz(-1.524205) q[3];
sx q[3];
rz(-1.5882322) q[3];
sx q[3];
rz(1.2113781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9059432) q[0];
sx q[0];
rz(-0.88053572) q[0];
sx q[0];
rz(-0.5916416) q[0];
rz(-0.94332424) q[1];
sx q[1];
rz(-1.0572409) q[1];
sx q[1];
rz(1.0234157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34020326) q[0];
sx q[0];
rz(-1.4638299) q[0];
sx q[0];
rz(1.68856) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7696174) q[2];
sx q[2];
rz(-2.0947273) q[2];
sx q[2];
rz(-1.1481448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4893787) q[1];
sx q[1];
rz(-1.6632329) q[1];
sx q[1];
rz(-3.1104065) q[1];
x q[2];
rz(2.4571506) q[3];
sx q[3];
rz(-1.3108284) q[3];
sx q[3];
rz(2.3459159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8817875) q[2];
sx q[2];
rz(-1.2299579) q[2];
sx q[2];
rz(1.7477431) q[2];
rz(-1.7143837) q[3];
sx q[3];
rz(-2.2673159) q[3];
sx q[3];
rz(-2.8425596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77089906) q[0];
sx q[0];
rz(-2.735266) q[0];
sx q[0];
rz(-2.4942177) q[0];
rz(2.8992843) q[1];
sx q[1];
rz(-1.7055885) q[1];
sx q[1];
rz(-1.4586331) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9779309) q[0];
sx q[0];
rz(-2.7898698) q[0];
sx q[0];
rz(2.7427384) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2307441) q[2];
sx q[2];
rz(-1.0648515) q[2];
sx q[2];
rz(-2.8252223) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.200176) q[1];
sx q[1];
rz(-1.5343018) q[1];
sx q[1];
rz(-1.0851651) q[1];
rz(0.8816054) q[3];
sx q[3];
rz(-0.91920105) q[3];
sx q[3];
rz(-2.6736265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4316537) q[2];
sx q[2];
rz(-2.2971051) q[2];
sx q[2];
rz(-1.9169774) q[2];
rz(-2.8905408) q[3];
sx q[3];
rz(-0.71728388) q[3];
sx q[3];
rz(2.203598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5421903) q[0];
sx q[0];
rz(-2.0090065) q[0];
sx q[0];
rz(2.9436924) q[0];
rz(1.9056162) q[1];
sx q[1];
rz(-0.37207741) q[1];
sx q[1];
rz(-0.0045675357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5008847) q[0];
sx q[0];
rz(-1.3237778) q[0];
sx q[0];
rz(2.5480854) q[0];
rz(-pi) q[1];
rz(1.4312237) q[2];
sx q[2];
rz(-2.6032631) q[2];
sx q[2];
rz(-2.0135422) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45426863) q[1];
sx q[1];
rz(-0.31336323) q[1];
sx q[1];
rz(-0.50468495) q[1];
x q[2];
rz(2.7907333) q[3];
sx q[3];
rz(-1.3727194) q[3];
sx q[3];
rz(0.7367062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4602451) q[2];
sx q[2];
rz(-0.494151) q[2];
sx q[2];
rz(-0.34026185) q[2];
rz(-1.5334689) q[3];
sx q[3];
rz(-1.3932649) q[3];
sx q[3];
rz(-1.6218012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.915864) q[0];
sx q[0];
rz(-2.4440843) q[0];
sx q[0];
rz(1.9541784) q[0];
rz(1.7200708) q[1];
sx q[1];
rz(-2.7233796) q[1];
sx q[1];
rz(3.0112867) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4069945) q[0];
sx q[0];
rz(-0.13665527) q[0];
sx q[0];
rz(-1.1278485) q[0];
rz(-2.2557507) q[2];
sx q[2];
rz(-0.81007137) q[2];
sx q[2];
rz(-2.5145234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3686982) q[1];
sx q[1];
rz(-0.16783585) q[1];
sx q[1];
rz(-0.0040199587) q[1];
rz(-pi) q[2];
rz(3.1316787) q[3];
sx q[3];
rz(-1.3080721) q[3];
sx q[3];
rz(0.81068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80798951) q[2];
sx q[2];
rz(-1.5895867) q[2];
sx q[2];
rz(1.0398593) q[2];
rz(0.17397675) q[3];
sx q[3];
rz(-2.0128553) q[3];
sx q[3];
rz(-2.4719293) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.533605) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(0.84130353) q[0];
rz(1.3249506) q[1];
sx q[1];
rz(-2.1957896) q[1];
sx q[1];
rz(-1.673505) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3260088) q[0];
sx q[0];
rz(-1.9202613) q[0];
sx q[0];
rz(-2.7563226) q[0];
x q[1];
rz(1.05215) q[2];
sx q[2];
rz(-1.5607578) q[2];
sx q[2];
rz(-1.280637) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.12097479) q[1];
sx q[1];
rz(-1.0134032) q[1];
sx q[1];
rz(-1.0110028) q[1];
x q[2];
rz(-1.928059) q[3];
sx q[3];
rz(-1.284984) q[3];
sx q[3];
rz(2.2605726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86772052) q[2];
sx q[2];
rz(-0.90110675) q[2];
sx q[2];
rz(0.64485288) q[2];
rz(1.0817179) q[3];
sx q[3];
rz(-2.7223301) q[3];
sx q[3];
rz(-2.4206415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69407392) q[0];
sx q[0];
rz(-2.9677291) q[0];
sx q[0];
rz(2.7287667) q[0];
rz(-1.5326356) q[1];
sx q[1];
rz(-0.91333476) q[1];
sx q[1];
rz(-2.434381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0577257) q[0];
sx q[0];
rz(-2.3207156) q[0];
sx q[0];
rz(-0.42695257) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17036713) q[2];
sx q[2];
rz(-0.74935407) q[2];
sx q[2];
rz(-0.10766497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.78145786) q[1];
sx q[1];
rz(-2.4943922) q[1];
sx q[1];
rz(-2.5639046) q[1];
rz(-pi) q[2];
rz(-2.1096346) q[3];
sx q[3];
rz(-1.4488359) q[3];
sx q[3];
rz(-1.7820047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0088719) q[2];
sx q[2];
rz(-2.8421695) q[2];
sx q[2];
rz(-2.6297165) q[2];
rz(2.4608608) q[3];
sx q[3];
rz(-1.778089) q[3];
sx q[3];
rz(-2.9752922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39636382) q[0];
sx q[0];
rz(-1.5487211) q[0];
sx q[0];
rz(0.45528278) q[0];
rz(2.0782754) q[1];
sx q[1];
rz(-2.2472007) q[1];
sx q[1];
rz(0.071970073) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44350478) q[0];
sx q[0];
rz(-3.1060954) q[0];
sx q[0];
rz(-2.7867691) q[0];
rz(-pi) q[1];
rz(2.7006702) q[2];
sx q[2];
rz(-1.5917695) q[2];
sx q[2];
rz(1.6787337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2769952) q[1];
sx q[1];
rz(-2.2114843) q[1];
sx q[1];
rz(-0.44257852) q[1];
rz(-pi) q[2];
rz(-0.25983475) q[3];
sx q[3];
rz(-1.2502115) q[3];
sx q[3];
rz(1.1449006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7887743) q[2];
sx q[2];
rz(-2.7249551) q[2];
sx q[2];
rz(2.519506) q[2];
rz(-0.82428733) q[3];
sx q[3];
rz(-2.0485179) q[3];
sx q[3];
rz(-2.2881959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033009919) q[0];
sx q[0];
rz(-1.1065296) q[0];
sx q[0];
rz(-0.95091096) q[0];
rz(-1.4190326) q[1];
sx q[1];
rz(-0.483069) q[1];
sx q[1];
rz(2.5249265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717903) q[0];
sx q[0];
rz(-1.5464725) q[0];
sx q[0];
rz(-2.5014244) q[0];
rz(-pi) q[1];
rz(2.6237216) q[2];
sx q[2];
rz(-0.67735891) q[2];
sx q[2];
rz(-3.1186539) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.93393713) q[1];
sx q[1];
rz(-0.21156921) q[1];
sx q[1];
rz(3.1283698) q[1];
rz(-2.1415878) q[3];
sx q[3];
rz(-2.1943008) q[3];
sx q[3];
rz(0.68736651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0906494) q[2];
sx q[2];
rz(-1.9885149) q[2];
sx q[2];
rz(2.0743745) q[2];
rz(-1.0956592) q[3];
sx q[3];
rz(-1.2376384) q[3];
sx q[3];
rz(2.1777976) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6296366) q[0];
sx q[0];
rz(-0.98573276) q[0];
sx q[0];
rz(-1.0712256) q[0];
rz(1.4818954) q[1];
sx q[1];
rz(-1.0922468) q[1];
sx q[1];
rz(-2.560871) q[1];
rz(-2.8790375) q[2];
sx q[2];
rz(-1.5947744) q[2];
sx q[2];
rz(-2.6626432) q[2];
rz(-2.7392503) q[3];
sx q[3];
rz(-0.33807031) q[3];
sx q[3];
rz(-2.8344179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
