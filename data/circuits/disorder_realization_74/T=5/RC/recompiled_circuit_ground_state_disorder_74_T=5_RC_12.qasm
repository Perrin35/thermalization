OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5689019) q[0];
sx q[0];
rz(-0.98839086) q[0];
sx q[0];
rz(2.6925777) q[0];
rz(2.9721337) q[1];
sx q[1];
rz(-0.13154498) q[1];
sx q[1];
rz(2.0102672) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5857475) q[0];
sx q[0];
rz(-2.1679584) q[0];
sx q[0];
rz(-0.35564977) q[0];
rz(-pi) q[1];
rz(-2.5070081) q[2];
sx q[2];
rz(-2.0014844) q[2];
sx q[2];
rz(0.3542977) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1234963) q[1];
sx q[1];
rz(-0.61401788) q[1];
sx q[1];
rz(-1.2662925) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2901957) q[3];
sx q[3];
rz(-0.56125703) q[3];
sx q[3];
rz(-0.11229501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.955287) q[2];
sx q[2];
rz(-0.39262843) q[2];
sx q[2];
rz(0.10360959) q[2];
rz(-2.7747532) q[3];
sx q[3];
rz(-1.47374) q[3];
sx q[3];
rz(1.8240671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1845301) q[0];
sx q[0];
rz(-2.6657031) q[0];
sx q[0];
rz(-2.6153508) q[0];
rz(-1.8980252) q[1];
sx q[1];
rz(-1.7310111) q[1];
sx q[1];
rz(-0.19613656) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84700245) q[0];
sx q[0];
rz(-0.46200141) q[0];
sx q[0];
rz(-1.4205971) q[0];
rz(-0.7095269) q[2];
sx q[2];
rz(-1.4440184) q[2];
sx q[2];
rz(1.7293872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1795219) q[1];
sx q[1];
rz(-1.7129494) q[1];
sx q[1];
rz(-1.1198055) q[1];
rz(-pi) q[2];
rz(-0.72153458) q[3];
sx q[3];
rz(-0.01577687) q[3];
sx q[3];
rz(2.0023014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5474995) q[2];
sx q[2];
rz(-0.27442351) q[2];
sx q[2];
rz(-0.37626949) q[2];
rz(-2.2606692) q[3];
sx q[3];
rz(-1.8062402) q[3];
sx q[3];
rz(-2.3295565) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6905717) q[0];
sx q[0];
rz(-2.8132827) q[0];
sx q[0];
rz(2.1300533) q[0];
rz(2.4729589) q[1];
sx q[1];
rz(-1.9790383) q[1];
sx q[1];
rz(2.5943601) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035346383) q[0];
sx q[0];
rz(-1.5515455) q[0];
sx q[0];
rz(-0.20142844) q[0];
x q[1];
rz(0.34788068) q[2];
sx q[2];
rz(-0.82150412) q[2];
sx q[2];
rz(2.2583928) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95067945) q[1];
sx q[1];
rz(-2.4362794) q[1];
sx q[1];
rz(-1.1937792) q[1];
rz(-2.9487398) q[3];
sx q[3];
rz(-0.82369419) q[3];
sx q[3];
rz(0.76092645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4554567) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(-1.6395462) q[2];
rz(-2.5655668) q[3];
sx q[3];
rz(-1.6582158) q[3];
sx q[3];
rz(-0.36147931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6308924) q[0];
sx q[0];
rz(-0.80131131) q[0];
sx q[0];
rz(-2.8019688) q[0];
rz(-0.54061186) q[1];
sx q[1];
rz(-2.4410591) q[1];
sx q[1];
rz(-0.9300173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.697111) q[0];
sx q[0];
rz(-0.15938317) q[0];
sx q[0];
rz(2.0305994) q[0];
x q[1];
rz(2.5581237) q[2];
sx q[2];
rz(-1.8975048) q[2];
sx q[2];
rz(-2.938156) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8992638) q[1];
sx q[1];
rz(-1.5643969) q[1];
sx q[1];
rz(1.0914299) q[1];
rz(2.2100979) q[3];
sx q[3];
rz(-2.3911797) q[3];
sx q[3];
rz(1.9281333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0129619) q[2];
sx q[2];
rz(-1.4350812) q[2];
sx q[2];
rz(-1.5677412) q[2];
rz(-0.74357998) q[3];
sx q[3];
rz(-0.79135528) q[3];
sx q[3];
rz(-2.1775406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6745233) q[0];
sx q[0];
rz(-2.2860797) q[0];
sx q[0];
rz(-0.056644406) q[0];
rz(-1.665834) q[1];
sx q[1];
rz(-2.00878) q[1];
sx q[1];
rz(1.7830361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751793) q[0];
sx q[0];
rz(-0.73690276) q[0];
sx q[0];
rz(2.3794258) q[0];
rz(-pi) q[1];
rz(-2.4841294) q[2];
sx q[2];
rz(-1.7320447) q[2];
sx q[2];
rz(-0.25875124) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1093692) q[1];
sx q[1];
rz(-2.0521161) q[1];
sx q[1];
rz(-0.66415031) q[1];
x q[2];
rz(0.56157063) q[3];
sx q[3];
rz(-1.4709146) q[3];
sx q[3];
rz(2.9390854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8654827) q[2];
sx q[2];
rz(-1.411974) q[2];
sx q[2];
rz(1.126368) q[2];
rz(-1.3364835) q[3];
sx q[3];
rz(-2.0650605) q[3];
sx q[3];
rz(2.0520463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4514076) q[0];
sx q[0];
rz(-1.0253588) q[0];
sx q[0];
rz(-3.0080646) q[0];
rz(-2.1639157) q[1];
sx q[1];
rz(-0.97594273) q[1];
sx q[1];
rz(2.8547063) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2539352) q[0];
sx q[0];
rz(-2.2713221) q[0];
sx q[0];
rz(-3.0791686) q[0];
rz(-pi) q[1];
rz(1.6966218) q[2];
sx q[2];
rz(-0.33060961) q[2];
sx q[2];
rz(2.5839407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0930909) q[1];
sx q[1];
rz(-1.7939006) q[1];
sx q[1];
rz(0.99324171) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.828015) q[3];
sx q[3];
rz(-0.99204274) q[3];
sx q[3];
rz(-0.46427765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7202683) q[2];
sx q[2];
rz(-1.4807533) q[2];
sx q[2];
rz(1.486091) q[2];
rz(0.36763516) q[3];
sx q[3];
rz(-2.1143819) q[3];
sx q[3];
rz(-2.0462842) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7169645) q[0];
sx q[0];
rz(-2.3724738) q[0];
sx q[0];
rz(-2.6677483) q[0];
rz(2.4773856) q[1];
sx q[1];
rz(-1.217548) q[1];
sx q[1];
rz(-1.7582105) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5067399) q[0];
sx q[0];
rz(-2.228745) q[0];
sx q[0];
rz(0.94292504) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1138713) q[2];
sx q[2];
rz(-1.1528735) q[2];
sx q[2];
rz(2.6270285) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5210628) q[1];
sx q[1];
rz(-2.8179666) q[1];
sx q[1];
rz(-0.77038295) q[1];
x q[2];
rz(-0.17857213) q[3];
sx q[3];
rz(-1.7925486) q[3];
sx q[3];
rz(2.5933883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15709269) q[2];
sx q[2];
rz(-1.7326771) q[2];
sx q[2];
rz(0.3248997) q[2];
rz(0.38069185) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(-2.0438173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070504524) q[0];
sx q[0];
rz(-2.4779713) q[0];
sx q[0];
rz(2.2027503) q[0];
rz(-2.4414869) q[1];
sx q[1];
rz(-2.5036) q[1];
sx q[1];
rz(-0.28900388) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2719921) q[0];
sx q[0];
rz(-1.648252) q[0];
sx q[0];
rz(-1.3415706) q[0];
rz(-pi) q[1];
rz(0.49792011) q[2];
sx q[2];
rz(-0.98518956) q[2];
sx q[2];
rz(-2.917054) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4314193) q[1];
sx q[1];
rz(-2.0906716) q[1];
sx q[1];
rz(-0.10837491) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38095946) q[3];
sx q[3];
rz(-1.5039615) q[3];
sx q[3];
rz(1.9240954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2603904) q[2];
sx q[2];
rz(-1.6198879) q[2];
sx q[2];
rz(2.728906) q[2];
rz(0.9203426) q[3];
sx q[3];
rz(-2.6350382) q[3];
sx q[3];
rz(-1.9119561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6270139) q[0];
sx q[0];
rz(-0.98015061) q[0];
sx q[0];
rz(2.8421616) q[0];
rz(-2.0191655) q[1];
sx q[1];
rz(-1.9405245) q[1];
sx q[1];
rz(1.0312414) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.497488) q[0];
sx q[0];
rz(-1.9369164) q[0];
sx q[0];
rz(-0.011542277) q[0];
rz(-1.769936) q[2];
sx q[2];
rz(-2.0250892) q[2];
sx q[2];
rz(0.21541883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.448062) q[1];
sx q[1];
rz(-1.8603357) q[1];
sx q[1];
rz(0.84622835) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93017756) q[3];
sx q[3];
rz(-1.0856461) q[3];
sx q[3];
rz(-1.5826011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1624182) q[2];
sx q[2];
rz(-1.3434429) q[2];
sx q[2];
rz(0.24270414) q[2];
rz(1.3001214) q[3];
sx q[3];
rz(-2.0827115) q[3];
sx q[3];
rz(0.21568957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31371394) q[0];
sx q[0];
rz(-0.21610459) q[0];
sx q[0];
rz(-1.4208273) q[0];
rz(2.8315262) q[1];
sx q[1];
rz(-2.2646751) q[1];
sx q[1];
rz(-1.5440595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89429796) q[0];
sx q[0];
rz(-1.9820807) q[0];
sx q[0];
rz(2.5687587) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9126911) q[2];
sx q[2];
rz(-1.4871305) q[2];
sx q[2];
rz(-2.8852706) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3812755) q[1];
sx q[1];
rz(-0.3418498) q[1];
sx q[1];
rz(2.443497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4214823) q[3];
sx q[3];
rz(-2.930899) q[3];
sx q[3];
rz(-1.9062476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0073504) q[2];
sx q[2];
rz(-1.4126567) q[2];
sx q[2];
rz(-1.7751815) q[2];
rz(-2.2640696) q[3];
sx q[3];
rz(-1.4656504) q[3];
sx q[3];
rz(0.88959488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14961814) q[0];
sx q[0];
rz(-1.5562973) q[0];
sx q[0];
rz(1.656116) q[0];
rz(-0.86434518) q[1];
sx q[1];
rz(-2.1280011) q[1];
sx q[1];
rz(-0.43307532) q[1];
rz(-0.36467411) q[2];
sx q[2];
rz(-0.52959792) q[2];
sx q[2];
rz(-2.9858225) q[2];
rz(1.4644365) q[3];
sx q[3];
rz(-0.48135664) q[3];
sx q[3];
rz(2.4873747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
