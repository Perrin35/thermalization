OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(0.983239) q[1];
sx q[1];
rz(2.6020738) q[1];
sx q[1];
rz(10.625216) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96785986) q[0];
sx q[0];
rz(-1.9779441) q[0];
sx q[0];
rz(-1.3074387) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0317694) q[2];
sx q[2];
rz(-1.6845778) q[2];
sx q[2];
rz(-0.10345085) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28509724) q[1];
sx q[1];
rz(-0.76438099) q[1];
sx q[1];
rz(0.38715036) q[1];
rz(-pi) q[2];
rz(-0.31906268) q[3];
sx q[3];
rz(-2.3205119) q[3];
sx q[3];
rz(2.0882437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0212705) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(1.5585287) q[2];
rz(-2.1448686) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(-0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834171) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(3.0875207) q[0];
rz(1.1955098) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(2.6057459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0756404) q[0];
sx q[0];
rz(-1.4464805) q[0];
sx q[0];
rz(2.1527704) q[0];
rz(-pi) q[1];
rz(0.29157721) q[2];
sx q[2];
rz(-0.73086408) q[2];
sx q[2];
rz(2.4183194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0332196) q[1];
sx q[1];
rz(-2.6544826) q[1];
sx q[1];
rz(0.56652041) q[1];
rz(1.348098) q[3];
sx q[3];
rz(-0.14651146) q[3];
sx q[3];
rz(0.033586249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(3.0751394) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(-2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41985837) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(0.15047519) q[0];
rz(-2.6843605) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(-3.1157852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3993527) q[0];
sx q[0];
rz(-2.8909546) q[0];
sx q[0];
rz(1.6514732) q[0];
rz(-pi) q[1];
rz(1.843156) q[2];
sx q[2];
rz(-2.8189427) q[2];
sx q[2];
rz(1.9793561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7681231) q[1];
sx q[1];
rz(-2.3489531) q[1];
sx q[1];
rz(0.5975238) q[1];
rz(-3.0098626) q[3];
sx q[3];
rz(-1.9506491) q[3];
sx q[3];
rz(-2.6704138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(0.5853816) q[2];
rz(-2.9600926) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90081763) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(0.88090849) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-0.53612971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25585184) q[0];
sx q[0];
rz(-0.8937853) q[0];
sx q[0];
rz(2.6141502) q[0];
rz(-pi) q[1];
rz(0.23303194) q[2];
sx q[2];
rz(-0.2430025) q[2];
sx q[2];
rz(-0.88569966) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0587479) q[1];
sx q[1];
rz(-1.1943598) q[1];
sx q[1];
rz(2.46545) q[1];
rz(-0.36066182) q[3];
sx q[3];
rz(-2.4324527) q[3];
sx q[3];
rz(-3.0879471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46999103) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(-0.70703834) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(-2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.2351284) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(1.0513603) q[0];
rz(1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-3.0984745) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2152104) q[0];
sx q[0];
rz(-0.97391093) q[0];
sx q[0];
rz(-0.2290639) q[0];
rz(-pi) q[1];
rz(-2.0448858) q[2];
sx q[2];
rz(-0.57069639) q[2];
sx q[2];
rz(1.2793465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4823618) q[1];
sx q[1];
rz(-1.2995969) q[1];
sx q[1];
rz(-1.8455452) q[1];
x q[2];
rz(1.4210012) q[3];
sx q[3];
rz(-2.3547223) q[3];
sx q[3];
rz(-2.8926135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0126426) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(-0.35219231) q[2];
rz(2.5514065) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(-0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234289) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(2.39134) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(1.0587143) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569865) q[0];
sx q[0];
rz(-1.3625506) q[0];
sx q[0];
rz(-0.44021846) q[0];
x q[1];
rz(1.9361587) q[2];
sx q[2];
rz(-2.0632671) q[2];
sx q[2];
rz(-1.1577215) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8844879) q[1];
sx q[1];
rz(-0.37933644) q[1];
sx q[1];
rz(0.92909716) q[1];
rz(-2.9162354) q[3];
sx q[3];
rz(-0.72031027) q[3];
sx q[3];
rz(-2.5845137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6713312) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(-1.9227825) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0999775) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(-1.6301427) q[0];
rz(-1.3776243) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(-2.2999433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8182897) q[0];
sx q[0];
rz(-2.0879732) q[0];
sx q[0];
rz(-1.3099758) q[0];
rz(-pi) q[1];
rz(-1.3782578) q[2];
sx q[2];
rz(-0.95270573) q[2];
sx q[2];
rz(0.75380177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48982606) q[1];
sx q[1];
rz(-3.0511599) q[1];
sx q[1];
rz(-1.0033146) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8089201) q[3];
sx q[3];
rz(-2.3861285) q[3];
sx q[3];
rz(1.2681703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(3.0916396) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426303) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(1.4021953) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-2.7239674) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0324875) q[0];
sx q[0];
rz(-2.0628477) q[0];
sx q[0];
rz(-1.0743272) q[0];
rz(-pi) q[1];
rz(-1.9365963) q[2];
sx q[2];
rz(-1.6779643) q[2];
sx q[2];
rz(2.0049713) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8466134) q[1];
sx q[1];
rz(-1.1049005) q[1];
sx q[1];
rz(0.93211517) q[1];
rz(-pi) q[2];
rz(-2.8475259) q[3];
sx q[3];
rz(-1.2122279) q[3];
sx q[3];
rz(-0.29953526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3461356) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(1.0428838) q[2];
rz(-2.4677094) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(0.90014443) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(-0.62966627) q[0];
rz(0.57485238) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(-2.1946857) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7787951) q[0];
sx q[0];
rz(-0.4501833) q[0];
sx q[0];
rz(-2.3953526) q[0];
x q[1];
rz(1.6106748) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(-3.1181042) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35518256) q[1];
sx q[1];
rz(-1.2372412) q[1];
sx q[1];
rz(1.4817609) q[1];
x q[2];
rz(0.19091786) q[3];
sx q[3];
rz(-2.0058161) q[3];
sx q[3];
rz(-2.2866979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(1.8927195) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(0.26556382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749851) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(-0.65504909) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5376741) q[0];
sx q[0];
rz(-0.18225056) q[0];
sx q[0];
rz(-2.8331579) q[0];
x q[1];
rz(0.23004736) q[2];
sx q[2];
rz(-0.86798475) q[2];
sx q[2];
rz(-0.41781296) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6380784) q[1];
sx q[1];
rz(-2.6942309) q[1];
sx q[1];
rz(2.4380986) q[1];
x q[2];
rz(-0.44390042) q[3];
sx q[3];
rz(-1.1392986) q[3];
sx q[3];
rz(-2.3865226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3506938) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(2.541686) q[2];
rz(0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(2.9121493) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(3.0145666) q[2];
sx q[2];
rz(-2.1913678) q[2];
sx q[2];
rz(0.36507228) q[2];
rz(-2.2453528) q[3];
sx q[3];
rz(-1.2590209) q[3];
sx q[3];
rz(-0.38929064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
