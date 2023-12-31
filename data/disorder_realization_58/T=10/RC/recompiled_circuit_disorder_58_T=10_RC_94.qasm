OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(0.82011861) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(-2.3109205) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18618628) q[0];
sx q[0];
rz(-0.49476981) q[0];
sx q[0];
rz(-2.4525989) q[0];
rz(-pi) q[1];
rz(-2.0482424) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(-0.29104656) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1066061) q[1];
sx q[1];
rz(-0.88284661) q[1];
sx q[1];
rz(-2.7290542) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44736638) q[3];
sx q[3];
rz(-1.9379741) q[3];
sx q[3];
rz(-2.6255053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(-0.1201771) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(-0.78951019) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(2.8149014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625898) q[0];
sx q[0];
rz(-2.0741182) q[0];
sx q[0];
rz(-2.3110564) q[0];
x q[1];
rz(-1.2454883) q[2];
sx q[2];
rz(-2.6763958) q[2];
sx q[2];
rz(-2.543769) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3152299) q[1];
sx q[1];
rz(-0.76347199) q[1];
sx q[1];
rz(2.1748494) q[1];
x q[2];
rz(-2.2585906) q[3];
sx q[3];
rz(-1.5044754) q[3];
sx q[3];
rz(-1.4955213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6313173) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(-0.10989799) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(-1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(-2.0842016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7085416) q[0];
sx q[0];
rz(-1.9110702) q[0];
sx q[0];
rz(0.021854594) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34122841) q[2];
sx q[2];
rz(-1.382302) q[2];
sx q[2];
rz(0.54268062) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96506572) q[1];
sx q[1];
rz(-0.9170734) q[1];
sx q[1];
rz(-2.7214126) q[1];
rz(-pi) q[2];
rz(-0.4226513) q[3];
sx q[3];
rz(-1.395441) q[3];
sx q[3];
rz(0.56817504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.7791629) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-0.16539703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9497313) q[0];
sx q[0];
rz(-2.9182069) q[0];
sx q[0];
rz(2.1641157) q[0];
rz(1.3894765) q[2];
sx q[2];
rz(-1.4099979) q[2];
sx q[2];
rz(0.60418512) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0927825) q[1];
sx q[1];
rz(-2.4895992) q[1];
sx q[1];
rz(-0.55744967) q[1];
x q[2];
rz(-0.27407077) q[3];
sx q[3];
rz(-1.032853) q[3];
sx q[3];
rz(0.9243954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(2.9361434) q[2];
rz(-2.0139587) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797392) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(0.35476312) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(-2.9096471) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35271586) q[0];
sx q[0];
rz(-1.0838325) q[0];
sx q[0];
rz(0.48077521) q[0];
rz(-2.4460692) q[2];
sx q[2];
rz(-1.4624274) q[2];
sx q[2];
rz(-2.7930789) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1272808) q[1];
sx q[1];
rz(-0.31186549) q[1];
sx q[1];
rz(0.77906268) q[1];
rz(-pi) q[2];
rz(0.36851818) q[3];
sx q[3];
rz(-1.5889865) q[3];
sx q[3];
rz(-2.8858678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0971138) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(-0.20275673) q[0];
rz(2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(2.1441377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049012262) q[0];
sx q[0];
rz(-2.2708587) q[0];
sx q[0];
rz(-2.9655365) q[0];
rz(-pi) q[1];
rz(1.1987655) q[2];
sx q[2];
rz(-2.2673006) q[2];
sx q[2];
rz(-0.55559413) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3012078) q[1];
sx q[1];
rz(-1.1494698) q[1];
sx q[1];
rz(-0.2627443) q[1];
rz(1.0783844) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(-1.2315962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9138907) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(2.690199) q[2];
rz(-2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8354427) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(2.5278032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928771) q[0];
sx q[0];
rz(-2.2971417) q[0];
sx q[0];
rz(2.8918299) q[0];
rz(2.8579312) q[2];
sx q[2];
rz(-1.6322871) q[2];
sx q[2];
rz(-2.9279857) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9871414) q[1];
sx q[1];
rz(-3.0220991) q[1];
sx q[1];
rz(-2.6277072) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11717637) q[3];
sx q[3];
rz(-2.0332608) q[3];
sx q[3];
rz(-1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(-1.6237367) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(2.1571295) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7229066) q[0];
sx q[0];
rz(-1.6081928) q[0];
sx q[0];
rz(2.5953369) q[0];
x q[1];
rz(-0.77163561) q[2];
sx q[2];
rz(-1.8688335) q[2];
sx q[2];
rz(0.76921295) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22673785) q[1];
sx q[1];
rz(-1.2054772) q[1];
sx q[1];
rz(1.564333) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6650388) q[3];
sx q[3];
rz(-1.8503975) q[3];
sx q[3];
rz(-2.9001146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70665923) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(0.53517503) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-0.6859268) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(0.92591441) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5092963) q[0];
sx q[0];
rz(-2.8863781) q[0];
sx q[0];
rz(-0.98310982) q[0];
rz(-pi) q[1];
rz(0.7943031) q[2];
sx q[2];
rz(-1.0216733) q[2];
sx q[2];
rz(1.4505475) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8402108) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(2.0135897) q[1];
rz(-1.4521452) q[3];
sx q[3];
rz(-1.7735529) q[3];
sx q[3];
rz(-2.6076536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9459076) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(0.53608981) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.6528116) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(0.39500239) q[0];
rz(-1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(-2.1283456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7618338) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(0.14302111) q[0];
rz(-pi) q[1];
rz(-0.35372325) q[2];
sx q[2];
rz(-1.5295267) q[2];
sx q[2];
rz(-1.8281787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71868616) q[1];
sx q[1];
rz(-1.4740853) q[1];
sx q[1];
rz(2.7194517) q[1];
x q[2];
rz(2.2447484) q[3];
sx q[3];
rz(-2.2138811) q[3];
sx q[3];
rz(0.53900063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(1.7761207) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41291819) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(-1.7309932) q[2];
sx q[2];
rz(-2.6916531) q[2];
sx q[2];
rz(2.0315363) q[2];
rz(-1.6462973) q[3];
sx q[3];
rz(-0.69926881) q[3];
sx q[3];
rz(0.26509501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
