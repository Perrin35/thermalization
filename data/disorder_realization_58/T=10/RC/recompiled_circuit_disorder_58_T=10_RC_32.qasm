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
rz(2.7352754) q[0];
sx q[0];
rz(8.6046594) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(0.83067218) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2031189) q[0];
sx q[0];
rz(-1.9460558) q[0];
sx q[0];
rz(-1.9012326) q[0];
rz(0.73859282) q[2];
sx q[2];
rz(-1.2054218) q[2];
sx q[2];
rz(-1.5473168) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1066061) q[1];
sx q[1];
rz(-0.88284661) q[1];
sx q[1];
rz(0.41253849) q[1];
rz(-pi) q[2];
rz(2.6942263) q[3];
sx q[3];
rz(-1.2036185) q[3];
sx q[3];
rz(0.51608738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(0.1201771) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(-0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08081089) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-2.2297915) q[0];
rz(-0.78951019) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(0.3266913) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806843) q[0];
sx q[0];
rz(-0.93870367) q[0];
sx q[0];
rz(2.500781) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1268483) q[2];
sx q[2];
rz(-1.4269281) q[2];
sx q[2];
rz(2.4614046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3152299) q[1];
sx q[1];
rz(-2.3781207) q[1];
sx q[1];
rz(-0.96674322) q[1];
rz(-0.88300206) q[3];
sx q[3];
rz(-1.5044754) q[3];
sx q[3];
rz(1.4955213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6313173) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(-3.0316947) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-1.057391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43305106) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(-3.1197381) q[0];
x q[1];
rz(0.51809394) q[2];
sx q[2];
rz(-0.38803852) q[2];
sx q[2];
rz(2.5990017) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8101465) q[1];
sx q[1];
rz(-0.76008893) q[1];
sx q[1];
rz(2.0600832) q[1];
rz(-pi) q[2];
rz(-2.7189414) q[3];
sx q[3];
rz(-1.395441) q[3];
sx q[3];
rz(-0.56817504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(-2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(0.5350565) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-2.9761956) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1918614) q[0];
sx q[0];
rz(-0.22338578) q[0];
sx q[0];
rz(-0.97747691) q[0];
x q[1];
rz(-0.16343127) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(-0.93726678) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0488102) q[1];
sx q[1];
rz(-2.4895992) q[1];
sx q[1];
rz(-0.55744967) q[1];
rz(-pi) q[2];
rz(-0.27407077) q[3];
sx q[3];
rz(-2.1087397) q[3];
sx q[3];
rz(2.2171973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7455204) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(0.20544927) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(1.8849461) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(2.9096471) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1922798) q[0];
sx q[0];
rz(-0.67040196) q[0];
sx q[0];
rz(0.85286661) q[0];
rz(0.16817981) q[2];
sx q[2];
rz(-0.7025223) q[2];
sx q[2];
rz(2.0481734) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0143118) q[1];
sx q[1];
rz(-0.31186549) q[1];
sx q[1];
rz(-0.77906268) q[1];
x q[2];
rz(0.050458126) q[3];
sx q[3];
rz(-0.36894635) q[3];
sx q[3];
rz(1.3621393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.3034472) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(2.1525106) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4076685) q[0];
sx q[0];
rz(-1.7051538) q[0];
sx q[0];
rz(-0.863048) q[0];
x q[1];
rz(1.9428271) q[2];
sx q[2];
rz(-2.2673006) q[2];
sx q[2];
rz(-2.5859985) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4221229) q[1];
sx q[1];
rz(-0.49233961) q[1];
sx q[1];
rz(1.0455529) q[1];
rz(-pi) q[2];
rz(-1.0783844) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(1.2315962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9138907) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(2.690199) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-0.62414449) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(-2.5278032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95425883) q[0];
sx q[0];
rz(-1.7566534) q[0];
sx q[0];
rz(2.3128187) q[0];
x q[1];
rz(-0.28366144) q[2];
sx q[2];
rz(-1.5093056) q[2];
sx q[2];
rz(-0.21360699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63749667) q[1];
sx q[1];
rz(-1.6747961) q[1];
sx q[1];
rz(-1.5118447) q[1];
x q[2];
rz(1.3404487) q[3];
sx q[3];
rz(-2.6655572) q[3];
sx q[3];
rz(1.7649094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.9667352) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(-1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.90010086) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(-0.98446313) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1293837) q[0];
sx q[0];
rz(-2.1166271) q[0];
sx q[0];
rz(1.6145541) q[0];
rz(-0.41495277) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(-2.6331537) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7998432) q[1];
sx q[1];
rz(-1.5647596) q[1];
sx q[1];
rz(0.36532613) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47655388) q[3];
sx q[3];
rz(-1.2911951) q[3];
sx q[3];
rz(0.24147803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4349334) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(-0.6859268) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(-2.2156782) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1122702) q[0];
sx q[0];
rz(-1.7824714) q[0];
sx q[0];
rz(2.9979343) q[0];
x q[1];
rz(0.8530059) q[2];
sx q[2];
rz(-0.91663137) q[2];
sx q[2];
rz(-0.6086364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.029286413) q[1];
sx q[1];
rz(-2.6767113) q[1];
sx q[1];
rz(1.2390562) q[1];
x q[2];
rz(2.6191606) q[3];
sx q[3];
rz(-0.2345095) q[3];
sx q[3];
rz(2.0731376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9459076) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(2.6055028) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(2.1283456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2241867) q[0];
sx q[0];
rz(-1.7099483) q[0];
sx q[0];
rz(1.805472) q[0];
rz(-2.7878694) q[2];
sx q[2];
rz(-1.5295267) q[2];
sx q[2];
rz(1.8281787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4229065) q[1];
sx q[1];
rz(-1.4740853) q[1];
sx q[1];
rz(2.7194517) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3771044) q[3];
sx q[3];
rz(-1.0478684) q[3];
sx q[3];
rz(-0.58512277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(-2.3821793) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(2.571648) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41291819) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(-0.57327523) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(-1.1258833) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(-1.4952954) q[3];
sx q[3];
rz(-2.4423238) q[3];
sx q[3];
rz(-2.8764976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
