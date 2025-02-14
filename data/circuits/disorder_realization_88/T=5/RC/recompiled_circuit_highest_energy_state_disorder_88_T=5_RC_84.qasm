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
rz(-2.5992844) q[0];
sx q[0];
rz(-3.0071654) q[0];
sx q[0];
rz(-2.0943213) q[0];
rz(2.8174227) q[1];
sx q[1];
rz(-2.9919762) q[1];
sx q[1];
rz(-0.9447929) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0271929) q[0];
sx q[0];
rz(-1.2978221) q[0];
sx q[0];
rz(-1.9045715) q[0];
rz(-pi) q[1];
rz(-1.8763674) q[2];
sx q[2];
rz(-2.4726119) q[2];
sx q[2];
rz(-1.389735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.289669) q[1];
sx q[1];
rz(-2.0229368) q[1];
sx q[1];
rz(2.585212) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0534001) q[3];
sx q[3];
rz(-1.6783829) q[3];
sx q[3];
rz(1.7501329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2904539) q[2];
sx q[2];
rz(-1.6893427) q[2];
sx q[2];
rz(1.5645082) q[2];
rz(-2.5751298) q[3];
sx q[3];
rz(-2.5201859) q[3];
sx q[3];
rz(1.8433146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31747776) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(-2.4216006) q[0];
rz(0.27436817) q[1];
sx q[1];
rz(-0.73148483) q[1];
sx q[1];
rz(-0.55363399) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8648656) q[0];
sx q[0];
rz(-1.6867406) q[0];
sx q[0];
rz(1.9254294) q[0];
x q[1];
rz(-1.881734) q[2];
sx q[2];
rz(-1.9568517) q[2];
sx q[2];
rz(0.60901035) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6362379) q[1];
sx q[1];
rz(-0.84213187) q[1];
sx q[1];
rz(1.0716754) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52218826) q[3];
sx q[3];
rz(-2.0117674) q[3];
sx q[3];
rz(0.14312927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1729892) q[2];
sx q[2];
rz(-1.9393238) q[2];
sx q[2];
rz(2.7776862) q[2];
rz(0.11161741) q[3];
sx q[3];
rz(-2.2026187) q[3];
sx q[3];
rz(2.3811471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.1641721) q[0];
sx q[0];
rz(-0.82472473) q[0];
sx q[0];
rz(0.0053996276) q[0];
rz(-2.3258356) q[1];
sx q[1];
rz(-1.1382256) q[1];
sx q[1];
rz(0.38591787) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9078131) q[0];
sx q[0];
rz(-0.44615567) q[0];
sx q[0];
rz(-2.4491007) q[0];
rz(-2.1109525) q[2];
sx q[2];
rz(-0.6331501) q[2];
sx q[2];
rz(-1.4667222) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94066835) q[1];
sx q[1];
rz(-1.6977807) q[1];
sx q[1];
rz(-1.9534355) q[1];
x q[2];
rz(-1.8688275) q[3];
sx q[3];
rz(-2.7103817) q[3];
sx q[3];
rz(-1.9995156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58891121) q[2];
sx q[2];
rz(-0.83325714) q[2];
sx q[2];
rz(1.7595278) q[2];
rz(2.9803993) q[3];
sx q[3];
rz(-1.4313982) q[3];
sx q[3];
rz(2.5126357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2564119) q[0];
sx q[0];
rz(-1.8653402) q[0];
sx q[0];
rz(1.8026344) q[0];
rz(-1.8244052) q[1];
sx q[1];
rz(-2.1664797) q[1];
sx q[1];
rz(0.29979527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6000344) q[0];
sx q[0];
rz(-0.93342121) q[0];
sx q[0];
rz(2.0079028) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6671403) q[2];
sx q[2];
rz(-1.0336733) q[2];
sx q[2];
rz(-0.65309292) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.055017415) q[1];
sx q[1];
rz(-2.7624532) q[1];
sx q[1];
rz(-1.192548) q[1];
x q[2];
rz(-0.13238975) q[3];
sx q[3];
rz(-0.93919884) q[3];
sx q[3];
rz(-1.6211444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51912159) q[2];
sx q[2];
rz(-0.93834472) q[2];
sx q[2];
rz(-2.8224714) q[2];
rz(3.0506813) q[3];
sx q[3];
rz(-0.14486434) q[3];
sx q[3];
rz(1.3566141) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4153445) q[0];
sx q[0];
rz(-0.81979668) q[0];
sx q[0];
rz(-0.0023284624) q[0];
rz(1.5706515) q[1];
sx q[1];
rz(-2.3410773) q[1];
sx q[1];
rz(-1.194582) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1462458) q[0];
sx q[0];
rz(-1.4527391) q[0];
sx q[0];
rz(-2.9857568) q[0];
x q[1];
rz(-0.04045893) q[2];
sx q[2];
rz(-2.2881857) q[2];
sx q[2];
rz(2.3641912) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0280605) q[1];
sx q[1];
rz(-0.36021566) q[1];
sx q[1];
rz(1.2768145) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55641716) q[3];
sx q[3];
rz(-1.8117944) q[3];
sx q[3];
rz(-0.3410546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2583367) q[2];
sx q[2];
rz(-2.363435) q[2];
sx q[2];
rz(-3.0641595) q[2];
rz(0.94610131) q[3];
sx q[3];
rz(-2.6587722) q[3];
sx q[3];
rz(0.4536804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2133863) q[0];
sx q[0];
rz(-3.0551857) q[0];
sx q[0];
rz(-1.7267831) q[0];
rz(2.4746223) q[1];
sx q[1];
rz(-1.357115) q[1];
sx q[1];
rz(0.40474969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20807438) q[0];
sx q[0];
rz(-0.41941038) q[0];
sx q[0];
rz(-2.389877) q[0];
x q[1];
rz(0.91333977) q[2];
sx q[2];
rz(-1.7796675) q[2];
sx q[2];
rz(-3.0830736) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8903146) q[1];
sx q[1];
rz(-1.108544) q[1];
sx q[1];
rz(-3.0733215) q[1];
rz(-pi) q[2];
rz(-0.94958206) q[3];
sx q[3];
rz(-2.5646665) q[3];
sx q[3];
rz(-2.2872137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3177967) q[2];
sx q[2];
rz(-0.78936374) q[2];
sx q[2];
rz(2.6981567) q[2];
rz(-2.7277842) q[3];
sx q[3];
rz(-0.2897073) q[3];
sx q[3];
rz(-0.88576353) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39171788) q[0];
sx q[0];
rz(-1.3626008) q[0];
sx q[0];
rz(-2.474127) q[0];
rz(3.093847) q[1];
sx q[1];
rz(-2.0011438) q[1];
sx q[1];
rz(-1.1183636) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3781609) q[0];
sx q[0];
rz(-1.973828) q[0];
sx q[0];
rz(-0.58149882) q[0];
rz(-0.12130731) q[2];
sx q[2];
rz(-2.1459142) q[2];
sx q[2];
rz(2.2030769) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0065919) q[1];
sx q[1];
rz(-1.1524832) q[1];
sx q[1];
rz(-2.5043284) q[1];
x q[2];
rz(-1.4468071) q[3];
sx q[3];
rz(-2.638804) q[3];
sx q[3];
rz(0.10956746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.095801) q[2];
sx q[2];
rz(-1.2597151) q[2];
sx q[2];
rz(-3.1328787) q[2];
rz(-1.9376612) q[3];
sx q[3];
rz(-2.7976076) q[3];
sx q[3];
rz(3.1194527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0225723) q[0];
sx q[0];
rz(-2.75596) q[0];
sx q[0];
rz(-2.2517396) q[0];
rz(1.8564557) q[1];
sx q[1];
rz(-1.0533918) q[1];
sx q[1];
rz(3.0056675) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1401299) q[0];
sx q[0];
rz(-2.9076932) q[0];
sx q[0];
rz(0.20568307) q[0];
rz(1.0338342) q[2];
sx q[2];
rz(-1.0091678) q[2];
sx q[2];
rz(1.461535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.071819788) q[1];
sx q[1];
rz(-1.7679224) q[1];
sx q[1];
rz(0.33918753) q[1];
rz(-1.9210757) q[3];
sx q[3];
rz(-1.7052584) q[3];
sx q[3];
rz(2.2148481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4385684) q[2];
sx q[2];
rz(-0.88006222) q[2];
sx q[2];
rz(1.0620037) q[2];
rz(-2.2117173) q[3];
sx q[3];
rz(-2.3558741) q[3];
sx q[3];
rz(0.78833956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7787665) q[0];
sx q[0];
rz(-2.7500948) q[0];
sx q[0];
rz(2.7407001) q[0];
rz(2.3800384) q[1];
sx q[1];
rz(-1.4925894) q[1];
sx q[1];
rz(-2.3525995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3112199) q[0];
sx q[0];
rz(-1.5548176) q[0];
sx q[0];
rz(-0.1275087) q[0];
rz(-1.7640186) q[2];
sx q[2];
rz(-0.68422645) q[2];
sx q[2];
rz(2.2598337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4915062) q[1];
sx q[1];
rz(-1.4515001) q[1];
sx q[1];
rz(-0.59081282) q[1];
rz(-pi) q[2];
x q[2];
rz(0.085461334) q[3];
sx q[3];
rz(-2.6151867) q[3];
sx q[3];
rz(1.5861024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57364982) q[2];
sx q[2];
rz(-1.2929792) q[2];
sx q[2];
rz(2.4207777) q[2];
rz(1.4912262) q[3];
sx q[3];
rz(-2.3590915) q[3];
sx q[3];
rz(-1.8318374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1184621) q[0];
sx q[0];
rz(-3.0266302) q[0];
sx q[0];
rz(0.7777099) q[0];
rz(0.65981162) q[1];
sx q[1];
rz(-0.86645627) q[1];
sx q[1];
rz(2.7515817) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4068425) q[0];
sx q[0];
rz(-1.5336541) q[0];
sx q[0];
rz(1.7057306) q[0];
x q[1];
rz(-0.71163746) q[2];
sx q[2];
rz(-0.88274985) q[2];
sx q[2];
rz(-1.2104863) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46040487) q[1];
sx q[1];
rz(-1.8680267) q[1];
sx q[1];
rz(1.1497208) q[1];
x q[2];
rz(-0.50647363) q[3];
sx q[3];
rz(-0.9829384) q[3];
sx q[3];
rz(-1.7755058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4593792) q[2];
sx q[2];
rz(-0.43872681) q[2];
sx q[2];
rz(-2.6518346) q[2];
rz(-0.30719906) q[3];
sx q[3];
rz(-2.2178853) q[3];
sx q[3];
rz(-1.7507621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78140344) q[0];
sx q[0];
rz(-1.6506945) q[0];
sx q[0];
rz(-1.6682464) q[0];
rz(1.0745984) q[1];
sx q[1];
rz(-2.2886724) q[1];
sx q[1];
rz(1.2235175) q[1];
rz(0.75847738) q[2];
sx q[2];
rz(-0.52634326) q[2];
sx q[2];
rz(1.1545622) q[2];
rz(-0.60081595) q[3];
sx q[3];
rz(-2.0471341) q[3];
sx q[3];
rz(2.9772191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
