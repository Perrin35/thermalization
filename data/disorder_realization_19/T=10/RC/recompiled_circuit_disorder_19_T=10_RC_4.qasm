OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0269545) q[0];
sx q[0];
rz(4.5933525) q[0];
sx q[0];
rz(10.070355) q[0];
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(-1.6436613) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49255532) q[0];
sx q[0];
rz(-1.7151924) q[0];
sx q[0];
rz(-1.9306246) q[0];
rz(-pi) q[1];
rz(1.6662007) q[2];
sx q[2];
rz(-2.5089052) q[2];
sx q[2];
rz(0.22104095) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9420535) q[1];
sx q[1];
rz(-2.374211) q[1];
sx q[1];
rz(2.8295423) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50819355) q[3];
sx q[3];
rz(-0.71532202) q[3];
sx q[3];
rz(0.040166044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92007414) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-2.8011838) q[2];
rz(0.83299625) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.3566646) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67265636) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(0.54498589) q[0];
rz(-2.2333721) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(-2.3166336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.724052) q[0];
sx q[0];
rz(-1.1678956) q[0];
sx q[0];
rz(2.9379662) q[0];
x q[1];
rz(2.7386875) q[2];
sx q[2];
rz(-1.5803792) q[2];
sx q[2];
rz(2.8275037) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35926871) q[1];
sx q[1];
rz(-1.6869079) q[1];
sx q[1];
rz(-0.5639204) q[1];
x q[2];
rz(-1.2689781) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(-0.9580108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5543582) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-2.6718455) q[0];
rz(-1.5247955) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(-3.1030531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32917133) q[0];
sx q[0];
rz(-3.1338131) q[0];
sx q[0];
rz(1.6532941) q[0];
rz(-pi) q[1];
rz(1.4171717) q[2];
sx q[2];
rz(-1.3426174) q[2];
sx q[2];
rz(-1.4128026) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5765499) q[1];
sx q[1];
rz(-1.0817263) q[1];
sx q[1];
rz(-3.0719041) q[1];
rz(0.58979804) q[3];
sx q[3];
rz(-1.0538978) q[3];
sx q[3];
rz(0.043881744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(-2.2290686) q[2];
rz(1.8330666) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7581166) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(-0.34969774) q[0];
rz(1.2448467) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(0.19515881) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53699025) q[0];
sx q[0];
rz(-1.6911117) q[0];
sx q[0];
rz(-0.093219482) q[0];
rz(-pi) q[1];
rz(-0.85907016) q[2];
sx q[2];
rz(-2.4675183) q[2];
sx q[2];
rz(2.7328343) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8198204) q[1];
sx q[1];
rz(-2.0211377) q[1];
sx q[1];
rz(0.26178534) q[1];
rz(-2.5024662) q[3];
sx q[3];
rz(-1.5838606) q[3];
sx q[3];
rz(-0.55707219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.267743) q[2];
rz(1.0686482) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0347663) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(3.0019794) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-2.7635014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17079167) q[0];
sx q[0];
rz(-2.1149201) q[0];
sx q[0];
rz(0.32969726) q[0];
rz(-0.81850448) q[2];
sx q[2];
rz(-2.2193925) q[2];
sx q[2];
rz(-1.4865781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8564344) q[1];
sx q[1];
rz(-2.563623) q[1];
sx q[1];
rz(-2.2846089) q[1];
rz(-1.3729172) q[3];
sx q[3];
rz(-1.1117001) q[3];
sx q[3];
rz(1.0517694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2698764) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(-2.2536229) q[2];
rz(-2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(-2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75509214) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(2.8185484) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(-0.91517085) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.299303) q[0];
sx q[0];
rz(-1.248484) q[0];
sx q[0];
rz(-1.6567985) q[0];
x q[1];
rz(-0.50778163) q[2];
sx q[2];
rz(-1.0630597) q[2];
sx q[2];
rz(2.7024262) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9050161) q[1];
sx q[1];
rz(-1.502115) q[1];
sx q[1];
rz(0.31355365) q[1];
x q[2];
rz(-0.1359453) q[3];
sx q[3];
rz(-2.4424986) q[3];
sx q[3];
rz(-1.618315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.53720981) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(-2.3510695) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(-0.61736068) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73615605) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(-0.36703584) q[0];
rz(1.2332747) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6293684) q[0];
sx q[0];
rz(-2.6656796) q[0];
sx q[0];
rz(-0.41643629) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13397127) q[2];
sx q[2];
rz(-2.6920762) q[2];
sx q[2];
rz(1.9884895) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.047445) q[1];
sx q[1];
rz(-0.74456753) q[1];
sx q[1];
rz(-0.57569699) q[1];
rz(-pi) q[2];
rz(1.5464877) q[3];
sx q[3];
rz(-0.50384249) q[3];
sx q[3];
rz(0.77038308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1620862) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(-0.0018421729) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2974671) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(0.076518245) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(0.75235596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9936179) q[0];
sx q[0];
rz(-0.60950845) q[0];
sx q[0];
rz(0.57172914) q[0];
rz(-pi) q[1];
rz(-1.760375) q[2];
sx q[2];
rz(-0.65738064) q[2];
sx q[2];
rz(-1.1084523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68554316) q[1];
sx q[1];
rz(-1.6838264) q[1];
sx q[1];
rz(-3.0526524) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8389143) q[3];
sx q[3];
rz(-0.61121002) q[3];
sx q[3];
rz(-1.4509033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(3.1414462) q[2];
rz(-2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(1.0060271) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(1.3148274) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00025230322) q[0];
sx q[0];
rz(-2.25657) q[0];
sx q[0];
rz(1.89639) q[0];
rz(2.7776412) q[2];
sx q[2];
rz(-1.4462573) q[2];
sx q[2];
rz(-0.044791128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41373738) q[1];
sx q[1];
rz(-1.1579885) q[1];
sx q[1];
rz(-0.84819838) q[1];
rz(-pi) q[2];
rz(-1.5789088) q[3];
sx q[3];
rz(-0.74339657) q[3];
sx q[3];
rz(-0.5589561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43977794) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4196639) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(-0.28668177) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(-0.39189664) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.648801) q[0];
sx q[0];
rz(-2.1282882) q[0];
sx q[0];
rz(3.1317943) q[0];
x q[1];
rz(1.7461807) q[2];
sx q[2];
rz(-0.66329623) q[2];
sx q[2];
rz(1.3822008) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0185768) q[1];
sx q[1];
rz(-1.3850817) q[1];
sx q[1];
rz(-0.8155483) q[1];
x q[2];
rz(-0.55962555) q[3];
sx q[3];
rz(-1.2979753) q[3];
sx q[3];
rz(2.9469957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5320756) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.1575451) q[2];
rz(-2.5907497) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75523238) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(-1.339284) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(1.7427313) q[2];
sx q[2];
rz(-2.3430765) q[2];
sx q[2];
rz(0.55745468) q[2];
rz(-0.34655456) q[3];
sx q[3];
rz(-2.0934436) q[3];
sx q[3];
rz(2.2128076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];