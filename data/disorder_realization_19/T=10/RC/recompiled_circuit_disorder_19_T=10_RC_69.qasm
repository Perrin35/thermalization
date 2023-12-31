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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49255532) q[0];
sx q[0];
rz(-1.7151924) q[0];
sx q[0];
rz(1.210968) q[0];
rz(-pi) q[1];
rz(0.94028084) q[2];
sx q[2];
rz(-1.6271546) q[2];
sx q[2];
rz(1.2727357) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59939304) q[1];
sx q[1];
rz(-1.7855872) q[1];
sx q[1];
rz(-2.3989124) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1708158) q[3];
sx q[3];
rz(-0.96066374) q[3];
sx q[3];
rz(0.67584544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-0.34040889) q[2];
rz(0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(2.5966068) q[0];
rz(-0.90822059) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(2.3166336) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93281125) q[0];
sx q[0];
rz(-0.44890807) q[0];
sx q[0];
rz(-1.1277898) q[0];
rz(-3.1171563) q[2];
sx q[2];
rz(-0.40301286) q[2];
sx q[2];
rz(-1.8624061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3927041) q[1];
sx q[1];
rz(-0.57447937) q[1];
sx q[1];
rz(-2.9267465) q[1];
x q[2];
rz(-2.1373848) q[3];
sx q[3];
rz(-1.7363747) q[3];
sx q[3];
rz(-0.86629888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58723441) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(-2.9651802) q[2];
rz(3.0107064) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(3.1089354) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762887) q[0];
sx q[0];
rz(-0.58350199) q[0];
sx q[0];
rz(-2.6718455) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(-3.1030531) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32917133) q[0];
sx q[0];
rz(-0.0077795452) q[0];
sx q[0];
rz(1.4882985) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5589587) q[2];
sx q[2];
rz(-0.27432549) q[2];
sx q[2];
rz(2.3290616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.027027834) q[1];
sx q[1];
rz(-1.5092883) q[1];
sx q[1];
rz(2.0608749) q[1];
rz(-pi) q[2];
rz(0.97088082) q[3];
sx q[3];
rz(-1.0661134) q[3];
sx q[3];
rz(1.2952627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(-2.2290686) q[2];
rz(-1.3085261) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(1.7002038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38347605) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(1.8967459) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-0.19515881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0965658) q[0];
sx q[0];
rz(-1.4782527) q[0];
sx q[0];
rz(1.4499614) q[0];
rz(-0.85907016) q[2];
sx q[2];
rz(-2.4675183) q[2];
sx q[2];
rz(-0.40875834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2745167) q[1];
sx q[1];
rz(-1.80596) q[1];
sx q[1];
rz(-1.1067252) q[1];
rz(-0.6391265) q[3];
sx q[3];
rz(-1.5838606) q[3];
sx q[3];
rz(0.55707219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9106456) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(-1.267743) q[2];
rz(2.0729444) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(-3.0019794) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(0.37809125) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9169086) q[0];
sx q[0];
rz(-1.2901257) q[0];
sx q[0];
rz(-2.1397489) q[0];
rz(2.3374704) q[2];
sx q[2];
rz(-2.1466549) q[2];
sx q[2];
rz(-2.5428307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34199076) q[1];
sx q[1];
rz(-1.2050036) q[1];
sx q[1];
rz(2.0288543) q[1];
rz(-2.7630745) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(0.88796973) q[2];
rz(-0.97638431) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(-2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75509214) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(-0.37762541) q[0];
rz(2.8185484) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(-2.2264218) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.299303) q[0];
sx q[0];
rz(-1.8931086) q[0];
sx q[0];
rz(1.6567985) q[0];
rz(-pi) q[1];
x q[1];
rz(2.633811) q[2];
sx q[2];
rz(-1.0630597) q[2];
sx q[2];
rz(-0.43916647) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8296216) q[1];
sx q[1];
rz(-1.8835856) q[1];
sx q[1];
rz(1.4986067) q[1];
rz(-0.69453199) q[3];
sx q[3];
rz(-1.6581222) q[3];
sx q[3];
rz(2.9897523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(-2.3510695) q[2];
rz(-2.8516155) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(0.61736068) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(0.36703584) q[0];
rz(-1.2332747) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(1.7162011) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7084301) q[0];
sx q[0];
rz(-1.3843952) q[0];
sx q[0];
rz(0.44048803) q[0];
rz(-3.0076214) q[2];
sx q[2];
rz(-2.6920762) q[2];
sx q[2];
rz(-1.9884895) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9219626) q[1];
sx q[1];
rz(-1.9486517) q[1];
sx q[1];
rz(2.4835543) q[1];
x q[2];
rz(2.0745139) q[3];
sx q[3];
rz(-1.5590612) q[3];
sx q[3];
rz(-2.3624682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(0.0018421729) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(-0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2974671) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(-0.076518245) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(2.3892367) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233326) q[0];
sx q[0];
rz(-1.8857297) q[0];
sx q[0];
rz(-0.53091913) q[0];
rz(-1.3812177) q[2];
sx q[2];
rz(-2.484212) q[2];
sx q[2];
rz(2.0331403) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89531089) q[1];
sx q[1];
rz(-1.4824251) q[1];
sx q[1];
rz(-1.6842711) q[1];
rz(1.7767056) q[3];
sx q[3];
rz(-2.1504953) q[3];
sx q[3];
rz(1.0866144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2032808) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(-1.1095307) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(-1.0060271) q[0];
rz(-0.26793119) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.3148274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7816759) q[0];
sx q[0];
rz(-1.3206375) q[0];
sx q[0];
rz(-0.71235384) q[0];
x q[1];
rz(-2.8034231) q[2];
sx q[2];
rz(-0.38376946) q[2];
sx q[2];
rz(1.9308117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4115574) q[1];
sx q[1];
rz(-2.3282603) q[1];
sx q[1];
rz(0.98585339) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3141765) q[3];
sx q[3];
rz(-1.5653059) q[3];
sx q[3];
rz(2.1237802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4196639) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(0.28668177) q[0];
rz(0.87896705) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(2.749696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.648801) q[0];
sx q[0];
rz(-2.1282882) q[0];
sx q[0];
rz(0.0097983629) q[0];
rz(-pi) q[1];
rz(1.395412) q[2];
sx q[2];
rz(-0.66329623) q[2];
sx q[2];
rz(-1.3822008) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3956986) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(1.3032773) q[1];
rz(1.2519022) q[3];
sx q[3];
rz(-2.1074169) q[3];
sx q[3];
rz(1.5434138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(-1.9840476) q[2];
rz(0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75523238) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-2.9677283) q[2];
sx q[2];
rz(-0.78730135) q[2];
sx q[2];
rz(0.80136328) q[2];
rz(-2.1203534) q[3];
sx q[3];
rz(-1.2720576) q[3];
sx q[3];
rz(-2.6779327) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
