OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(2.0164665) q[0];
rz(-5.2611051) q[1];
sx q[1];
rz(2.4740969) q[1];
sx q[1];
rz(11.752887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047377) q[0];
sx q[0];
rz(-1.9362402) q[0];
sx q[0];
rz(-2.9904537) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6935319) q[2];
sx q[2];
rz(-1.8984814) q[2];
sx q[2];
rz(2.0304012) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58373815) q[1];
sx q[1];
rz(-1.1314794) q[1];
sx q[1];
rz(-1.735461) q[1];
rz(1.068068) q[3];
sx q[3];
rz(-0.18618551) q[3];
sx q[3];
rz(-1.26621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4890613) q[2];
sx q[2];
rz(-1.4514613) q[2];
sx q[2];
rz(0.92864621) q[2];
rz(-1.5993902) q[3];
sx q[3];
rz(-1.8079115) q[3];
sx q[3];
rz(1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20392513) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(3.0145338) q[0];
rz(2.1584885) q[1];
sx q[1];
rz(-1.3652722) q[1];
sx q[1];
rz(-0.7712706) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5708904) q[0];
sx q[0];
rz(-1.9237907) q[0];
sx q[0];
rz(0.81309005) q[0];
rz(-pi) q[1];
rz(-0.060734435) q[2];
sx q[2];
rz(-1.2337451) q[2];
sx q[2];
rz(0.13523808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3202618) q[1];
sx q[1];
rz(-1.2148569) q[1];
sx q[1];
rz(-1.9411646) q[1];
rz(2.1089606) q[3];
sx q[3];
rz(-2.0797044) q[3];
sx q[3];
rz(-2.5050688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1473006) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(3.1331983) q[2];
rz(-0.66347915) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(-2.8765163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0405149) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(2.7000632) q[0];
rz(2.1532374) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(-3.006014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7429333) q[0];
sx q[0];
rz(-0.017983111) q[0];
sx q[0];
rz(-0.44264408) q[0];
rz(-2.709948) q[2];
sx q[2];
rz(-0.76485357) q[2];
sx q[2];
rz(-0.60827309) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24568711) q[1];
sx q[1];
rz(-1.0789597) q[1];
sx q[1];
rz(-0.40426429) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3531923) q[3];
sx q[3];
rz(-0.74624589) q[3];
sx q[3];
rz(-1.3775795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93379891) q[2];
sx q[2];
rz(-1.7094882) q[2];
sx q[2];
rz(0.03820339) q[2];
rz(0.52538747) q[3];
sx q[3];
rz(-0.62756687) q[3];
sx q[3];
rz(0.6853404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4811089) q[0];
sx q[0];
rz(-2.3479192) q[0];
sx q[0];
rz(-0.61087459) q[0];
rz(1.8065709) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(-0.23922051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51620519) q[0];
sx q[0];
rz(-0.62407485) q[0];
sx q[0];
rz(1.1136057) q[0];
rz(-0.68065721) q[2];
sx q[2];
rz(-0.58154642) q[2];
sx q[2];
rz(-0.30738632) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4799616) q[1];
sx q[1];
rz(-2.0212272) q[1];
sx q[1];
rz(-1.9305139) q[1];
x q[2];
rz(2.6801706) q[3];
sx q[3];
rz(-1.7506101) q[3];
sx q[3];
rz(-1.622004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4227582) q[2];
sx q[2];
rz(-1.6966635) q[2];
sx q[2];
rz(0.0017496721) q[2];
rz(2.8478029) q[3];
sx q[3];
rz(-1.9521451) q[3];
sx q[3];
rz(0.26688117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2413498) q[0];
sx q[0];
rz(-1.3768063) q[0];
sx q[0];
rz(2.2985261) q[0];
rz(-0.33755606) q[1];
sx q[1];
rz(-1.866021) q[1];
sx q[1];
rz(1.6036124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3030515) q[0];
sx q[0];
rz(-1.3796796) q[0];
sx q[0];
rz(-2.3421846) q[0];
rz(2.847371) q[2];
sx q[2];
rz(-0.18202848) q[2];
sx q[2];
rz(0.95532571) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5462436) q[1];
sx q[1];
rz(-1.3145295) q[1];
sx q[1];
rz(1.7728189) q[1];
rz(-pi) q[2];
rz(1.4902275) q[3];
sx q[3];
rz(-1.779883) q[3];
sx q[3];
rz(-2.5029456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.43188492) q[2];
sx q[2];
rz(-1.0872492) q[2];
sx q[2];
rz(-1.0464) q[2];
rz(2.8228068) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(-2.5939202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980314) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(2.4936254) q[0];
rz(-1.7421534) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(-1.8125777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15547046) q[0];
sx q[0];
rz(-0.24674812) q[0];
sx q[0];
rz(1.107649) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5115254) q[2];
sx q[2];
rz(-0.35379256) q[2];
sx q[2];
rz(1.1902863) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15550286) q[1];
sx q[1];
rz(-2.2515319) q[1];
sx q[1];
rz(0.40386856) q[1];
rz(-pi) q[2];
rz(2.5446114) q[3];
sx q[3];
rz(-0.59904811) q[3];
sx q[3];
rz(-2.0839276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38137388) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(1.5927429) q[2];
rz(0.069843944) q[3];
sx q[3];
rz(-1.2589688) q[3];
sx q[3];
rz(0.86863345) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0233362) q[0];
sx q[0];
rz(-1.9255487) q[0];
sx q[0];
rz(-0.94183952) q[0];
rz(-2.9815004) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(-0.19217415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9401682) q[0];
sx q[0];
rz(-1.7946825) q[0];
sx q[0];
rz(-0.036947771) q[0];
x q[1];
rz(3.1015119) q[2];
sx q[2];
rz(-2.0922086) q[2];
sx q[2];
rz(-1.5837976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77171626) q[1];
sx q[1];
rz(-1.2315005) q[1];
sx q[1];
rz(0.75668613) q[1];
rz(-pi) q[2];
rz(-3.1195379) q[3];
sx q[3];
rz(-2.0682749) q[3];
sx q[3];
rz(-0.22145222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75752246) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(2.7080217) q[2];
rz(-1.5844257) q[3];
sx q[3];
rz(-1.6470563) q[3];
sx q[3];
rz(2.7092194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(1.3442511) q[0];
rz(1.9380219) q[1];
sx q[1];
rz(-1.9529587) q[1];
sx q[1];
rz(-1.7787836) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2479073) q[0];
sx q[0];
rz(-3.1046668) q[0];
sx q[0];
rz(-0.93494995) q[0];
rz(2.9672181) q[2];
sx q[2];
rz(-0.83372901) q[2];
sx q[2];
rz(2.9437609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92719936) q[1];
sx q[1];
rz(-1.3424557) q[1];
sx q[1];
rz(-1.4424999) q[1];
rz(-pi) q[2];
rz(0.31764389) q[3];
sx q[3];
rz(-0.90932019) q[3];
sx q[3];
rz(2.3931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5726996) q[2];
sx q[2];
rz(-2.3685444) q[2];
sx q[2];
rz(2.1577238) q[2];
rz(0.24108663) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(-0.69055313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47138658) q[0];
sx q[0];
rz(-2.3455878) q[0];
sx q[0];
rz(-2.8669226) q[0];
rz(-1.7999016) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(-2.0955657) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90993308) q[0];
sx q[0];
rz(-2.446736) q[0];
sx q[0];
rz(0.52082638) q[0];
x q[1];
rz(3.0606424) q[2];
sx q[2];
rz(-1.3737203) q[2];
sx q[2];
rz(1.0287036) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6089692) q[1];
sx q[1];
rz(-1.6358422) q[1];
sx q[1];
rz(-2.5576127) q[1];
rz(1.0677797) q[3];
sx q[3];
rz(-1.5340337) q[3];
sx q[3];
rz(-0.063677841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86299738) q[2];
sx q[2];
rz(-1.3906761) q[2];
sx q[2];
rz(-2.7421303) q[2];
rz(0.56600371) q[3];
sx q[3];
rz(-0.96650201) q[3];
sx q[3];
rz(-2.32617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852916) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(2.5592819) q[0];
rz(-2.9583926) q[1];
sx q[1];
rz(-2.2661426) q[1];
sx q[1];
rz(2.3351672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9574653) q[0];
sx q[0];
rz(-0.13617198) q[0];
sx q[0];
rz(2.2720112) q[0];
rz(-0.33558947) q[2];
sx q[2];
rz(-1.3177455) q[2];
sx q[2];
rz(-2.0584681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88397056) q[1];
sx q[1];
rz(-0.7695573) q[1];
sx q[1];
rz(-2.8474924) q[1];
rz(-pi) q[2];
rz(-2.8424758) q[3];
sx q[3];
rz(-0.36206216) q[3];
sx q[3];
rz(2.5821834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23239423) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(2.3363028) q[2];
rz(2.9946839) q[3];
sx q[3];
rz(-0.78607905) q[3];
sx q[3];
rz(2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2522226) q[0];
sx q[0];
rz(-1.3852373) q[0];
sx q[0];
rz(-1.2607384) q[0];
rz(-1.5421142) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(0.11453828) q[2];
sx q[2];
rz(-2.2907612) q[2];
sx q[2];
rz(-1.9775122) q[2];
rz(-2.4965548) q[3];
sx q[3];
rz(-1.2923664) q[3];
sx q[3];
rz(-0.016135767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
