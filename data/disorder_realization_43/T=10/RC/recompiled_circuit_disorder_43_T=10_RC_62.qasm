OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2551978) q[0];
sx q[0];
rz(-1.891341) q[0];
sx q[0];
rz(1.3347081) q[0];
rz(2.788738) q[1];
sx q[1];
rz(-2.9810413) q[1];
sx q[1];
rz(-0.97595739) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87603509) q[0];
sx q[0];
rz(-1.8119438) q[0];
sx q[0];
rz(-0.59224706) q[0];
rz(1.0078148) q[2];
sx q[2];
rz(-1.8554167) q[2];
sx q[2];
rz(2.9993338) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0426892) q[1];
sx q[1];
rz(-0.96682036) q[1];
sx q[1];
rz(2.9788115) q[1];
rz(-pi) q[2];
rz(-1.1638155) q[3];
sx q[3];
rz(-2.3325936) q[3];
sx q[3];
rz(-2.80796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(-1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587104) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(0.006342412) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7550678) q[0];
sx q[0];
rz(-1.3773943) q[0];
sx q[0];
rz(0.13271876) q[0];
rz(2.6132934) q[2];
sx q[2];
rz(-1.0353147) q[2];
sx q[2];
rz(0.9165802) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7087242) q[1];
sx q[1];
rz(-1.8808865) q[1];
sx q[1];
rz(-1.6354145) q[1];
rz(-pi) q[2];
rz(1.8216324) q[3];
sx q[3];
rz(-2.1187966) q[3];
sx q[3];
rz(2.8743924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94770849) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-0.87810278) q[2];
rz(-2.2795423) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(-3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(-3.1306144) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(-0.095741622) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84155267) q[0];
sx q[0];
rz(-1.7471572) q[0];
sx q[0];
rz(-2.9862613) q[0];
rz(-pi) q[1];
rz(2.6355866) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(0.72088748) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1167736) q[1];
sx q[1];
rz(-2.1110592) q[1];
sx q[1];
rz(-2.5953672) q[1];
x q[2];
rz(-2.5903969) q[3];
sx q[3];
rz(-1.6008198) q[3];
sx q[3];
rz(-0.72202819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0720955) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(-2.3068008) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(-0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9469706) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(-2.6158781) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(-2.1077164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945275) q[0];
sx q[0];
rz(-2.3405502) q[0];
sx q[0];
rz(2.4497776) q[0];
rz(-2.6150377) q[2];
sx q[2];
rz(-1.0609846) q[2];
sx q[2];
rz(1.9499792) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8139207) q[1];
sx q[1];
rz(-1.3797626) q[1];
sx q[1];
rz(0.066285985) q[1];
x q[2];
rz(1.8978118) q[3];
sx q[3];
rz(-1.838284) q[3];
sx q[3];
rz(0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(0.41473266) q[0];
rz(1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(-2.5674852) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5363279) q[0];
sx q[0];
rz(-0.61075532) q[0];
sx q[0];
rz(-1.7954134) q[0];
rz(-pi) q[1];
rz(0.026002361) q[2];
sx q[2];
rz(-0.71547316) q[2];
sx q[2];
rz(2.0040087) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9628323) q[1];
sx q[1];
rz(-1.292359) q[1];
sx q[1];
rz(0.36414418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3247213) q[3];
sx q[3];
rz(-2.4613791) q[3];
sx q[3];
rz(-2.084793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.85270143) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.627702) q[2];
rz(-0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(0.16695887) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(0.73227698) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9376611) q[0];
sx q[0];
rz(-1.4121778) q[0];
sx q[0];
rz(-0.21477867) q[0];
x q[1];
rz(-2.2339091) q[2];
sx q[2];
rz(-2.216202) q[2];
sx q[2];
rz(1.4484608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69387586) q[1];
sx q[1];
rz(-0.69732053) q[1];
sx q[1];
rz(-2.7867774) q[1];
x q[2];
rz(-0.332571) q[3];
sx q[3];
rz(-0.59623527) q[3];
sx q[3];
rz(2.3308844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(1.9801271) q[2];
rz(-2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(2.9851595) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(2.3715473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8806481) q[0];
sx q[0];
rz(-2.5136607) q[0];
sx q[0];
rz(0.29675608) q[0];
rz(-pi) q[1];
rz(-0.40600834) q[2];
sx q[2];
rz(-1.0205262) q[2];
sx q[2];
rz(-1.5232435) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5321977) q[1];
sx q[1];
rz(-2.1708793) q[1];
sx q[1];
rz(-0.9432015) q[1];
rz(-pi) q[2];
rz(-0.50290147) q[3];
sx q[3];
rz(-1.8637878) q[3];
sx q[3];
rz(-0.46616947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(2.4027951) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-2.7752005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05224932) q[0];
sx q[0];
rz(-1.0567259) q[0];
sx q[0];
rz(0.71787562) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62944062) q[2];
sx q[2];
rz(-1.2120314) q[2];
sx q[2];
rz(-1.4251054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7494292) q[1];
sx q[1];
rz(-2.2656419) q[1];
sx q[1];
rz(1.3543345) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6493158) q[3];
sx q[3];
rz(-2.588387) q[3];
sx q[3];
rz(-1.5403403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1239132) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(-0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(-2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(2.4813095) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(-2.3502137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7991379) q[0];
sx q[0];
rz(-0.90478169) q[0];
sx q[0];
rz(-0.40823437) q[0];
rz(0.80957885) q[2];
sx q[2];
rz(-2.5123345) q[2];
sx q[2];
rz(1.6639683) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3190805) q[1];
sx q[1];
rz(-1.1573536) q[1];
sx q[1];
rz(2.6507069) q[1];
x q[2];
rz(-1.295624) q[3];
sx q[3];
rz(-0.55579805) q[3];
sx q[3];
rz(2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(1.0212612) q[2];
rz(-2.7630473) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(2.5436201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1148949) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(-2.5337906) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(1.7609319) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62771195) q[0];
sx q[0];
rz(-0.6427592) q[0];
sx q[0];
rz(-1.7895262) q[0];
x q[1];
rz(0.025823921) q[2];
sx q[2];
rz(-2.0061473) q[2];
sx q[2];
rz(3.0961852) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7541411) q[1];
sx q[1];
rz(-1.2536465) q[1];
sx q[1];
rz(2.0003358) q[1];
x q[2];
rz(1.7396183) q[3];
sx q[3];
rz(-1.3843378) q[3];
sx q[3];
rz(0.25690119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7878788) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(0.33622462) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040141) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(-2.5972988) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(-0.82472807) q[2];
sx q[2];
rz(-1.8246973) q[2];
sx q[2];
rz(1.4551103) q[2];
rz(2.1401134) q[3];
sx q[3];
rz(-1.7703198) q[3];
sx q[3];
rz(1.9386335) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];