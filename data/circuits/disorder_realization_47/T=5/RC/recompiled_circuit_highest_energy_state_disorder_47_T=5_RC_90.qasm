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
rz(-1.4909622) q[0];
sx q[0];
rz(-1.7234252) q[0];
sx q[0];
rz(-1.6562847) q[0];
rz(-2.0909042) q[1];
sx q[1];
rz(4.8839999) q[1];
sx q[1];
rz(13.304741) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4838388) q[0];
sx q[0];
rz(-1.8138348) q[0];
sx q[0];
rz(-1.7263822) q[0];
rz(-pi) q[1];
rz(-1.1777056) q[2];
sx q[2];
rz(-1.7121268) q[2];
sx q[2];
rz(1.3648975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25827894) q[1];
sx q[1];
rz(-1.4313233) q[1];
sx q[1];
rz(-0.44513227) q[1];
rz(-pi) q[2];
rz(-1.2655696) q[3];
sx q[3];
rz(-2.1681227) q[3];
sx q[3];
rz(1.0553774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38306132) q[2];
sx q[2];
rz(-2.8585275) q[2];
sx q[2];
rz(2.7281249) q[2];
rz(-2.5773898) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(-1.9570501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41459945) q[0];
sx q[0];
rz(-2.0818384) q[0];
sx q[0];
rz(-0.10398908) q[0];
rz(-0.14101401) q[1];
sx q[1];
rz(-0.68599373) q[1];
sx q[1];
rz(2.7242421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7433421) q[0];
sx q[0];
rz(-2.5210327) q[0];
sx q[0];
rz(-2.1055566) q[0];
x q[1];
rz(1.5386469) q[2];
sx q[2];
rz(-2.7776383) q[2];
sx q[2];
rz(-0.049843069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5896259) q[1];
sx q[1];
rz(-1.5291277) q[1];
sx q[1];
rz(1.2677626) q[1];
x q[2];
rz(2.6127698) q[3];
sx q[3];
rz(-2.072087) q[3];
sx q[3];
rz(0.59893417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61791164) q[2];
sx q[2];
rz(-2.1126426) q[2];
sx q[2];
rz(2.9376612) q[2];
rz(-1.6265053) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6025036) q[0];
sx q[0];
rz(-1.4668377) q[0];
sx q[0];
rz(2.7432826) q[0];
rz(2.0643945) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(-0.55346742) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578149) q[0];
sx q[0];
rz(-2.4071065) q[0];
sx q[0];
rz(-0.015588394) q[0];
x q[1];
rz(1.3559504) q[2];
sx q[2];
rz(-2.2707377) q[2];
sx q[2];
rz(0.46589771) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3543713) q[1];
sx q[1];
rz(-1.0007035) q[1];
sx q[1];
rz(1.1524423) q[1];
x q[2];
rz(1.3952399) q[3];
sx q[3];
rz(-1.7581098) q[3];
sx q[3];
rz(-2.8009149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.074097721) q[2];
sx q[2];
rz(-0.40253887) q[2];
sx q[2];
rz(-1.2158166) q[2];
rz(2.1408234) q[3];
sx q[3];
rz(-1.7583022) q[3];
sx q[3];
rz(-1.8892939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85553402) q[0];
sx q[0];
rz(-2.0834041) q[0];
sx q[0];
rz(1.704294) q[0];
rz(-1.8057436) q[1];
sx q[1];
rz(-2.5156486) q[1];
sx q[1];
rz(0.45509532) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.06741) q[0];
sx q[0];
rz(-0.95522308) q[0];
sx q[0];
rz(-1.4727946) q[0];
rz(-pi) q[1];
rz(1.9432151) q[2];
sx q[2];
rz(-1.7657868) q[2];
sx q[2];
rz(2.7253828) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14070732) q[1];
sx q[1];
rz(-1.8155087) q[1];
sx q[1];
rz(-0.34414704) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20858553) q[3];
sx q[3];
rz(-1.6634395) q[3];
sx q[3];
rz(1.1535597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8670292) q[2];
sx q[2];
rz(-2.6919591) q[2];
sx q[2];
rz(0.8320128) q[2];
rz(2.4882107) q[3];
sx q[3];
rz(-0.91975206) q[3];
sx q[3];
rz(0.67460361) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5117383) q[0];
sx q[0];
rz(-1.8269704) q[0];
sx q[0];
rz(-2.4526556) q[0];
rz(0.2116994) q[1];
sx q[1];
rz(-1.4392821) q[1];
sx q[1];
rz(-0.45417085) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70063299) q[0];
sx q[0];
rz(-1.129375) q[0];
sx q[0];
rz(-0.9671797) q[0];
x q[1];
rz(-2.7709097) q[2];
sx q[2];
rz(-0.96763583) q[2];
sx q[2];
rz(2.2070845) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59133021) q[1];
sx q[1];
rz(-2.2095517) q[1];
sx q[1];
rz(-0.62229054) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0096142) q[3];
sx q[3];
rz(-0.30721634) q[3];
sx q[3];
rz(-0.12061435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5357369) q[2];
sx q[2];
rz(-1.3593295) q[2];
sx q[2];
rz(-2.6306756) q[2];
rz(-1.2153252) q[3];
sx q[3];
rz(-1.5646224) q[3];
sx q[3];
rz(0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11584347) q[0];
sx q[0];
rz(-2.8307493) q[0];
sx q[0];
rz(2.6281443) q[0];
rz(-0.91649857) q[1];
sx q[1];
rz(-1.9648353) q[1];
sx q[1];
rz(-1.7880012) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62666639) q[0];
sx q[0];
rz(-1.1276704) q[0];
sx q[0];
rz(1.4791545) q[0];
x q[1];
rz(1.9802092) q[2];
sx q[2];
rz(-1.9099897) q[2];
sx q[2];
rz(-1.8095176) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60823764) q[1];
sx q[1];
rz(-1.5441193) q[1];
sx q[1];
rz(2.6538626) q[1];
rz(-pi) q[2];
rz(-1.4818125) q[3];
sx q[3];
rz(-1.6355733) q[3];
sx q[3];
rz(1.3353867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4765656) q[2];
sx q[2];
rz(-2.5883784) q[2];
sx q[2];
rz(0.033585699) q[2];
rz(-2.636886) q[3];
sx q[3];
rz(-1.7028156) q[3];
sx q[3];
rz(2.3869042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018709239) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(-1.1705742) q[0];
rz(-0.19078828) q[1];
sx q[1];
rz(-0.46563322) q[1];
sx q[1];
rz(-0.64291397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87827728) q[0];
sx q[0];
rz(-1.8366792) q[0];
sx q[0];
rz(-3.0356867) q[0];
rz(-pi) q[1];
rz(2.8823938) q[2];
sx q[2];
rz(-1.2507696) q[2];
sx q[2];
rz(2.7519403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9906147) q[1];
sx q[1];
rz(-2.7508368) q[1];
sx q[1];
rz(-0.67024605) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6899818) q[3];
sx q[3];
rz(-1.9728807) q[3];
sx q[3];
rz(-1.5086255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99870318) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(-1.9996803) q[2];
rz(3.111908) q[3];
sx q[3];
rz(-1.6073062) q[3];
sx q[3];
rz(0.77965411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2340045) q[0];
sx q[0];
rz(-1.9588082) q[0];
sx q[0];
rz(1.3080066) q[0];
rz(-1.1169149) q[1];
sx q[1];
rz(-1.4968548) q[1];
sx q[1];
rz(-0.21211472) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4329047) q[0];
sx q[0];
rz(-1.2554662) q[0];
sx q[0];
rz(-0.84724119) q[0];
rz(-pi) q[1];
rz(2.8372357) q[2];
sx q[2];
rz(-1.3469014) q[2];
sx q[2];
rz(2.9606113) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6740283) q[1];
sx q[1];
rz(-1.0564359) q[1];
sx q[1];
rz(0.5089203) q[1];
rz(-2.2961286) q[3];
sx q[3];
rz(-2.0975631) q[3];
sx q[3];
rz(-0.23108069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9461296) q[2];
sx q[2];
rz(-1.3529494) q[2];
sx q[2];
rz(-1.8606261) q[2];
rz(-1.403275) q[3];
sx q[3];
rz(-0.91126982) q[3];
sx q[3];
rz(-2.4173071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69287777) q[0];
sx q[0];
rz(-2.4249478) q[0];
sx q[0];
rz(0.88248673) q[0];
rz(0.29300434) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(1.1262456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4871656) q[0];
sx q[0];
rz(-1.2058655) q[0];
sx q[0];
rz(0.69661822) q[0];
rz(-0.59992591) q[2];
sx q[2];
rz(-1.5740105) q[2];
sx q[2];
rz(-2.7300138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3346904) q[1];
sx q[1];
rz(-0.77065361) q[1];
sx q[1];
rz(-1.4937964) q[1];
rz(-0.4057986) q[3];
sx q[3];
rz(-1.8631336) q[3];
sx q[3];
rz(-0.97989156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5249411) q[2];
sx q[2];
rz(-0.61351073) q[2];
sx q[2];
rz(2.9743527) q[2];
rz(-3.1104769) q[3];
sx q[3];
rz(-1.1789221) q[3];
sx q[3];
rz(-0.64605609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-0.50877082) q[0];
sx q[0];
rz(-1.8080067) q[0];
sx q[0];
rz(-0.81047812) q[0];
rz(-1.3088016) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(-3.0442309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33695147) q[0];
sx q[0];
rz(-1.4048249) q[0];
sx q[0];
rz(2.354855) q[0];
rz(-pi) q[1];
rz(3.112744) q[2];
sx q[2];
rz(-1.9518153) q[2];
sx q[2];
rz(0.25298564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2557981) q[1];
sx q[1];
rz(-2.5558639) q[1];
sx q[1];
rz(0.08759193) q[1];
rz(-pi) q[2];
rz(-3.024289) q[3];
sx q[3];
rz(-1.8364834) q[3];
sx q[3];
rz(-2.9151288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29791609) q[2];
sx q[2];
rz(-1.0138136) q[2];
sx q[2];
rz(1.8664912) q[2];
rz(-2.7013333) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(-0.27537235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1210099) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(-0.62200017) q[1];
sx q[1];
rz(-1.2983464) q[1];
sx q[1];
rz(1.3141528) q[1];
rz(-1.4595819) q[2];
sx q[2];
rz(-1.4930475) q[2];
sx q[2];
rz(0.14676506) q[2];
rz(1.3986599) q[3];
sx q[3];
rz(-0.93440104) q[3];
sx q[3];
rz(1.222119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
