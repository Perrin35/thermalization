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
rz(2.4861205) q[0];
sx q[0];
rz(-0.95215005) q[0];
sx q[0];
rz(3.0478391) q[0];
rz(-1.7733511) q[1];
sx q[1];
rz(-1.5612839) q[1];
sx q[1];
rz(-2.4743647) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.102948) q[0];
sx q[0];
rz(-1.1955452) q[0];
sx q[0];
rz(1.0374336) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66114963) q[2];
sx q[2];
rz(-0.3232269) q[2];
sx q[2];
rz(1.8913392) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2591483) q[1];
sx q[1];
rz(-1.2899152) q[1];
sx q[1];
rz(2.2876863) q[1];
rz(2.2556858) q[3];
sx q[3];
rz(-2.0768171) q[3];
sx q[3];
rz(-2.9850849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5253456) q[2];
sx q[2];
rz(-1.5364001) q[2];
sx q[2];
rz(-0.023539143) q[2];
rz(-0.25447887) q[3];
sx q[3];
rz(-1.1602217) q[3];
sx q[3];
rz(0.75828534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.94754058) q[0];
sx q[0];
rz(-1.7451311) q[0];
sx q[0];
rz(-0.61554712) q[0];
rz(-0.60000348) q[1];
sx q[1];
rz(-1.871385) q[1];
sx q[1];
rz(0.34872762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6541512) q[0];
sx q[0];
rz(-2.2117227) q[0];
sx q[0];
rz(1.6665672) q[0];
x q[1];
rz(0.28962692) q[2];
sx q[2];
rz(-1.2420601) q[2];
sx q[2];
rz(-1.5266974) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80417577) q[1];
sx q[1];
rz(-1.2722208) q[1];
sx q[1];
rz(-2.0470624) q[1];
rz(-1.0828937) q[3];
sx q[3];
rz(-0.55888218) q[3];
sx q[3];
rz(-0.36039613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49150026) q[2];
sx q[2];
rz(-1.8365752) q[2];
sx q[2];
rz(2.4864062) q[2];
rz(-0.16264597) q[3];
sx q[3];
rz(-0.9442257) q[3];
sx q[3];
rz(-0.20774016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96838897) q[0];
sx q[0];
rz(-2.0252616) q[0];
sx q[0];
rz(3.047347) q[0];
rz(-1.8415797) q[1];
sx q[1];
rz(-0.899122) q[1];
sx q[1];
rz(-1.947044) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31904063) q[0];
sx q[0];
rz(-1.8739432) q[0];
sx q[0];
rz(0.52355741) q[0];
rz(3.0897806) q[2];
sx q[2];
rz(-1.9345043) q[2];
sx q[2];
rz(1.6339982) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.647741) q[1];
sx q[1];
rz(-1.5268699) q[1];
sx q[1];
rz(-1.8274587) q[1];
rz(-0.76285513) q[3];
sx q[3];
rz(-1.8486946) q[3];
sx q[3];
rz(-0.40464532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3208348) q[2];
sx q[2];
rz(-1.6944378) q[2];
sx q[2];
rz(-0.40795946) q[2];
rz(-1.7631433) q[3];
sx q[3];
rz(-0.89865509) q[3];
sx q[3];
rz(-3.0234226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57678643) q[0];
sx q[0];
rz(-1.7498359) q[0];
sx q[0];
rz(0.80192649) q[0];
rz(0.39464125) q[1];
sx q[1];
rz(-1.0486187) q[1];
sx q[1];
rz(0.035042979) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5749767) q[0];
sx q[0];
rz(-1.0270976) q[0];
sx q[0];
rz(-2.7046142) q[0];
rz(2.6357627) q[2];
sx q[2];
rz(-0.37908812) q[2];
sx q[2];
rz(-1.0729147) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9525682) q[1];
sx q[1];
rz(-3.0810297) q[1];
sx q[1];
rz(-2.3409055) q[1];
x q[2];
rz(0.35538843) q[3];
sx q[3];
rz(-1.3878763) q[3];
sx q[3];
rz(-1.8754043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0798215) q[2];
sx q[2];
rz(-2.0786736) q[2];
sx q[2];
rz(1.7064095) q[2];
rz(1.7363413) q[3];
sx q[3];
rz(-1.0900213) q[3];
sx q[3];
rz(-1.1100356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.6768796) q[0];
sx q[0];
rz(-1.1486624) q[0];
sx q[0];
rz(-1.1418463) q[0];
rz(-2.6308718) q[1];
sx q[1];
rz(-1.9074214) q[1];
sx q[1];
rz(1.7150735) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1349484) q[0];
sx q[0];
rz(-1.0492968) q[0];
sx q[0];
rz(-2.5000076) q[0];
rz(0.22682206) q[2];
sx q[2];
rz(-1.2661627) q[2];
sx q[2];
rz(-0.90639988) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50496265) q[1];
sx q[1];
rz(-0.98082029) q[1];
sx q[1];
rz(3.0061199) q[1];
rz(-0.1673836) q[3];
sx q[3];
rz(-1.1684844) q[3];
sx q[3];
rz(1.2237807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4003754) q[2];
sx q[2];
rz(-1.3072689) q[2];
sx q[2];
rz(-2.884414) q[2];
rz(2.2687965) q[3];
sx q[3];
rz(-1.3215439) q[3];
sx q[3];
rz(-2.0040373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1207101) q[0];
sx q[0];
rz(-2.6873984) q[0];
sx q[0];
rz(1.0783476) q[0];
rz(2.2348166) q[1];
sx q[1];
rz(-0.6858784) q[1];
sx q[1];
rz(-0.041898601) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2611194) q[0];
sx q[0];
rz(-2.6676237) q[0];
sx q[0];
rz(1.649545) q[0];
rz(-pi) q[1];
rz(1.4153775) q[2];
sx q[2];
rz(-2.4684484) q[2];
sx q[2];
rz(2.239925) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4087569) q[1];
sx q[1];
rz(-0.46161595) q[1];
sx q[1];
rz(2.2240586) q[1];
x q[2];
rz(-2.6073805) q[3];
sx q[3];
rz(-0.90972661) q[3];
sx q[3];
rz(0.40701696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8580253) q[2];
sx q[2];
rz(-2.4391386) q[2];
sx q[2];
rz(-2.2502327) q[2];
rz(-0.71409613) q[3];
sx q[3];
rz(-0.73914206) q[3];
sx q[3];
rz(2.0785296) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5080268) q[0];
sx q[0];
rz(-1.897568) q[0];
sx q[0];
rz(-0.6749534) q[0];
rz(-1.5114816) q[1];
sx q[1];
rz(-0.62050301) q[1];
sx q[1];
rz(2.564548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2733611) q[0];
sx q[0];
rz(-1.8087862) q[0];
sx q[0];
rz(1.0704051) q[0];
rz(-1.4680732) q[2];
sx q[2];
rz(-1.1396164) q[2];
sx q[2];
rz(-1.9959677) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31957175) q[1];
sx q[1];
rz(-0.66950018) q[1];
sx q[1];
rz(-1.7029525) q[1];
rz(-pi) q[2];
rz(2.5222328) q[3];
sx q[3];
rz(-2.1652048) q[3];
sx q[3];
rz(2.847282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.453489) q[2];
sx q[2];
rz(-1.1223531) q[2];
sx q[2];
rz(-0.92448676) q[2];
rz(-1.2201355) q[3];
sx q[3];
rz(-0.28885463) q[3];
sx q[3];
rz(-3.0815304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.39153758) q[0];
sx q[0];
rz(-2.4704762) q[0];
sx q[0];
rz(1.3439939) q[0];
rz(-1.2229819) q[1];
sx q[1];
rz(-1.5593301) q[1];
sx q[1];
rz(-1.5917684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90541461) q[0];
sx q[0];
rz(-1.2981725) q[0];
sx q[0];
rz(0.18698606) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7666367) q[2];
sx q[2];
rz(-1.6757312) q[2];
sx q[2];
rz(-2.8255759) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2651974) q[1];
sx q[1];
rz(-1.205447) q[1];
sx q[1];
rz(2.4202034) q[1];
rz(1.6652558) q[3];
sx q[3];
rz(-2.3324663) q[3];
sx q[3];
rz(2.4868696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.082077114) q[2];
sx q[2];
rz(-1.0249219) q[2];
sx q[2];
rz(0.22641851) q[2];
rz(-3.1021127) q[3];
sx q[3];
rz(-1.4791146) q[3];
sx q[3];
rz(1.1994908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3969642) q[0];
sx q[0];
rz(-1.2122943) q[0];
sx q[0];
rz(-0.82897559) q[0];
rz(-0.12116155) q[1];
sx q[1];
rz(-2.1794901) q[1];
sx q[1];
rz(1.6896348) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73614299) q[0];
sx q[0];
rz(-2.2847957) q[0];
sx q[0];
rz(-1.3881486) q[0];
rz(-2.317886) q[2];
sx q[2];
rz(-1.888477) q[2];
sx q[2];
rz(1.4448656) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1157997) q[1];
sx q[1];
rz(-0.83459548) q[1];
sx q[1];
rz(0.92206149) q[1];
rz(-pi) q[2];
rz(0.39912721) q[3];
sx q[3];
rz(-1.0525956) q[3];
sx q[3];
rz(-0.46405989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58763233) q[2];
sx q[2];
rz(-2.3162737) q[2];
sx q[2];
rz(-2.2770605) q[2];
rz(-0.65166059) q[3];
sx q[3];
rz(-2.103002) q[3];
sx q[3];
rz(-0.017875044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5038576) q[0];
sx q[0];
rz(-0.88541579) q[0];
sx q[0];
rz(0.46022415) q[0];
rz(1.7962615) q[1];
sx q[1];
rz(-2.5049152) q[1];
sx q[1];
rz(-1.3705137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84708285) q[0];
sx q[0];
rz(-2.6381603) q[0];
sx q[0];
rz(-2.2398021) q[0];
rz(-pi) q[1];
rz(-2.5039191) q[2];
sx q[2];
rz(-1.8013926) q[2];
sx q[2];
rz(1.892145) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3440086) q[1];
sx q[1];
rz(-1.8086149) q[1];
sx q[1];
rz(1.306755) q[1];
x q[2];
rz(-1.3636144) q[3];
sx q[3];
rz(-1.6231114) q[3];
sx q[3];
rz(2.4567043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3687849) q[2];
sx q[2];
rz(-1.3489172) q[2];
sx q[2];
rz(-2.9985912) q[2];
rz(-2.9755106) q[3];
sx q[3];
rz(-2.4487285) q[3];
sx q[3];
rz(-2.3486229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76846692) q[0];
sx q[0];
rz(-1.65092) q[0];
sx q[0];
rz(-1.7478818) q[0];
rz(1.985818) q[1];
sx q[1];
rz(-2.0230237) q[1];
sx q[1];
rz(0.3442234) q[1];
rz(-2.406698) q[2];
sx q[2];
rz(-1.9857039) q[2];
sx q[2];
rz(-2.3891941) q[2];
rz(2.803346) q[3];
sx q[3];
rz(-2.1460642) q[3];
sx q[3];
rz(1.2301302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
