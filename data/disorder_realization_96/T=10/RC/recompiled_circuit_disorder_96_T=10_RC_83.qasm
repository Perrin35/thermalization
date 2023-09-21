OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(0.81737104) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(-1.2004381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5759597) q[0];
sx q[0];
rz(-2.6607249) q[0];
sx q[0];
rz(2.5984882) q[0];
rz(-pi) q[1];
rz(-0.13237662) q[2];
sx q[2];
rz(-2.1059603) q[2];
sx q[2];
rz(-1.3995427) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5719205) q[1];
sx q[1];
rz(-1.8351646) q[1];
sx q[1];
rz(-2.4155248) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3463763) q[3];
sx q[3];
rz(-1.3391558) q[3];
sx q[3];
rz(0.29602805) q[3];
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
rz(-1.583064) q[2];
rz(0.99672404) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834171) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(1.9460829) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(-0.53584677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5553404) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(0.14848407) q[0];
x q[1];
rz(0.29157721) q[2];
sx q[2];
rz(-0.73086408) q[2];
sx q[2];
rz(2.4183194) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1670926) q[1];
sx q[1];
rz(-1.316861) q[1];
sx q[1];
rz(0.42029917) q[1];
rz(-pi) q[2];
rz(0.032580839) q[3];
sx q[3];
rz(-1.4279281) q[3];
sx q[3];
rz(-2.8829765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0216996) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(-0.066453233) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(-2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41985837) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(2.9911175) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(3.1157852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8255071) q[0];
sx q[0];
rz(-1.320991) q[0];
sx q[0];
rz(-3.120963) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8823207) q[2];
sx q[2];
rz(-1.6561964) q[2];
sx q[2];
rz(2.4740919) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64297134) q[1];
sx q[1];
rz(-1.9830623) q[1];
sx q[1];
rz(2.4436414) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0098626) q[3];
sx q[3];
rz(-1.1909435) q[3];
sx q[3];
rz(-0.47117885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1228483) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(0.5853816) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(-0.88090849) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-2.6054629) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8857408) q[0];
sx q[0];
rz(-0.8937853) q[0];
sx q[0];
rz(-2.6141502) q[0];
rz(-1.5136112) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(1.1255217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0587479) q[1];
sx q[1];
rz(-1.9472329) q[1];
sx q[1];
rz(-2.46545) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4651277) q[3];
sx q[3];
rz(-1.3389265) q[3];
sx q[3];
rz(1.7959309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6716016) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(-1.0446576) q[2];
rz(-2.4345543) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(-2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2351284) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(-1.4936739) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(0.043118127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9263822) q[0];
sx q[0];
rz(-0.97391093) q[0];
sx q[0];
rz(2.9125288) q[0];
rz(0.28508614) q[2];
sx q[2];
rz(-1.0694155) q[2];
sx q[2];
rz(-1.3146871) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67152126) q[1];
sx q[1];
rz(-0.38362353) q[1];
sx q[1];
rz(0.77312153) q[1];
rz(-0.14857265) q[3];
sx q[3];
rz(-2.3464977) q[3];
sx q[3];
rz(3.1032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(-0.35219231) q[2];
rz(-2.5514065) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5181638) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(1.9859001) q[0];
rz(2.39134) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(-2.0828784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9583225) q[0];
sx q[0];
rz(-2.000862) q[0];
sx q[0];
rz(-1.3413315) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52144737) q[2];
sx q[2];
rz(-1.8910742) q[2];
sx q[2];
rz(-2.5495868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8844879) q[1];
sx q[1];
rz(-0.37933644) q[1];
sx q[1];
rz(0.92909716) q[1];
rz(-pi) q[2];
rz(-2.433957) q[3];
sx q[3];
rz(-1.4228627) q[3];
sx q[3];
rz(-1.1843475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6713312) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.2188101) q[2];
rz(-1.9865215) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.7003805) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041615151) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(-1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(-0.84164936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.323303) q[0];
sx q[0];
rz(-1.0536195) q[0];
sx q[0];
rz(-1.8316168) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62692554) q[2];
sx q[2];
rz(-1.4142087) q[2];
sx q[2];
rz(0.92948929) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0591653) q[1];
sx q[1];
rz(-1.6470243) q[1];
sx q[1];
rz(0.04870292) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3326725) q[3];
sx q[3];
rz(-2.3861285) q[3];
sx q[3];
rz(-1.8734224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1402309) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(3.0916396) q[2];
rz(-0.66155457) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(-2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89896232) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(1.7393973) q[0];
rz(-0.095480355) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4293489) q[0];
sx q[0];
rz(-2.0040383) q[0];
sx q[0];
rz(-0.54746763) q[0];
x q[1];
rz(-1.278644) q[2];
sx q[2];
rz(-2.7610965) q[2];
sx q[2];
rz(2.4351062) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3216746) q[1];
sx q[1];
rz(-0.770861) q[1];
sx q[1];
rz(-0.8701156) q[1];
x q[2];
rz(2.8475259) q[3];
sx q[3];
rz(-1.9293647) q[3];
sx q[3];
rz(2.8420574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(-2.0987089) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(-2.5119264) q[0];
rz(-0.57485238) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-2.1946857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3627975) q[0];
sx q[0];
rz(-0.4501833) q[0];
sx q[0];
rz(2.3953526) q[0];
rz(-pi) q[1];
rz(2.7526555) q[2];
sx q[2];
rz(-3.0367594) q[2];
sx q[2];
rz(-2.7742085) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0526035) q[1];
sx q[1];
rz(-0.34480428) q[1];
sx q[1];
rz(0.25119541) q[1];
x q[2];
rz(-0.19091786) q[3];
sx q[3];
rz(-1.1357765) q[3];
sx q[3];
rz(-2.2866979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56069121) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(-2.8760288) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749851) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(0.65504909) q[0];
rz(-2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(-0.64430976) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5376741) q[0];
sx q[0];
rz(-2.9593421) q[0];
sx q[0];
rz(-2.8331579) q[0];
rz(-pi) q[1];
rz(2.2868025) q[2];
sx q[2];
rz(-1.7456747) q[2];
sx q[2];
rz(1.0027494) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6380784) q[1];
sx q[1];
rz(-2.6942309) q[1];
sx q[1];
rz(0.7034941) q[1];
rz(-pi) q[2];
rz(-0.82018567) q[3];
sx q[3];
rz(-2.5327442) q[3];
sx q[3];
rz(-1.604515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3506938) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(2.541686) q[2];
rz(-0.89896262) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(-1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29466378) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(0.22944336) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(2.1951998) q[2];
sx q[2];
rz(-1.6740435) q[2];
sx q[2];
rz(-1.2798535) q[2];
rz(-0.39137822) q[3];
sx q[3];
rz(-2.2073675) q[3];
sx q[3];
rz(-1.7195306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];