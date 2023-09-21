OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(4.5259024) q[0];
sx q[0];
rz(10.685267) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7605654) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(-2.7678124) q[0];
rz(0.14416868) q[2];
sx q[2];
rz(-1.8509794) q[2];
sx q[2];
rz(0.35137128) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94581374) q[1];
sx q[1];
rz(-2.0683156) q[1];
sx q[1];
rz(-1.0832018) q[1];
rz(0.63763036) q[3];
sx q[3];
rz(-2.5507567) q[3];
sx q[3];
rz(1.6319815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(-2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(-1.4555567) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992422) q[0];
sx q[0];
rz(-0.99827168) q[0];
sx q[0];
rz(0.19897977) q[0];
rz(-pi) q[1];
rz(-2.935264) q[2];
sx q[2];
rz(-2.7596139) q[2];
sx q[2];
rz(2.1260335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5160408) q[1];
sx q[1];
rz(-1.5572773) q[1];
sx q[1];
rz(2.0936172) q[1];
x q[2];
rz(1.5442113) q[3];
sx q[3];
rz(-1.1592602) q[3];
sx q[3];
rz(-2.6708024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(-2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-2.341111) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(-1.172539) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55484178) q[0];
sx q[0];
rz(-1.8790073) q[0];
sx q[0];
rz(-2.5462333) q[0];
rz(-pi) q[1];
rz(2.795479) q[2];
sx q[2];
rz(-0.78595224) q[2];
sx q[2];
rz(1.9217938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8923924) q[1];
sx q[1];
rz(-1.8141659) q[1];
sx q[1];
rz(1.0358441) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3451194) q[3];
sx q[3];
rz(-0.25697069) q[3];
sx q[3];
rz(-0.32303177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(2.2303936) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38055414) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(-0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(2.6180843) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9273705) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(-2.5582696) q[0];
rz(-0.24387118) q[2];
sx q[2];
rz(-0.3393617) q[2];
sx q[2];
rz(1.7983758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3829141) q[1];
sx q[1];
rz(-2.3944693) q[1];
sx q[1];
rz(1.1927356) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0047699) q[3];
sx q[3];
rz(-2.1577583) q[3];
sx q[3];
rz(0.50859355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52544242) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(-2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48150912) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(-2.1496444) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46734738) q[0];
sx q[0];
rz(-1.6773946) q[0];
sx q[0];
rz(-2.9664413) q[0];
rz(-pi) q[1];
rz(3.132658) q[2];
sx q[2];
rz(-1.2679456) q[2];
sx q[2];
rz(2.6645899) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3757513) q[1];
sx q[1];
rz(-1.6110794) q[1];
sx q[1];
rz(2.8326616) q[1];
rz(-3.0117412) q[3];
sx q[3];
rz(-0.81749812) q[3];
sx q[3];
rz(2.0714456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(0.27080718) q[2];
rz(-2.9233542) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(-2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.7927992) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(-1.4250925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0483635) q[0];
sx q[0];
rz(-1.3472124) q[0];
sx q[0];
rz(2.6327052) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0503057) q[2];
sx q[2];
rz(-1.8249776) q[2];
sx q[2];
rz(-2.9448178) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6188366) q[1];
sx q[1];
rz(-2.6869876) q[1];
sx q[1];
rz(1.7230117) q[1];
rz(-pi) q[2];
rz(2.607843) q[3];
sx q[3];
rz(-1.0373877) q[3];
sx q[3];
rz(-2.651754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1340593) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(-2.5777204) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(-0.82180506) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069698378) q[0];
sx q[0];
rz(-2.0928203) q[0];
sx q[0];
rz(2.4126023) q[0];
rz(-pi) q[1];
rz(-0.49634883) q[2];
sx q[2];
rz(-2.2349572) q[2];
sx q[2];
rz(1.6830483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8824749) q[1];
sx q[1];
rz(-1.3438517) q[1];
sx q[1];
rz(2.4005753) q[1];
rz(0.58737289) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(-1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(-1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41641763) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(0.82292557) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.8364505) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37874052) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(0.61600323) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4370286) q[2];
sx q[2];
rz(-1.8578055) q[2];
sx q[2];
rz(-1.9412083) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1888684) q[1];
sx q[1];
rz(-2.3844068) q[1];
sx q[1];
rz(0.56307478) q[1];
rz(0.023530258) q[3];
sx q[3];
rz(-1.2315893) q[3];
sx q[3];
rz(-2.6587405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1398853) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(1.6513599) q[2];
rz(-1.0772609) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-0.3219147) q[0];
rz(1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-2.4386491) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212595) q[0];
sx q[0];
rz(-0.54176211) q[0];
sx q[0];
rz(-0.74777491) q[0];
rz(-pi) q[1];
rz(2.0800955) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(-2.408037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1785537) q[1];
sx q[1];
rz(-2.5759765) q[1];
sx q[1];
rz(-1.8815243) q[1];
rz(0.33396696) q[3];
sx q[3];
rz(-0.98200646) q[3];
sx q[3];
rz(-0.86563084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90074173) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(-0.30570269) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(-1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(-1.9158069) q[0];
rz(-2.2380791) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(0.46863619) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7077431) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(1.7953403) q[0];
rz(-pi) q[1];
rz(-1.956316) q[2];
sx q[2];
rz(-1.8744933) q[2];
sx q[2];
rz(2.7591443) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84506449) q[1];
sx q[1];
rz(-1.6390641) q[1];
sx q[1];
rz(1.3045842) q[1];
rz(-pi) q[2];
rz(-1.1115132) q[3];
sx q[3];
rz(-1.1772403) q[3];
sx q[3];
rz(0.72898385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(-1.998385) q[2];
rz(0.11463595) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464012) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(0.62190965) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(-0.68998228) q[2];
sx q[2];
rz(-2.1524515) q[2];
sx q[2];
rz(-0.099302789) q[2];
rz(-0.078483742) q[3];
sx q[3];
rz(-0.92072903) q[3];
sx q[3];
rz(-2.846684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];