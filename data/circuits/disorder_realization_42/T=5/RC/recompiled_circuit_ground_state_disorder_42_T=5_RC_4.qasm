OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.2494994) q[0];
sx q[0];
rz(4.3262859) q[0];
sx q[0];
rz(10.560137) q[0];
rz(-0.42179498) q[1];
sx q[1];
rz(-2.3611958) q[1];
sx q[1];
rz(2.5193522) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93995433) q[0];
sx q[0];
rz(-1.7169984) q[0];
sx q[0];
rz(-0.12473434) q[0];
rz(-pi) q[1];
rz(-0.038921629) q[2];
sx q[2];
rz(-1.1238255) q[2];
sx q[2];
rz(-2.1291422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2127004) q[1];
sx q[1];
rz(-1.5824741) q[1];
sx q[1];
rz(2.547193) q[1];
rz(-pi) q[2];
rz(-0.99986003) q[3];
sx q[3];
rz(-1.3821162) q[3];
sx q[3];
rz(-0.41485559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3368095) q[2];
sx q[2];
rz(-1.4562891) q[2];
sx q[2];
rz(3.0554904) q[2];
rz(-2.6297249) q[3];
sx q[3];
rz(-2.0972926) q[3];
sx q[3];
rz(0.69387236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7030199) q[0];
sx q[0];
rz(-0.95412552) q[0];
sx q[0];
rz(0.91617209) q[0];
rz(2.84962) q[1];
sx q[1];
rz(-2.5501854) q[1];
sx q[1];
rz(2.3672262) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2949897) q[0];
sx q[0];
rz(-1.2520391) q[0];
sx q[0];
rz(1.2002719) q[0];
rz(-pi) q[1];
rz(1.9161445) q[2];
sx q[2];
rz(-1.5651476) q[2];
sx q[2];
rz(-0.7330837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7978002) q[1];
sx q[1];
rz(-1.4693854) q[1];
sx q[1];
rz(0.12354688) q[1];
rz(2.1470469) q[3];
sx q[3];
rz(-2.3062097) q[3];
sx q[3];
rz(1.6540838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3463717) q[2];
sx q[2];
rz(-0.56196153) q[2];
sx q[2];
rz(0.58491659) q[2];
rz(1.6631205) q[3];
sx q[3];
rz(-1.9641179) q[3];
sx q[3];
rz(-1.6519206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65603489) q[0];
sx q[0];
rz(-1.0233044) q[0];
sx q[0];
rz(-2.5265332) q[0];
rz(-0.3282322) q[1];
sx q[1];
rz(-1.0098207) q[1];
sx q[1];
rz(-2.097791) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29662499) q[0];
sx q[0];
rz(-0.052694885) q[0];
sx q[0];
rz(-0.096676306) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1310668) q[2];
sx q[2];
rz(-0.42834706) q[2];
sx q[2];
rz(2.7017252) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9658326) q[1];
sx q[1];
rz(-1.9632259) q[1];
sx q[1];
rz(-0.85939851) q[1];
x q[2];
rz(-2.7451594) q[3];
sx q[3];
rz(-1.2370438) q[3];
sx q[3];
rz(-2.9891356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7528167) q[2];
sx q[2];
rz(-2.7478168) q[2];
sx q[2];
rz(0.30889312) q[2];
rz(2.654352) q[3];
sx q[3];
rz(-1.4332708) q[3];
sx q[3];
rz(-2.2545831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76172817) q[0];
sx q[0];
rz(-0.74493113) q[0];
sx q[0];
rz(0.47759011) q[0];
rz(-1.2848805) q[1];
sx q[1];
rz(-1.4297337) q[1];
sx q[1];
rz(-2.6755948) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1168524) q[0];
sx q[0];
rz(-2.5821487) q[0];
sx q[0];
rz(2.9970905) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45134194) q[2];
sx q[2];
rz(-1.1548496) q[2];
sx q[2];
rz(-1.1385663) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0424846) q[1];
sx q[1];
rz(-2.3900552) q[1];
sx q[1];
rz(-1.621975) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68434207) q[3];
sx q[3];
rz(-2.0668233) q[3];
sx q[3];
rz(-3.0471714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8183225) q[2];
sx q[2];
rz(-2.6939836) q[2];
sx q[2];
rz(-1.6935879) q[2];
rz(1.7285534) q[3];
sx q[3];
rz(-1.5331242) q[3];
sx q[3];
rz(-1.9938699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84561658) q[0];
sx q[0];
rz(-0.71346658) q[0];
sx q[0];
rz(3.0699068) q[0];
rz(1.6638311) q[1];
sx q[1];
rz(-1.3641337) q[1];
sx q[1];
rz(2.51547) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7132162) q[0];
sx q[0];
rz(-0.74031836) q[0];
sx q[0];
rz(-1.4367144) q[0];
rz(-1.4914054) q[2];
sx q[2];
rz(-0.65020665) q[2];
sx q[2];
rz(-0.84480594) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0873634) q[1];
sx q[1];
rz(-2.3647671) q[1];
sx q[1];
rz(-1.8190977) q[1];
rz(-pi) q[2];
rz(2.7887129) q[3];
sx q[3];
rz(-1.7418752) q[3];
sx q[3];
rz(1.7069579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46586299) q[2];
sx q[2];
rz(-0.53368038) q[2];
sx q[2];
rz(-1.6443058) q[2];
rz(-2.7346482) q[3];
sx q[3];
rz(-0.99138433) q[3];
sx q[3];
rz(-0.86010325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96754718) q[0];
sx q[0];
rz(-1.8438735) q[0];
sx q[0];
rz(0.26200727) q[0];
rz(1.920248) q[1];
sx q[1];
rz(-1.6294934) q[1];
sx q[1];
rz(1.9409723) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32258275) q[0];
sx q[0];
rz(-2.2457425) q[0];
sx q[0];
rz(-2.2133166) q[0];
x q[1];
rz(2.6083857) q[2];
sx q[2];
rz(-2.2099127) q[2];
sx q[2];
rz(2.5653735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.756304) q[1];
sx q[1];
rz(-1.0693764) q[1];
sx q[1];
rz(2.8629567) q[1];
rz(3.1061078) q[3];
sx q[3];
rz(-1.4789797) q[3];
sx q[3];
rz(-2.6694856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1205341) q[2];
sx q[2];
rz(-0.48212442) q[2];
sx q[2];
rz(-0.47312197) q[2];
rz(2.1038697) q[3];
sx q[3];
rz(-2.7389052) q[3];
sx q[3];
rz(2.0033526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4859908) q[0];
sx q[0];
rz(-0.36360535) q[0];
sx q[0];
rz(0.77014357) q[0];
rz(-2.7262402) q[1];
sx q[1];
rz(-1.5130006) q[1];
sx q[1];
rz(-1.8619246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7466542) q[0];
sx q[0];
rz(-1.5956559) q[0];
sx q[0];
rz(2.3034873) q[0];
rz(0.12721956) q[2];
sx q[2];
rz(-1.9992277) q[2];
sx q[2];
rz(-0.7920533) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7918167) q[1];
sx q[1];
rz(-2.8703051) q[1];
sx q[1];
rz(2.1955745) q[1];
x q[2];
rz(1.5804251) q[3];
sx q[3];
rz(-1.0917679) q[3];
sx q[3];
rz(2.0183394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27560774) q[2];
sx q[2];
rz(-1.393968) q[2];
sx q[2];
rz(2.7371791) q[2];
rz(1.9870997) q[3];
sx q[3];
rz(-1.6136074) q[3];
sx q[3];
rz(1.9563458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6605717) q[0];
sx q[0];
rz(-2.4206929) q[0];
sx q[0];
rz(-0.50203669) q[0];
rz(-1.0188811) q[1];
sx q[1];
rz(-1.6682245) q[1];
sx q[1];
rz(0.92207164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88387892) q[0];
sx q[0];
rz(-1.3259616) q[0];
sx q[0];
rz(1.6689945) q[0];
x q[1];
rz(-2.5174705) q[2];
sx q[2];
rz(-1.3692453) q[2];
sx q[2];
rz(2.1132642) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82220158) q[1];
sx q[1];
rz(-1.5272045) q[1];
sx q[1];
rz(1.6258214) q[1];
x q[2];
rz(-1.5516547) q[3];
sx q[3];
rz(-0.12775207) q[3];
sx q[3];
rz(0.16523339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.33793882) q[2];
sx q[2];
rz(-1.2793845) q[2];
sx q[2];
rz(-3.0214018) q[2];
rz(-3.1386612) q[3];
sx q[3];
rz(-0.35792297) q[3];
sx q[3];
rz(-3.1357989) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1015162) q[0];
sx q[0];
rz(-2.4637971) q[0];
sx q[0];
rz(1.5090322) q[0];
rz(-0.82935968) q[1];
sx q[1];
rz(-0.50193915) q[1];
sx q[1];
rz(-1.1306184) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42989995) q[0];
sx q[0];
rz(-1.9112753) q[0];
sx q[0];
rz(0.32993028) q[0];
rz(-pi) q[1];
rz(1.0605325) q[2];
sx q[2];
rz(-2.2203373) q[2];
sx q[2];
rz(-2.5819957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3106416) q[1];
sx q[1];
rz(-2.0071908) q[1];
sx q[1];
rz(-0.92943014) q[1];
x q[2];
rz(-0.25205125) q[3];
sx q[3];
rz(-1.0672616) q[3];
sx q[3];
rz(-1.2331215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.074038) q[2];
sx q[2];
rz(-2.270547) q[2];
sx q[2];
rz(-1.7769495) q[2];
rz(0.015297628) q[3];
sx q[3];
rz(-2.3099895) q[3];
sx q[3];
rz(-1.3219705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44052723) q[0];
sx q[0];
rz(-2.4897713) q[0];
sx q[0];
rz(0.1524674) q[0];
rz(-1.7561779) q[1];
sx q[1];
rz(-1.0453753) q[1];
sx q[1];
rz(-0.32858953) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4824959) q[0];
sx q[0];
rz(-1.5799111) q[0];
sx q[0];
rz(0.013633058) q[0];
x q[1];
rz(0.25224884) q[2];
sx q[2];
rz(-2.3636732) q[2];
sx q[2];
rz(-1.8053448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8161502) q[1];
sx q[1];
rz(-0.24005664) q[1];
sx q[1];
rz(1.4453959) q[1];
rz(2.4963107) q[3];
sx q[3];
rz(-1.2835763) q[3];
sx q[3];
rz(2.9336799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80455989) q[2];
sx q[2];
rz(-1.1351265) q[2];
sx q[2];
rz(1.0767153) q[2];
rz(-2.1870901) q[3];
sx q[3];
rz(-1.4408305) q[3];
sx q[3];
rz(0.29630989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3267219) q[0];
sx q[0];
rz(-0.9878511) q[0];
sx q[0];
rz(-1.7350154) q[0];
rz(1.8234491) q[1];
sx q[1];
rz(-1.6271918) q[1];
sx q[1];
rz(0.79961332) q[1];
rz(3.0982207) q[2];
sx q[2];
rz(-1.5159173) q[2];
sx q[2];
rz(2.2295502) q[2];
rz(-1.7794505) q[3];
sx q[3];
rz(-0.91051523) q[3];
sx q[3];
rz(-0.12403535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
