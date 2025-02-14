OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0291542) q[0];
sx q[0];
rz(-1.6310863) q[0];
sx q[0];
rz(-0.36614585) q[0];
rz(2.0868299) q[1];
sx q[1];
rz(4.4259405) q[1];
sx q[1];
rz(8.6851867) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30692431) q[0];
sx q[0];
rz(-2.5471304) q[0];
sx q[0];
rz(-2.2365946) q[0];
rz(0.74313626) q[2];
sx q[2];
rz(-2.3541267) q[2];
sx q[2];
rz(-0.95215381) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.242567) q[1];
sx q[1];
rz(-2.0787313) q[1];
sx q[1];
rz(-2.6655156) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2991907) q[3];
sx q[3];
rz(-1.259262) q[3];
sx q[3];
rz(-0.91115644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4313878) q[2];
sx q[2];
rz(-1.4375765) q[2];
sx q[2];
rz(1.0817184) q[2];
rz(1.2181351) q[3];
sx q[3];
rz(-1.5617322) q[3];
sx q[3];
rz(-1.4712099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1350988) q[0];
sx q[0];
rz(-0.83469892) q[0];
sx q[0];
rz(0.83719069) q[0];
rz(0.93765014) q[1];
sx q[1];
rz(-1.7559914) q[1];
sx q[1];
rz(-0.12486501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77878419) q[0];
sx q[0];
rz(-1.9904525) q[0];
sx q[0];
rz(-2.5843689) q[0];
rz(2.6972448) q[2];
sx q[2];
rz(-0.62558657) q[2];
sx q[2];
rz(-1.4252961) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8229577) q[1];
sx q[1];
rz(-2.1698916) q[1];
sx q[1];
rz(-1.4692719) q[1];
rz(-pi) q[2];
rz(3.0286507) q[3];
sx q[3];
rz(-0.30755755) q[3];
sx q[3];
rz(-1.8932193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79895926) q[2];
sx q[2];
rz(-1.066076) q[2];
sx q[2];
rz(2.7349045) q[2];
rz(1.1325356) q[3];
sx q[3];
rz(-1.5187902) q[3];
sx q[3];
rz(0.86743814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7750074) q[0];
sx q[0];
rz(-2.5560515) q[0];
sx q[0];
rz(1.0803692) q[0];
rz(-0.22251546) q[1];
sx q[1];
rz(-2.3718926) q[1];
sx q[1];
rz(-2.5710107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8195651) q[0];
sx q[0];
rz(-1.8796258) q[0];
sx q[0];
rz(-1.1701128) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81265575) q[2];
sx q[2];
rz(-0.77746292) q[2];
sx q[2];
rz(1.2948546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5160061) q[1];
sx q[1];
rz(-1.592844) q[1];
sx q[1];
rz(-2.3384678) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5665652) q[3];
sx q[3];
rz(-1.7238657) q[3];
sx q[3];
rz(-1.1014155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3065765) q[2];
sx q[2];
rz(-1.6604275) q[2];
sx q[2];
rz(-1.906685) q[2];
rz(-2.0645449) q[3];
sx q[3];
rz(-1.5826694) q[3];
sx q[3];
rz(2.7243015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1398337) q[0];
sx q[0];
rz(-1.3707021) q[0];
sx q[0];
rz(0.68029252) q[0];
rz(1.3975877) q[1];
sx q[1];
rz(-2.0060507) q[1];
sx q[1];
rz(0.23522338) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7021735) q[0];
sx q[0];
rz(-1.2105976) q[0];
sx q[0];
rz(-2.5605306) q[0];
rz(-pi) q[1];
rz(2.806947) q[2];
sx q[2];
rz(-2.2586933) q[2];
sx q[2];
rz(0.84977023) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75896133) q[1];
sx q[1];
rz(-1.7899035) q[1];
sx q[1];
rz(-2.5576038) q[1];
x q[2];
rz(1.4621939) q[3];
sx q[3];
rz(-1.3050526) q[3];
sx q[3];
rz(2.7419529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.816421) q[2];
sx q[2];
rz(-2.0656526) q[2];
sx q[2];
rz(1.761033) q[2];
rz(2.374968) q[3];
sx q[3];
rz(-2.4216757) q[3];
sx q[3];
rz(-0.64002526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0945053) q[0];
sx q[0];
rz(-2.6386059) q[0];
sx q[0];
rz(1.8736396) q[0];
rz(2.2044334) q[1];
sx q[1];
rz(-0.93473923) q[1];
sx q[1];
rz(-2.8840205) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4801538) q[0];
sx q[0];
rz(-2.1262125) q[0];
sx q[0];
rz(2.9633949) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1721971) q[2];
sx q[2];
rz(-1.8276518) q[2];
sx q[2];
rz(-1.0207301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89119205) q[1];
sx q[1];
rz(-2.0124984) q[1];
sx q[1];
rz(1.3370727) q[1];
rz(-2.0553642) q[3];
sx q[3];
rz(-1.7811097) q[3];
sx q[3];
rz(0.14968256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.81689721) q[2];
sx q[2];
rz(-0.92528737) q[2];
sx q[2];
rz(-0.84358215) q[2];
rz(2.2364565) q[3];
sx q[3];
rz(-2.2456808) q[3];
sx q[3];
rz(2.8781273) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0596685) q[0];
sx q[0];
rz(-2.7984239) q[0];
sx q[0];
rz(-1.7053509) q[0];
rz(-1.8998655) q[1];
sx q[1];
rz(-1.714434) q[1];
sx q[1];
rz(2.5232975) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6795507) q[0];
sx q[0];
rz(-0.79663314) q[0];
sx q[0];
rz(2.3737143) q[0];
rz(-pi) q[1];
rz(-2.253153) q[2];
sx q[2];
rz(-1.6996146) q[2];
sx q[2];
rz(-1.0494378) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3028812) q[1];
sx q[1];
rz(-2.6225323) q[1];
sx q[1];
rz(-2.5005591) q[1];
x q[2];
rz(2.6446436) q[3];
sx q[3];
rz(-2.5922311) q[3];
sx q[3];
rz(-3.0086111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5228086) q[2];
sx q[2];
rz(-1.2455384) q[2];
sx q[2];
rz(-1.9469384) q[2];
rz(-1.4905802) q[3];
sx q[3];
rz(-1.7257163) q[3];
sx q[3];
rz(1.850261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2456583) q[0];
sx q[0];
rz(-2.7722562) q[0];
sx q[0];
rz(2.9781407) q[0];
rz(2.5864736) q[1];
sx q[1];
rz(-1.6682383) q[1];
sx q[1];
rz(0.82383531) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15292955) q[0];
sx q[0];
rz(-1.2548162) q[0];
sx q[0];
rz(-1.7220604) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26668875) q[2];
sx q[2];
rz(-0.80563918) q[2];
sx q[2];
rz(-2.0229111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2545611) q[1];
sx q[1];
rz(-1.7486835) q[1];
sx q[1];
rz(-2.6807941) q[1];
x q[2];
rz(1.2402522) q[3];
sx q[3];
rz(-1.9623266) q[3];
sx q[3];
rz(0.82465224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0550363) q[2];
sx q[2];
rz(-1.3056511) q[2];
sx q[2];
rz(1.6769064) q[2];
rz(3.0875409) q[3];
sx q[3];
rz(-0.87658221) q[3];
sx q[3];
rz(-2.7189253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61580491) q[0];
sx q[0];
rz(-2.7532888) q[0];
sx q[0];
rz(1.7098606) q[0];
rz(-1.8852385) q[1];
sx q[1];
rz(-1.6836555) q[1];
sx q[1];
rz(-1.8690522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1521069) q[0];
sx q[0];
rz(-2.4740334) q[0];
sx q[0];
rz(-0.24002536) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50813307) q[2];
sx q[2];
rz(-1.8197803) q[2];
sx q[2];
rz(2.5935885) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9965197) q[1];
sx q[1];
rz(-0.47059083) q[1];
sx q[1];
rz(-0.7820635) q[1];
rz(-pi) q[2];
x q[2];
rz(1.921659) q[3];
sx q[3];
rz(-2.2964356) q[3];
sx q[3];
rz(1.2540639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5566179) q[2];
sx q[2];
rz(-1.7375526) q[2];
sx q[2];
rz(-0.72648826) q[2];
rz(0.11897421) q[3];
sx q[3];
rz(-1.3930895) q[3];
sx q[3];
rz(-2.6208641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93100905) q[0];
sx q[0];
rz(-1.9887661) q[0];
sx q[0];
rz(-2.7880805) q[0];
rz(1.1734236) q[1];
sx q[1];
rz(-1.8559772) q[1];
sx q[1];
rz(1.8566424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0255844) q[0];
sx q[0];
rz(-2.0049529) q[0];
sx q[0];
rz(-2.9896319) q[0];
x q[1];
rz(1.1750269) q[2];
sx q[2];
rz(-1.149973) q[2];
sx q[2];
rz(1.2683587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7834397) q[1];
sx q[1];
rz(-2.846711) q[1];
sx q[1];
rz(1.4963368) q[1];
rz(-pi) q[2];
rz(1.4015163) q[3];
sx q[3];
rz(-0.96099412) q[3];
sx q[3];
rz(2.9823401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53359199) q[2];
sx q[2];
rz(-2.3054391) q[2];
sx q[2];
rz(-2.3261435) q[2];
rz(0.6472553) q[3];
sx q[3];
rz(-1.7630968) q[3];
sx q[3];
rz(-0.7816202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4746998) q[0];
sx q[0];
rz(-2.875138) q[0];
sx q[0];
rz(1.7099963) q[0];
rz(-2.73009) q[1];
sx q[1];
rz(-1.9930379) q[1];
sx q[1];
rz(-0.12849753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0602101) q[0];
sx q[0];
rz(-1.2527221) q[0];
sx q[0];
rz(2.4583865) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13984404) q[2];
sx q[2];
rz(-2.2763414) q[2];
sx q[2];
rz(-1.2024513) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5646657) q[1];
sx q[1];
rz(-0.91178545) q[1];
sx q[1];
rz(0.66522775) q[1];
x q[2];
rz(-1.9838215) q[3];
sx q[3];
rz(-1.0810269) q[3];
sx q[3];
rz(-0.44014376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0633885) q[2];
sx q[2];
rz(-1.8144763) q[2];
sx q[2];
rz(0.77888954) q[2];
rz(-1.286346) q[3];
sx q[3];
rz(-2.0467919) q[3];
sx q[3];
rz(-1.26545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7243937) q[0];
sx q[0];
rz(-1.9211641) q[0];
sx q[0];
rz(-1.7838508) q[0];
rz(-2.4667274) q[1];
sx q[1];
rz(-1.4874896) q[1];
sx q[1];
rz(1.0543324) q[1];
rz(0.24302287) q[2];
sx q[2];
rz(-0.41223787) q[2];
sx q[2];
rz(1.4282474) q[2];
rz(3.0822637) q[3];
sx q[3];
rz(-2.3151201) q[3];
sx q[3];
rz(-0.7636418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
