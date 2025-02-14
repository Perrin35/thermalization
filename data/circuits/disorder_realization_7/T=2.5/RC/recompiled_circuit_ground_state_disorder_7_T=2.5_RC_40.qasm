OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11243842) q[0];
sx q[0];
rz(4.772679) q[0];
sx q[0];
rz(9.7909238) q[0];
rz(2.0868299) q[1];
sx q[1];
rz(4.4259405) q[1];
sx q[1];
rz(8.6851867) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30692431) q[0];
sx q[0];
rz(-2.5471304) q[0];
sx q[0];
rz(2.2365946) q[0];
rz(-pi) q[1];
rz(0.6366836) q[2];
sx q[2];
rz(-1.0708059) q[2];
sx q[2];
rz(1.9477109) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.057159) q[1];
sx q[1];
rz(-2.4600303) q[1];
sx q[1];
rz(-0.8820028) q[1];
x q[2];
rz(1.2991907) q[3];
sx q[3];
rz(-1.259262) q[3];
sx q[3];
rz(2.2304362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71020484) q[2];
sx q[2];
rz(-1.4375765) q[2];
sx q[2];
rz(-1.0817184) q[2];
rz(-1.9234575) q[3];
sx q[3];
rz(-1.5617322) q[3];
sx q[3];
rz(1.6703828) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0064938) q[0];
sx q[0];
rz(-2.3068937) q[0];
sx q[0];
rz(2.304402) q[0];
rz(0.93765014) q[1];
sx q[1];
rz(-1.7559914) q[1];
sx q[1];
rz(-0.12486501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21289028) q[0];
sx q[0];
rz(-0.68395185) q[0];
sx q[0];
rz(0.70080832) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44434785) q[2];
sx q[2];
rz(-0.62558657) q[2];
sx q[2];
rz(1.7162965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8229577) q[1];
sx q[1];
rz(-0.97170107) q[1];
sx q[1];
rz(-1.4692719) q[1];
rz(2.8358744) q[3];
sx q[3];
rz(-1.5366712) q[3];
sx q[3];
rz(-0.21473884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3426334) q[2];
sx q[2];
rz(-2.0755167) q[2];
sx q[2];
rz(-2.7349045) q[2];
rz(-2.0090571) q[3];
sx q[3];
rz(-1.6228024) q[3];
sx q[3];
rz(-0.86743814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36658528) q[0];
sx q[0];
rz(-0.58554119) q[0];
sx q[0];
rz(-1.0803692) q[0];
rz(-0.22251546) q[1];
sx q[1];
rz(-0.76970005) q[1];
sx q[1];
rz(-0.57058191) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1207229) q[0];
sx q[0];
rz(-1.9515224) q[0];
sx q[0];
rz(2.8080567) q[0];
rz(-pi) q[1];
rz(0.9502842) q[2];
sx q[2];
rz(-2.0740905) q[2];
sx q[2];
rz(2.2719943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92237207) q[1];
sx q[1];
rz(-0.76792323) q[1];
sx q[1];
rz(1.6025387) q[1];
x q[2];
rz(0.57502745) q[3];
sx q[3];
rz(-1.7238657) q[3];
sx q[3];
rz(-2.0401772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83501619) q[2];
sx q[2];
rz(-1.6604275) q[2];
sx q[2];
rz(-1.2349077) q[2];
rz(2.0645449) q[3];
sx q[3];
rz(-1.5589232) q[3];
sx q[3];
rz(-0.41729116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1398337) q[0];
sx q[0];
rz(-1.7708906) q[0];
sx q[0];
rz(-2.4613001) q[0];
rz(-1.3975877) q[1];
sx q[1];
rz(-2.0060507) q[1];
sx q[1];
rz(-0.23522338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7021735) q[0];
sx q[0];
rz(-1.930995) q[0];
sx q[0];
rz(0.58106208) q[0];
rz(-pi) q[1];
rz(0.3346457) q[2];
sx q[2];
rz(-2.2586933) q[2];
sx q[2];
rz(2.2918224) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3826313) q[1];
sx q[1];
rz(-1.7899035) q[1];
sx q[1];
rz(0.58398881) q[1];
x q[2];
rz(-0.26724487) q[3];
sx q[3];
rz(-1.4660204) q[3];
sx q[3];
rz(1.9418093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.816421) q[2];
sx q[2];
rz(-1.07594) q[2];
sx q[2];
rz(-1.761033) q[2];
rz(-2.374968) q[3];
sx q[3];
rz(-0.719917) q[3];
sx q[3];
rz(2.5015674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945053) q[0];
sx q[0];
rz(-2.6386059) q[0];
sx q[0];
rz(-1.8736396) q[0];
rz(-0.9371593) q[1];
sx q[1];
rz(-2.2068534) q[1];
sx q[1];
rz(2.8840205) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66143888) q[0];
sx q[0];
rz(-2.1262125) q[0];
sx q[0];
rz(-0.17819779) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0054162) q[2];
sx q[2];
rz(-0.64766696) q[2];
sx q[2];
rz(-2.2369564) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2504006) q[1];
sx q[1];
rz(-1.1290943) q[1];
sx q[1];
rz(1.3370727) q[1];
rz(-2.9048728) q[3];
sx q[3];
rz(-1.0977931) q[3];
sx q[3];
rz(-1.6110171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81689721) q[2];
sx q[2];
rz(-0.92528737) q[2];
sx q[2];
rz(-2.2980105) q[2];
rz(-2.2364565) q[3];
sx q[3];
rz(-0.89591187) q[3];
sx q[3];
rz(2.8781273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0596685) q[0];
sx q[0];
rz(-0.34316871) q[0];
sx q[0];
rz(-1.4362417) q[0];
rz(1.8998655) q[1];
sx q[1];
rz(-1.4271586) q[1];
sx q[1];
rz(-0.61829511) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51793365) q[0];
sx q[0];
rz(-1.0305287) q[0];
sx q[0];
rz(-2.1884657) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88843966) q[2];
sx q[2];
rz(-1.6996146) q[2];
sx q[2];
rz(-1.0494378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8346295) q[1];
sx q[1];
rz(-1.8719881) q[1];
sx q[1];
rz(2.7121915) q[1];
rz(1.2868086) q[3];
sx q[3];
rz(-1.0939438) q[3];
sx q[3];
rz(-0.43343024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5228086) q[2];
sx q[2];
rz(-1.2455384) q[2];
sx q[2];
rz(1.9469384) q[2];
rz(-1.4905802) q[3];
sx q[3];
rz(-1.7257163) q[3];
sx q[3];
rz(1.850261) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8959344) q[0];
sx q[0];
rz(-0.36933649) q[0];
sx q[0];
rz(0.16345197) q[0];
rz(0.55511904) q[1];
sx q[1];
rz(-1.6682383) q[1];
sx q[1];
rz(-0.82383531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15292955) q[0];
sx q[0];
rz(-1.8867765) q[0];
sx q[0];
rz(1.7220604) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8386318) q[2];
sx q[2];
rz(-0.80129708) q[2];
sx q[2];
rz(-1.4943701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8870315) q[1];
sx q[1];
rz(-1.3929092) q[1];
sx q[1];
rz(-0.46079854) q[1];
rz(-pi) q[2];
rz(-1.2402522) q[3];
sx q[3];
rz(-1.179266) q[3];
sx q[3];
rz(0.82465224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0865563) q[2];
sx q[2];
rz(-1.8359416) q[2];
sx q[2];
rz(-1.6769064) q[2];
rz(3.0875409) q[3];
sx q[3];
rz(-2.2650104) q[3];
sx q[3];
rz(-0.42266735) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61580491) q[0];
sx q[0];
rz(-2.7532888) q[0];
sx q[0];
rz(1.4317321) q[0];
rz(1.2563541) q[1];
sx q[1];
rz(-1.6836555) q[1];
sx q[1];
rz(-1.8690522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4542101) q[0];
sx q[0];
rz(-2.2159528) q[0];
sx q[0];
rz(1.3855471) q[0];
rz(-pi) q[1];
rz(-2.6600443) q[2];
sx q[2];
rz(-2.5805743) q[2];
sx q[2];
rz(-2.53538) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8426395) q[1];
sx q[1];
rz(-1.2455518) q[1];
sx q[1];
rz(0.3463604) q[1];
rz(0.36964396) q[3];
sx q[3];
rz(-0.79189205) q[3];
sx q[3];
rz(1.7580851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5566179) q[2];
sx q[2];
rz(-1.40404) q[2];
sx q[2];
rz(2.4151044) q[2];
rz(-3.0226184) q[3];
sx q[3];
rz(-1.7485031) q[3];
sx q[3];
rz(-0.52072853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93100905) q[0];
sx q[0];
rz(-1.9887661) q[0];
sx q[0];
rz(-0.3535122) q[0];
rz(-1.1734236) q[1];
sx q[1];
rz(-1.2856154) q[1];
sx q[1];
rz(-1.2849503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1160082) q[0];
sx q[0];
rz(-2.0049529) q[0];
sx q[0];
rz(-0.15196073) q[0];
rz(-2.6899723) q[2];
sx q[2];
rz(-1.2112144) q[2];
sx q[2];
rz(-0.13338415) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8576926) q[1];
sx q[1];
rz(-1.5491747) q[1];
sx q[1];
rz(-1.8649072) q[1];
rz(-pi) q[2];
rz(-0.6165778) q[3];
sx q[3];
rz(-1.7093466) q[3];
sx q[3];
rz(-1.6324754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6080007) q[2];
sx q[2];
rz(-2.3054391) q[2];
sx q[2];
rz(0.81544915) q[2];
rz(0.6472553) q[3];
sx q[3];
rz(-1.7630968) q[3];
sx q[3];
rz(-0.7816202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4746998) q[0];
sx q[0];
rz(-2.875138) q[0];
sx q[0];
rz(1.4315963) q[0];
rz(0.41150269) q[1];
sx q[1];
rz(-1.1485547) q[1];
sx q[1];
rz(-3.0130951) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0813825) q[0];
sx q[0];
rz(-1.8888705) q[0];
sx q[0];
rz(-2.4583865) q[0];
rz(0.86040489) q[2];
sx q[2];
rz(-1.4644854) q[2];
sx q[2];
rz(0.27731976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.45396475) q[1];
sx q[1];
rz(-1.0610136) q[1];
sx q[1];
rz(-2.3483454) q[1];
rz(0.5271051) q[3];
sx q[3];
rz(-1.2087421) q[3];
sx q[3];
rz(1.3339588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0633885) q[2];
sx q[2];
rz(-1.3271164) q[2];
sx q[2];
rz(-2.3627031) q[2];
rz(-1.286346) q[3];
sx q[3];
rz(-2.0467919) q[3];
sx q[3];
rz(1.8761427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(2.7401926) q[2];
sx q[2];
rz(-1.6673604) q[2];
sx q[2];
rz(2.7756804) q[2];
rz(-1.6350843) q[3];
sx q[3];
rz(-0.74623204) q[3];
sx q[3];
rz(-0.85109477) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
