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
rz(-1.0547628) q[1];
sx q[1];
rz(-1.2843479) q[1];
sx q[1];
rz(0.73959124) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30692431) q[0];
sx q[0];
rz(-0.59446228) q[0];
sx q[0];
rz(-0.90499808) q[0];
rz(-pi) q[1];
rz(2.3984564) q[2];
sx q[2];
rz(-2.3541267) q[2];
sx q[2];
rz(0.95215381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.057159) q[1];
sx q[1];
rz(-2.4600303) q[1];
sx q[1];
rz(-2.2595899) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8189956) q[3];
sx q[3];
rz(-1.3125714) q[3];
sx q[3];
rz(-0.74479529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71020484) q[2];
sx q[2];
rz(-1.7040161) q[2];
sx q[2];
rz(-1.0817184) q[2];
rz(-1.2181351) q[3];
sx q[3];
rz(-1.5617322) q[3];
sx q[3];
rz(1.4712099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0064938) q[0];
sx q[0];
rz(-0.83469892) q[0];
sx q[0];
rz(2.304402) q[0];
rz(-2.2039425) q[1];
sx q[1];
rz(-1.3856013) q[1];
sx q[1];
rz(-3.0167276) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9287024) q[0];
sx q[0];
rz(-2.4576408) q[0];
sx q[0];
rz(0.70080832) q[0];
rz(-pi) q[1];
rz(0.44434785) q[2];
sx q[2];
rz(-0.62558657) q[2];
sx q[2];
rz(1.4252961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8320476) q[1];
sx q[1];
rz(-1.654594) q[1];
sx q[1];
rz(-0.60150163) q[1];
rz(-3.0286507) q[3];
sx q[3];
rz(-0.30755755) q[3];
sx q[3];
rz(-1.2483734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3426334) q[2];
sx q[2];
rz(-2.0755167) q[2];
sx q[2];
rz(-2.7349045) q[2];
rz(2.0090571) q[3];
sx q[3];
rz(-1.5187902) q[3];
sx q[3];
rz(-0.86743814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7750074) q[0];
sx q[0];
rz(-2.5560515) q[0];
sx q[0];
rz(-1.0803692) q[0];
rz(2.9190772) q[1];
sx q[1];
rz(-2.3718926) q[1];
sx q[1];
rz(0.57058191) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1207229) q[0];
sx q[0];
rz(-1.9515224) q[0];
sx q[0];
rz(0.33353593) q[0];
rz(-2.1913085) q[2];
sx q[2];
rz(-1.0675021) q[2];
sx q[2];
rz(0.86959834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92237207) q[1];
sx q[1];
rz(-0.76792323) q[1];
sx q[1];
rz(1.6025387) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3889844) q[3];
sx q[3];
rz(-1.003328) q[3];
sx q[3];
rz(0.37088567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83501619) q[2];
sx q[2];
rz(-1.6604275) q[2];
sx q[2];
rz(-1.2349077) q[2];
rz(1.0770477) q[3];
sx q[3];
rz(-1.5589232) q[3];
sx q[3];
rz(-2.7243015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0017589105) q[0];
sx q[0];
rz(-1.3707021) q[0];
sx q[0];
rz(-2.4613001) q[0];
rz(1.744005) q[1];
sx q[1];
rz(-2.0060507) q[1];
sx q[1];
rz(2.9063693) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0455189) q[0];
sx q[0];
rz(-2.1102724) q[0];
sx q[0];
rz(-1.9941313) q[0];
rz(0.3346457) q[2];
sx q[2];
rz(-0.88289936) q[2];
sx q[2];
rz(-2.2918224) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.75896133) q[1];
sx q[1];
rz(-1.3516892) q[1];
sx q[1];
rz(0.58398881) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7626183) q[3];
sx q[3];
rz(-0.28659209) q[3];
sx q[3];
rz(-3.135457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.816421) q[2];
sx q[2];
rz(-2.0656526) q[2];
sx q[2];
rz(1.3805597) q[2];
rz(0.76662463) q[3];
sx q[3];
rz(-0.719917) q[3];
sx q[3];
rz(2.5015674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470873) q[0];
sx q[0];
rz(-0.50298679) q[0];
sx q[0];
rz(-1.8736396) q[0];
rz(-2.2044334) q[1];
sx q[1];
rz(-2.2068534) q[1];
sx q[1];
rz(-2.8840205) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99058305) q[0];
sx q[0];
rz(-0.58043142) q[0];
sx q[0];
rz(1.849017) q[0];
x q[1];
rz(-2.1721971) q[2];
sx q[2];
rz(-1.3139408) q[2];
sx q[2];
rz(1.0207301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7422545) q[1];
sx q[1];
rz(-0.4961001) q[1];
sx q[1];
rz(0.45544405) q[1];
x q[2];
rz(1.0862284) q[3];
sx q[3];
rz(-1.3604829) q[3];
sx q[3];
rz(-0.14968256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81689721) q[2];
sx q[2];
rz(-2.2163053) q[2];
sx q[2];
rz(2.2980105) q[2];
rz(0.90513611) q[3];
sx q[3];
rz(-2.2456808) q[3];
sx q[3];
rz(-2.8781273) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0596685) q[0];
sx q[0];
rz(-0.34316871) q[0];
sx q[0];
rz(-1.7053509) q[0];
rz(-1.8998655) q[1];
sx q[1];
rz(-1.714434) q[1];
sx q[1];
rz(-0.61829511) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.462042) q[0];
sx q[0];
rz(-0.79663314) q[0];
sx q[0];
rz(2.3737143) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3682057) q[2];
sx q[2];
rz(-0.69249047) q[2];
sx q[2];
rz(2.4634374) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8346295) q[1];
sx q[1];
rz(-1.8719881) q[1];
sx q[1];
rz(0.42940112) q[1];
rz(-pi) q[2];
rz(1.854784) q[3];
sx q[3];
rz(-1.0939438) q[3];
sx q[3];
rz(-2.7081624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5228086) q[2];
sx q[2];
rz(-1.8960543) q[2];
sx q[2];
rz(1.1946542) q[2];
rz(-1.6510125) q[3];
sx q[3];
rz(-1.4158764) q[3];
sx q[3];
rz(-1.2913316) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2456583) q[0];
sx q[0];
rz(-0.36933649) q[0];
sx q[0];
rz(-2.9781407) q[0];
rz(0.55511904) q[1];
sx q[1];
rz(-1.4733543) q[1];
sx q[1];
rz(0.82383531) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6763944) q[0];
sx q[0];
rz(-1.4270743) q[0];
sx q[0];
rz(2.8222047) q[0];
rz(2.8749039) q[2];
sx q[2];
rz(-2.3359535) q[2];
sx q[2];
rz(2.0229111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1154309) q[1];
sx q[1];
rz(-2.6499601) q[1];
sx q[1];
rz(-2.7573654) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7300445) q[3];
sx q[3];
rz(-1.8754861) q[3];
sx q[3];
rz(-2.2652486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0865563) q[2];
sx q[2];
rz(-1.3056511) q[2];
sx q[2];
rz(-1.4646863) q[2];
rz(-3.0875409) q[3];
sx q[3];
rz(-0.87658221) q[3];
sx q[3];
rz(-0.42266735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5257877) q[0];
sx q[0];
rz(-0.38830385) q[0];
sx q[0];
rz(-1.7098606) q[0];
rz(1.2563541) q[1];
sx q[1];
rz(-1.6836555) q[1];
sx q[1];
rz(1.2725405) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4542101) q[0];
sx q[0];
rz(-2.2159528) q[0];
sx q[0];
rz(1.3855471) q[0];
rz(-pi) q[1];
x q[1];
rz(1.854004) q[2];
sx q[2];
rz(-1.079756) q[2];
sx q[2];
rz(-1.9824189) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9845674) q[1];
sx q[1];
rz(-1.8982982) q[1];
sx q[1];
rz(-1.9150325) q[1];
rz(-1.921659) q[3];
sx q[3];
rz(-2.2964356) q[3];
sx q[3];
rz(1.8875287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58497477) q[2];
sx q[2];
rz(-1.7375526) q[2];
sx q[2];
rz(-0.72648826) q[2];
rz(-0.11897421) q[3];
sx q[3];
rz(-1.3930895) q[3];
sx q[3];
rz(2.6208641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93100905) q[0];
sx q[0];
rz(-1.1528265) q[0];
sx q[0];
rz(2.7880805) q[0];
rz(1.1734236) q[1];
sx q[1];
rz(-1.2856154) q[1];
sx q[1];
rz(-1.8566424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1160082) q[0];
sx q[0];
rz(-1.1366397) q[0];
sx q[0];
rz(-0.15196073) q[0];
x q[1];
rz(1.9665657) q[2];
sx q[2];
rz(-1.9916196) q[2];
sx q[2];
rz(-1.873234) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8576926) q[1];
sx q[1];
rz(-1.5491747) q[1];
sx q[1];
rz(1.2766854) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23663123) q[3];
sx q[3];
rz(-0.62997216) q[3];
sx q[3];
rz(3.0108242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53359199) q[2];
sx q[2];
rz(-2.3054391) q[2];
sx q[2];
rz(-2.3261435) q[2];
rz(2.4943374) q[3];
sx q[3];
rz(-1.7630968) q[3];
sx q[3];
rz(-2.3599724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4746998) q[0];
sx q[0];
rz(-2.875138) q[0];
sx q[0];
rz(-1.4315963) q[0];
rz(0.41150269) q[1];
sx q[1];
rz(-1.9930379) q[1];
sx q[1];
rz(3.0130951) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14376613) q[0];
sx q[0];
rz(-0.7426922) q[0];
sx q[0];
rz(2.6608442) q[0];
rz(-2.2811878) q[2];
sx q[2];
rz(-1.4644854) q[2];
sx q[2];
rz(-2.8642729) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65662277) q[1];
sx q[1];
rz(-2.2421213) q[1];
sx q[1];
rz(-2.2436618) q[1];
rz(-pi) q[2];
rz(-1.9838215) q[3];
sx q[3];
rz(-1.0810269) q[3];
sx q[3];
rz(-0.44014376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0782042) q[2];
sx q[2];
rz(-1.3271164) q[2];
sx q[2];
rz(0.77888954) q[2];
rz(1.8552467) q[3];
sx q[3];
rz(-2.0467919) q[3];
sx q[3];
rz(-1.26545) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243937) q[0];
sx q[0];
rz(-1.2204285) q[0];
sx q[0];
rz(1.3577419) q[0];
rz(2.4667274) q[1];
sx q[1];
rz(-1.6541031) q[1];
sx q[1];
rz(-2.0872603) q[1];
rz(1.67564) q[2];
sx q[2];
rz(-1.9702199) q[2];
sx q[2];
rz(1.1639845) q[2];
rz(-1.5065083) q[3];
sx q[3];
rz(-2.3953606) q[3];
sx q[3];
rz(2.2904979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
