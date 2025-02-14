OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0535799) q[0];
sx q[0];
rz(-0.3126643) q[0];
sx q[0];
rz(3.0076658) q[0];
rz(3.1349831) q[1];
sx q[1];
rz(-1.5899038) q[1];
sx q[1];
rz(2.7166727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44533396) q[0];
sx q[0];
rz(-2.9888392) q[0];
sx q[0];
rz(2.2905718) q[0];
rz(-pi) q[1];
rz(-1.0966461) q[2];
sx q[2];
rz(-0.72410781) q[2];
sx q[2];
rz(0.5257789) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1462789) q[1];
sx q[1];
rz(-1.2321207) q[1];
sx q[1];
rz(-2.2168753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4185804) q[3];
sx q[3];
rz(-1.0499081) q[3];
sx q[3];
rz(2.3191444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46140823) q[2];
sx q[2];
rz(-2.7646061) q[2];
sx q[2];
rz(-3.027463) q[2];
rz(-0.38873172) q[3];
sx q[3];
rz(-2.8753493) q[3];
sx q[3];
rz(1.8069327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1330133) q[0];
sx q[0];
rz(-0.40224922) q[0];
sx q[0];
rz(-3.114793) q[0];
rz(1.1773479) q[1];
sx q[1];
rz(-0.64759308) q[1];
sx q[1];
rz(0.037638232) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0018432) q[0];
sx q[0];
rz(-1.5842373) q[0];
sx q[0];
rz(-2.3670908) q[0];
rz(-pi) q[1];
rz(-0.22372795) q[2];
sx q[2];
rz(-1.7729974) q[2];
sx q[2];
rz(0.0082317106) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37477949) q[1];
sx q[1];
rz(-2.1407003) q[1];
sx q[1];
rz(-2.5606696) q[1];
rz(-0.2554727) q[3];
sx q[3];
rz(-2.3049816) q[3];
sx q[3];
rz(-0.59630442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8893975) q[2];
sx q[2];
rz(-2.1571721) q[2];
sx q[2];
rz(-0.21749116) q[2];
rz(-0.43976954) q[3];
sx q[3];
rz(-1.7871126) q[3];
sx q[3];
rz(-0.10194889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.6169823) q[0];
sx q[0];
rz(-3.1341902) q[0];
sx q[0];
rz(-2.5819085) q[0];
rz(0.89645487) q[1];
sx q[1];
rz(-0.47210109) q[1];
sx q[1];
rz(-0.40759531) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893163) q[0];
sx q[0];
rz(-1.6860776) q[0];
sx q[0];
rz(-0.98457054) q[0];
x q[1];
rz(1.4846663) q[2];
sx q[2];
rz(-0.67170947) q[2];
sx q[2];
rz(-1.145878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0478974) q[1];
sx q[1];
rz(-2.0353248) q[1];
sx q[1];
rz(1.4471022) q[1];
x q[2];
rz(-0.30826195) q[3];
sx q[3];
rz(-2.097766) q[3];
sx q[3];
rz(1.7402203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1233623) q[2];
sx q[2];
rz(-1.1574278) q[2];
sx q[2];
rz(1.0250214) q[2];
rz(-1.7068663) q[3];
sx q[3];
rz(-0.44308174) q[3];
sx q[3];
rz(-2.4923435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822815) q[0];
sx q[0];
rz(-0.27114961) q[0];
sx q[0];
rz(0.49736381) q[0];
rz(-1.2936032) q[1];
sx q[1];
rz(-2.135364) q[1];
sx q[1];
rz(1.9088378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3022753) q[0];
sx q[0];
rz(-1.0091821) q[0];
sx q[0];
rz(-1.3297434) q[0];
x q[1];
rz(1.1511939) q[2];
sx q[2];
rz(-0.91165724) q[2];
sx q[2];
rz(-1.731038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9962483) q[1];
sx q[1];
rz(-1.8176907) q[1];
sx q[1];
rz(-1.7801784) q[1];
x q[2];
rz(-3.0972872) q[3];
sx q[3];
rz(-2.7409275) q[3];
sx q[3];
rz(1.6687499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5645912) q[2];
sx q[2];
rz(-0.61508721) q[2];
sx q[2];
rz(-2.9992529) q[2];
rz(1.9650991) q[3];
sx q[3];
rz(-2.1409972) q[3];
sx q[3];
rz(0.4308027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.38554) q[0];
sx q[0];
rz(-2.1591594) q[0];
sx q[0];
rz(0.2734215) q[0];
rz(2.7418819) q[1];
sx q[1];
rz(-2.3276261) q[1];
sx q[1];
rz(-0.25693691) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32460913) q[0];
sx q[0];
rz(-0.63623172) q[0];
sx q[0];
rz(1.0093635) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24946282) q[2];
sx q[2];
rz(-1.9469065) q[2];
sx q[2];
rz(0.94406908) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2013984) q[1];
sx q[1];
rz(-1.640854) q[1];
sx q[1];
rz(0.045658535) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2502348) q[3];
sx q[3];
rz(-1.3472424) q[3];
sx q[3];
rz(-2.4628061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5895245) q[2];
sx q[2];
rz(-2.7858211) q[2];
sx q[2];
rz(0.75999981) q[2];
rz(2.149557) q[3];
sx q[3];
rz(-1.9325247) q[3];
sx q[3];
rz(1.9973283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29042596) q[0];
sx q[0];
rz(-2.5100584) q[0];
sx q[0];
rz(-0.60612154) q[0];
rz(1.0468613) q[1];
sx q[1];
rz(-1.650834) q[1];
sx q[1];
rz(0.077233888) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7116878) q[0];
sx q[0];
rz(-3.0992134) q[0];
sx q[0];
rz(1.2689356) q[0];
rz(0.31862835) q[2];
sx q[2];
rz(-2.0214404) q[2];
sx q[2];
rz(-0.81540996) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97925767) q[1];
sx q[1];
rz(-2.013315) q[1];
sx q[1];
rz(0.44404946) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24155946) q[3];
sx q[3];
rz(-2.4643341) q[3];
sx q[3];
rz(-2.9100071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2250526) q[2];
sx q[2];
rz(-0.70987916) q[2];
sx q[2];
rz(-0.27099657) q[2];
rz(-2.5535876) q[3];
sx q[3];
rz(-0.80395144) q[3];
sx q[3];
rz(-2.7325381) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8005017) q[0];
sx q[0];
rz(-1.611447) q[0];
sx q[0];
rz(-0.28277582) q[0];
rz(-0.68280363) q[1];
sx q[1];
rz(-0.86137259) q[1];
sx q[1];
rz(0.85321325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.565292) q[0];
sx q[0];
rz(-1.1992264) q[0];
sx q[0];
rz(-2.316733) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47854418) q[2];
sx q[2];
rz(-0.87018064) q[2];
sx q[2];
rz(1.2532013) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0921923) q[1];
sx q[1];
rz(-2.8554648) q[1];
sx q[1];
rz(2.186568) q[1];
x q[2];
rz(2.1370579) q[3];
sx q[3];
rz(-2.9101203) q[3];
sx q[3];
rz(-1.8799409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0736488) q[2];
sx q[2];
rz(-0.51764071) q[2];
sx q[2];
rz(-0.29120564) q[2];
rz(-1.9368885) q[3];
sx q[3];
rz(-2.5681345) q[3];
sx q[3];
rz(-1.5706435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.7912927) q[0];
sx q[0];
rz(-1.6181823) q[0];
sx q[0];
rz(2.5147901) q[0];
rz(1.0621915) q[1];
sx q[1];
rz(-2.3419582) q[1];
sx q[1];
rz(-0.83745426) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8663794) q[0];
sx q[0];
rz(-1.4555706) q[0];
sx q[0];
rz(-1.1408477) q[0];
rz(-pi) q[1];
rz(1.3344572) q[2];
sx q[2];
rz(-1.2405985) q[2];
sx q[2];
rz(-0.2969674) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0272701) q[1];
sx q[1];
rz(-2.6114846) q[1];
sx q[1];
rz(-1.8370173) q[1];
rz(1.7193878) q[3];
sx q[3];
rz(-2.0839543) q[3];
sx q[3];
rz(-3.0363014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.31967878) q[2];
sx q[2];
rz(-1.2205114) q[2];
sx q[2];
rz(1.0411881) q[2];
rz(0.90211165) q[3];
sx q[3];
rz(-1.9769042) q[3];
sx q[3];
rz(2.5743217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5015471) q[0];
sx q[0];
rz(-2.8431659) q[0];
sx q[0];
rz(0.10064594) q[0];
rz(1.3238662) q[1];
sx q[1];
rz(-2.3034425) q[1];
sx q[1];
rz(0.59555882) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32257358) q[0];
sx q[0];
rz(-2.3710476) q[0];
sx q[0];
rz(0.84011232) q[0];
rz(1.3084437) q[2];
sx q[2];
rz(-2.4624918) q[2];
sx q[2];
rz(-2.6237112) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42587316) q[1];
sx q[1];
rz(-1.8173479) q[1];
sx q[1];
rz(-0.73664011) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.727739) q[3];
sx q[3];
rz(-2.3683386) q[3];
sx q[3];
rz(-2.2048304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14463921) q[2];
sx q[2];
rz(-2.3914631) q[2];
sx q[2];
rz(-1.3422802) q[2];
rz(-0.73623776) q[3];
sx q[3];
rz(-0.33595556) q[3];
sx q[3];
rz(-0.098522447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017224273) q[0];
sx q[0];
rz(-2.4898744) q[0];
sx q[0];
rz(2.8454054) q[0];
rz(-1.724297) q[1];
sx q[1];
rz(-0.92767757) q[1];
sx q[1];
rz(-0.16253026) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86315853) q[0];
sx q[0];
rz(-1.5774283) q[0];
sx q[0];
rz(1.5445821) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35803609) q[2];
sx q[2];
rz(-2.0581783) q[2];
sx q[2];
rz(1.0768989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49419241) q[1];
sx q[1];
rz(-1.7583349) q[1];
sx q[1];
rz(1.4306551) q[1];
rz(1.0034836) q[3];
sx q[3];
rz(-2.0334019) q[3];
sx q[3];
rz(-1.0779276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6235003) q[2];
sx q[2];
rz(-1.9237498) q[2];
sx q[2];
rz(-2.6859542) q[2];
rz(1.0738922) q[3];
sx q[3];
rz(-2.5116601) q[3];
sx q[3];
rz(2.5560801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0475273) q[0];
sx q[0];
rz(-1.8210664) q[0];
sx q[0];
rz(-0.6022712) q[0];
rz(-0.31996721) q[1];
sx q[1];
rz(-2.0260369) q[1];
sx q[1];
rz(1.8021348) q[1];
rz(-1.510965) q[2];
sx q[2];
rz(-1.2345805) q[2];
sx q[2];
rz(2.6474093) q[2];
rz(-2.3665041) q[3];
sx q[3];
rz(-1.6476484) q[3];
sx q[3];
rz(-2.021029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
