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
rz(0.1306611) q[0];
sx q[0];
rz(4.9424439) q[0];
sx q[0];
rz(10.688936) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(-1.0882508) q[1];
sx q[1];
rz(0.01297125) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54161763) q[0];
sx q[0];
rz(-2.0429039) q[0];
sx q[0];
rz(1.7469445) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4616724) q[2];
sx q[2];
rz(-1.4089917) q[2];
sx q[2];
rz(-1.8617804) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56264191) q[1];
sx q[1];
rz(-1.2516252) q[1];
sx q[1];
rz(-0.36960543) q[1];
x q[2];
rz(-2.8230328) q[3];
sx q[3];
rz(-0.6069285) q[3];
sx q[3];
rz(-1.7817117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26244792) q[2];
sx q[2];
rz(-1.4538572) q[2];
sx q[2];
rz(0.095414735) q[2];
rz(-2.6573507) q[3];
sx q[3];
rz(-2.3384194) q[3];
sx q[3];
rz(1.6835083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5937623) q[0];
sx q[0];
rz(-1.4365124) q[0];
sx q[0];
rz(-2.4875212) q[0];
rz(-1.3520799) q[1];
sx q[1];
rz(-2.4857931) q[1];
sx q[1];
rz(-0.10993122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2034775) q[0];
sx q[0];
rz(-3.118096) q[0];
sx q[0];
rz(1.0272988) q[0];
rz(-pi) q[1];
rz(2.6835126) q[2];
sx q[2];
rz(-0.9914248) q[2];
sx q[2];
rz(-1.149385) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0866894) q[1];
sx q[1];
rz(-2.6536021) q[1];
sx q[1];
rz(2.7665374) q[1];
rz(1.8255005) q[3];
sx q[3];
rz(-1.9557908) q[3];
sx q[3];
rz(0.32190548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8249417) q[2];
sx q[2];
rz(-1.2085088) q[2];
sx q[2];
rz(1.4671154) q[2];
rz(-1.0351099) q[3];
sx q[3];
rz(-2.3605774) q[3];
sx q[3];
rz(1.8234113) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255945) q[0];
sx q[0];
rz(-2.9075629) q[0];
sx q[0];
rz(-0.79063928) q[0];
rz(2.3979893) q[1];
sx q[1];
rz(-1.598282) q[1];
sx q[1];
rz(0.57949439) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4369219) q[0];
sx q[0];
rz(-2.1050859) q[0];
sx q[0];
rz(2.2493258) q[0];
rz(-pi) q[1];
rz(-2.8786009) q[2];
sx q[2];
rz(-1.8557669) q[2];
sx q[2];
rz(-0.44885744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.070755006) q[1];
sx q[1];
rz(-0.26076128) q[1];
sx q[1];
rz(-1.4255131) q[1];
x q[2];
rz(-1.0813339) q[3];
sx q[3];
rz(-2.0474985) q[3];
sx q[3];
rz(2.7923194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6774595) q[2];
sx q[2];
rz(-0.83659283) q[2];
sx q[2];
rz(0.38132384) q[2];
rz(-0.85121202) q[3];
sx q[3];
rz(-0.92654735) q[3];
sx q[3];
rz(0.73494953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6094991) q[0];
sx q[0];
rz(-0.61436009) q[0];
sx q[0];
rz(2.7771948) q[0];
rz(0.20026194) q[1];
sx q[1];
rz(-1.7297144) q[1];
sx q[1];
rz(-3.0416378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9635506) q[0];
sx q[0];
rz(-3.0638803) q[0];
sx q[0];
rz(1.8366209) q[0];
rz(-pi) q[1];
rz(1.1420629) q[2];
sx q[2];
rz(-1.4754646) q[2];
sx q[2];
rz(0.39458654) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8122319) q[1];
sx q[1];
rz(-1.6633004) q[1];
sx q[1];
rz(1.5483556) q[1];
rz(-0.57860472) q[3];
sx q[3];
rz(-2.1044995) q[3];
sx q[3];
rz(2.5924975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.072307) q[2];
sx q[2];
rz(-2.0733209) q[2];
sx q[2];
rz(3.0266673) q[2];
rz(-1.140444) q[3];
sx q[3];
rz(-1.8585669) q[3];
sx q[3];
rz(-2.1663402) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9370148) q[0];
sx q[0];
rz(-0.95252043) q[0];
sx q[0];
rz(2.9691147) q[0];
rz(-1.0109488) q[1];
sx q[1];
rz(-0.5609678) q[1];
sx q[1];
rz(-0.15484658) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5342543) q[0];
sx q[0];
rz(-1.3528498) q[0];
sx q[0];
rz(1.2675257) q[0];
x q[1];
rz(-2.2151383) q[2];
sx q[2];
rz(-1.7688892) q[2];
sx q[2];
rz(-1.4938172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78945827) q[1];
sx q[1];
rz(-2.704014) q[1];
sx q[1];
rz(0.34559135) q[1];
x q[2];
rz(2.9633767) q[3];
sx q[3];
rz(-1.4608188) q[3];
sx q[3];
rz(0.66457716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.054472063) q[2];
sx q[2];
rz(-2.8643769) q[2];
sx q[2];
rz(-2.7692914) q[2];
rz(-0.39792684) q[3];
sx q[3];
rz(-1.9979265) q[3];
sx q[3];
rz(0.00042375617) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18846866) q[0];
sx q[0];
rz(-0.57250452) q[0];
sx q[0];
rz(1.5484126) q[0];
rz(-2.7717223) q[1];
sx q[1];
rz(-1.3984171) q[1];
sx q[1];
rz(2.5028548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8014378) q[0];
sx q[0];
rz(-1.4169013) q[0];
sx q[0];
rz(1.3305386) q[0];
rz(-pi) q[1];
rz(2.9385826) q[2];
sx q[2];
rz(-2.087321) q[2];
sx q[2];
rz(-2.2882216) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.35013851) q[1];
sx q[1];
rz(-2.0220857) q[1];
sx q[1];
rz(2.7103488) q[1];
rz(-pi) q[2];
rz(-0.79213284) q[3];
sx q[3];
rz(-1.9895305) q[3];
sx q[3];
rz(1.2611539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0729735) q[2];
sx q[2];
rz(-2.0899453) q[2];
sx q[2];
rz(-0.2529141) q[2];
rz(2.2933293) q[3];
sx q[3];
rz(-1.4784644) q[3];
sx q[3];
rz(-0.084913582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0406168) q[0];
sx q[0];
rz(-0.210013) q[0];
sx q[0];
rz(-0.14895359) q[0];
rz(1.3111929) q[1];
sx q[1];
rz(-1.4102178) q[1];
sx q[1];
rz(-0.054100903) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.220517) q[0];
sx q[0];
rz(-2.069325) q[0];
sx q[0];
rz(-1.4709298) q[0];
rz(3.0981423) q[2];
sx q[2];
rz(-2.3386526) q[2];
sx q[2];
rz(0.96722764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91840345) q[1];
sx q[1];
rz(-1.6805873) q[1];
sx q[1];
rz(0.70007433) q[1];
x q[2];
rz(-1.0955174) q[3];
sx q[3];
rz(-1.5915119) q[3];
sx q[3];
rz(1.6894065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3397303) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(-0.96699634) q[2];
rz(2.4407834) q[3];
sx q[3];
rz(-0.98027027) q[3];
sx q[3];
rz(-0.35082671) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12857777) q[0];
sx q[0];
rz(-1.9676493) q[0];
sx q[0];
rz(1.5933734) q[0];
rz(-0.37823996) q[1];
sx q[1];
rz(-2.389237) q[1];
sx q[1];
rz(-2.7517448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3340095) q[0];
sx q[0];
rz(-0.54348677) q[0];
sx q[0];
rz(-2.7568211) q[0];
rz(-1.8748449) q[2];
sx q[2];
rz(-1.5549193) q[2];
sx q[2];
rz(2.7823663) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2521542) q[1];
sx q[1];
rz(-2.5233626) q[1];
sx q[1];
rz(-2.5095224) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7933361) q[3];
sx q[3];
rz(-1.8197818) q[3];
sx q[3];
rz(0.96446645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0491911) q[2];
sx q[2];
rz(-1.1213877) q[2];
sx q[2];
rz(-1.8828877) q[2];
rz(0.24426584) q[3];
sx q[3];
rz(-1.786307) q[3];
sx q[3];
rz(-2.145483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.32023892) q[0];
sx q[0];
rz(-1.2293674) q[0];
sx q[0];
rz(-1.356333) q[0];
rz(2.9755196) q[1];
sx q[1];
rz(-2.1898654) q[1];
sx q[1];
rz(-0.43620268) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8863731) q[0];
sx q[0];
rz(-0.56111911) q[0];
sx q[0];
rz(2.5291689) q[0];
x q[1];
rz(-2.8032254) q[2];
sx q[2];
rz(-1.1778579) q[2];
sx q[2];
rz(1.636105) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.055298565) q[1];
sx q[1];
rz(-2.3196967) q[1];
sx q[1];
rz(0.31676745) q[1];
rz(0.14472503) q[3];
sx q[3];
rz(-2.7008778) q[3];
sx q[3];
rz(-3.0851229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1205552) q[2];
sx q[2];
rz(-1.1641116) q[2];
sx q[2];
rz(2.5980921) q[2];
rz(-0.049526878) q[3];
sx q[3];
rz(-1.6560358) q[3];
sx q[3];
rz(-1.653695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5307584) q[0];
sx q[0];
rz(-0.13228358) q[0];
sx q[0];
rz(0.079205967) q[0];
rz(2.7414956) q[1];
sx q[1];
rz(-0.85405093) q[1];
sx q[1];
rz(-1.0265464) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3883822) q[0];
sx q[0];
rz(-2.3738656) q[0];
sx q[0];
rz(1.9110762) q[0];
rz(-pi) q[1];
rz(1.6658989) q[2];
sx q[2];
rz(-1.9355023) q[2];
sx q[2];
rz(-2.022418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0821918) q[1];
sx q[1];
rz(-2.804356) q[1];
sx q[1];
rz(-0.2880917) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21422106) q[3];
sx q[3];
rz(-0.35839265) q[3];
sx q[3];
rz(0.58979366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2058699) q[2];
sx q[2];
rz(-1.8368072) q[2];
sx q[2];
rz(-1.6892461) q[2];
rz(-0.82990372) q[3];
sx q[3];
rz(-2.2645576) q[3];
sx q[3];
rz(0.95480603) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6185388) q[0];
sx q[0];
rz(-1.1445615) q[0];
sx q[0];
rz(-1.6339697) q[0];
rz(1.2670831) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(-0.53955033) q[2];
sx q[2];
rz(-2.087941) q[2];
sx q[2];
rz(2.6826442) q[2];
rz(-0.55947727) q[3];
sx q[3];
rz(-1.8636129) q[3];
sx q[3];
rz(-1.158798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
