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
rz(2.7754468) q[0];
rz(2.0868299) q[1];
sx q[1];
rz(-1.8572448) q[1];
sx q[1];
rz(-0.73959124) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6897637) q[0];
sx q[0];
rz(-1.1146995) q[0];
sx q[0];
rz(-0.39556894) q[0];
x q[1];
rz(-2.3984564) q[2];
sx q[2];
rz(-0.78746591) q[2];
sx q[2];
rz(-2.1894388) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.057159) q[1];
sx q[1];
rz(-0.68156238) q[1];
sx q[1];
rz(-0.8820028) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2991907) q[3];
sx q[3];
rz(-1.8823307) q[3];
sx q[3];
rz(0.91115644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.71020484) q[2];
sx q[2];
rz(-1.4375765) q[2];
sx q[2];
rz(-2.0598742) q[2];
rz(-1.2181351) q[3];
sx q[3];
rz(-1.5617322) q[3];
sx q[3];
rz(-1.6703828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.1350988) q[0];
sx q[0];
rz(-2.3068937) q[0];
sx q[0];
rz(2.304402) q[0];
rz(0.93765014) q[1];
sx q[1];
rz(-1.7559914) q[1];
sx q[1];
rz(3.0167276) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9287024) q[0];
sx q[0];
rz(-2.4576408) q[0];
sx q[0];
rz(2.4407843) q[0];
x q[1];
rz(-1.2697095) q[2];
sx q[2];
rz(-1.0137179) q[2];
sx q[2];
rz(2.2473492) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13990046) q[1];
sx q[1];
rz(-0.6065953) q[1];
sx q[1];
rz(-0.14735518) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0286507) q[3];
sx q[3];
rz(-0.30755755) q[3];
sx q[3];
rz(1.2483734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79895926) q[2];
sx q[2];
rz(-2.0755167) q[2];
sx q[2];
rz(-0.40668818) q[2];
rz(-2.0090571) q[3];
sx q[3];
rz(-1.5187902) q[3];
sx q[3];
rz(0.86743814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7750074) q[0];
sx q[0];
rz(-0.58554119) q[0];
sx q[0];
rz(-1.0803692) q[0];
rz(0.22251546) q[1];
sx q[1];
rz(-2.3718926) q[1];
sx q[1];
rz(-0.57058191) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8195651) q[0];
sx q[0];
rz(-1.8796258) q[0];
sx q[0];
rz(1.9714799) q[0];
rz(0.81265575) q[2];
sx q[2];
rz(-0.77746292) q[2];
sx q[2];
rz(-1.2948546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92237207) q[1];
sx q[1];
rz(-2.3736694) q[1];
sx q[1];
rz(1.6025387) q[1];
x q[2];
rz(-2.5665652) q[3];
sx q[3];
rz(-1.7238657) q[3];
sx q[3];
rz(1.1014155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.83501619) q[2];
sx q[2];
rz(-1.6604275) q[2];
sx q[2];
rz(1.2349077) q[2];
rz(-1.0770477) q[3];
sx q[3];
rz(-1.5826694) q[3];
sx q[3];
rz(-2.7243015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0017589105) q[0];
sx q[0];
rz(-1.7708906) q[0];
sx q[0];
rz(0.68029252) q[0];
rz(1.744005) q[1];
sx q[1];
rz(-1.1355419) q[1];
sx q[1];
rz(0.23522338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0455189) q[0];
sx q[0];
rz(-1.0313202) q[0];
sx q[0];
rz(1.1474613) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3346457) q[2];
sx q[2];
rz(-2.2586933) q[2];
sx q[2];
rz(-2.2918224) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.75896133) q[1];
sx q[1];
rz(-1.3516892) q[1];
sx q[1];
rz(0.58398881) q[1];
rz(-pi) q[2];
rz(-1.6793988) q[3];
sx q[3];
rz(-1.83654) q[3];
sx q[3];
rz(0.39963978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32517165) q[2];
sx q[2];
rz(-1.07594) q[2];
sx q[2];
rz(-1.761033) q[2];
rz(-2.374968) q[3];
sx q[3];
rz(-2.4216757) q[3];
sx q[3];
rz(0.64002526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470873) q[0];
sx q[0];
rz(-2.6386059) q[0];
sx q[0];
rz(-1.8736396) q[0];
rz(2.2044334) q[1];
sx q[1];
rz(-2.2068534) q[1];
sx q[1];
rz(-0.25757214) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1510096) q[0];
sx q[0];
rz(-2.5611612) q[0];
sx q[0];
rz(-1.2925757) q[0];
x q[1];
rz(2.0054162) q[2];
sx q[2];
rz(-2.4939257) q[2];
sx q[2];
rz(2.2369564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57818164) q[1];
sx q[1];
rz(-1.7817307) q[1];
sx q[1];
rz(0.45244777) q[1];
rz(-pi) q[2];
rz(2.9048728) q[3];
sx q[3];
rz(-1.0977931) q[3];
sx q[3];
rz(1.6110171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.81689721) q[2];
sx q[2];
rz(-2.2163053) q[2];
sx q[2];
rz(2.2980105) q[2];
rz(2.2364565) q[3];
sx q[3];
rz(-2.2456808) q[3];
sx q[3];
rz(-0.26346537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0819241) q[0];
sx q[0];
rz(-0.34316871) q[0];
sx q[0];
rz(-1.4362417) q[0];
rz(-1.2417271) q[1];
sx q[1];
rz(-1.4271586) q[1];
sx q[1];
rz(2.5232975) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.462042) q[0];
sx q[0];
rz(-2.3449595) q[0];
sx q[0];
rz(2.3737143) q[0];
rz(-pi) q[1];
rz(-0.16538362) q[2];
sx q[2];
rz(-0.895154) q[2];
sx q[2];
rz(0.4173511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8346295) q[1];
sx q[1];
rz(-1.8719881) q[1];
sx q[1];
rz(0.42940112) q[1];
rz(-pi) q[2];
rz(-1.854784) q[3];
sx q[3];
rz(-2.0476488) q[3];
sx q[3];
rz(0.43343024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5228086) q[2];
sx q[2];
rz(-1.8960543) q[2];
sx q[2];
rz(1.9469384) q[2];
rz(1.4905802) q[3];
sx q[3];
rz(-1.4158764) q[3];
sx q[3];
rz(1.850261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8959344) q[0];
sx q[0];
rz(-2.7722562) q[0];
sx q[0];
rz(-2.9781407) q[0];
rz(-0.55511904) q[1];
sx q[1];
rz(-1.4733543) q[1];
sx q[1];
rz(2.3177573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6763944) q[0];
sx q[0];
rz(-1.4270743) q[0];
sx q[0];
rz(-2.8222047) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78764913) q[2];
sx q[2];
rz(-1.3795492) q[2];
sx q[2];
rz(-2.8764962) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.026161748) q[1];
sx q[1];
rz(-0.49163252) q[1];
sx q[1];
rz(-0.38422725) q[1];
x q[2];
rz(2.4753544) q[3];
sx q[3];
rz(-2.6347646) q[3];
sx q[3];
rz(3.0493065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0550363) q[2];
sx q[2];
rz(-1.3056511) q[2];
sx q[2];
rz(-1.4646863) q[2];
rz(0.054051789) q[3];
sx q[3];
rz(-2.2650104) q[3];
sx q[3];
rz(0.42266735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5257877) q[0];
sx q[0];
rz(-0.38830385) q[0];
sx q[0];
rz(-1.4317321) q[0];
rz(1.2563541) q[1];
sx q[1];
rz(-1.6836555) q[1];
sx q[1];
rz(-1.8690522) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22879951) q[0];
sx q[0];
rz(-1.4230898) q[0];
sx q[0];
rz(0.65346752) q[0];
rz(2.6600443) q[2];
sx q[2];
rz(-2.5805743) q[2];
sx q[2];
rz(2.53538) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8426395) q[1];
sx q[1];
rz(-1.2455518) q[1];
sx q[1];
rz(-0.3463604) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.921659) q[3];
sx q[3];
rz(-0.84515709) q[3];
sx q[3];
rz(1.2540639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58497477) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2105836) q[0];
sx q[0];
rz(-1.9887661) q[0];
sx q[0];
rz(0.3535122) q[0];
rz(1.968169) q[1];
sx q[1];
rz(-1.2856154) q[1];
sx q[1];
rz(1.8566424) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.374732) q[0];
sx q[0];
rz(-2.6832125) q[0];
sx q[0];
rz(-1.8863669) q[0];
x q[1];
rz(-1.1750269) q[2];
sx q[2];
rz(-1.9916196) q[2];
sx q[2];
rz(1.2683587) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.358153) q[1];
sx q[1];
rz(-0.29488161) q[1];
sx q[1];
rz(1.6452559) q[1];
rz(-pi) q[2];
rz(2.9049614) q[3];
sx q[3];
rz(-2.5116205) q[3];
sx q[3];
rz(-3.0108242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6080007) q[2];
sx q[2];
rz(-2.3054391) q[2];
sx q[2];
rz(-2.3261435) q[2];
rz(-2.4943374) q[3];
sx q[3];
rz(-1.7630968) q[3];
sx q[3];
rz(2.3599724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6668929) q[0];
sx q[0];
rz(-0.26645461) q[0];
sx q[0];
rz(-1.7099963) q[0];
rz(-2.73009) q[1];
sx q[1];
rz(-1.9930379) q[1];
sx q[1];
rz(3.0130951) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602101) q[0];
sx q[0];
rz(-1.8888705) q[0];
sx q[0];
rz(2.4583865) q[0];
rz(-pi) q[1];
rz(0.13984404) q[2];
sx q[2];
rz(-2.2763414) q[2];
sx q[2];
rz(-1.9391413) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.576927) q[1];
sx q[1];
rz(-0.91178545) q[1];
sx q[1];
rz(0.66522775) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5271051) q[3];
sx q[3];
rz(-1.2087421) q[3];
sx q[3];
rz(-1.3339588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0782042) q[2];
sx q[2];
rz(-1.3271164) q[2];
sx q[2];
rz(2.3627031) q[2];
rz(1.286346) q[3];
sx q[3];
rz(-2.0467919) q[3];
sx q[3];
rz(1.26545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.417199) q[0];
sx q[0];
rz(-1.2204285) q[0];
sx q[0];
rz(1.3577419) q[0];
rz(0.67486528) q[1];
sx q[1];
rz(-1.4874896) q[1];
sx q[1];
rz(1.0543324) q[1];
rz(2.7401926) q[2];
sx q[2];
rz(-1.6673604) q[2];
sx q[2];
rz(2.7756804) q[2];
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
