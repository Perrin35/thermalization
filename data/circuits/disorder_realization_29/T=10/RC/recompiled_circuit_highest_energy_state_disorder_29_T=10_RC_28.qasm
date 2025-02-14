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
rz(-2.2588377) q[0];
sx q[0];
rz(-0.14482276) q[0];
sx q[0];
rz(-0.75135279) q[0];
rz(-1.6169647) q[1];
sx q[1];
rz(5.3137988) q[1];
sx q[1];
rz(9.245524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061926024) q[0];
sx q[0];
rz(-2.5653337) q[0];
sx q[0];
rz(-1.040333) q[0];
x q[1];
rz(3.1044311) q[2];
sx q[2];
rz(-2.4183309) q[2];
sx q[2];
rz(-0.0497555) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7764531) q[1];
sx q[1];
rz(-1.6426597) q[1];
sx q[1];
rz(-0.30993575) q[1];
rz(1.4901913) q[3];
sx q[3];
rz(-2.5063519) q[3];
sx q[3];
rz(-2.1135534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90193191) q[2];
sx q[2];
rz(-1.8472981) q[2];
sx q[2];
rz(0.61061668) q[2];
rz(-2.2527952) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(0.73386598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4438542) q[0];
sx q[0];
rz(-2.839851) q[0];
sx q[0];
rz(-1.8119716) q[0];
rz(-1.4854206) q[1];
sx q[1];
rz(-1.679436) q[1];
sx q[1];
rz(0.86404538) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1631016) q[0];
sx q[0];
rz(-1.8231043) q[0];
sx q[0];
rz(-0.10575328) q[0];
rz(2.1435166) q[2];
sx q[2];
rz(-0.53600271) q[2];
sx q[2];
rz(-1.052945) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5031918) q[1];
sx q[1];
rz(-1.4309037) q[1];
sx q[1];
rz(1.3963455) q[1];
x q[2];
rz(-0.95324272) q[3];
sx q[3];
rz(-2.2526178) q[3];
sx q[3];
rz(0.98495959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0580505) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(2.9502499) q[2];
rz(2.8355016) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(-2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8323583) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(-2.7171296) q[0];
rz(-1.3407432) q[1];
sx q[1];
rz(-2.2258591) q[1];
sx q[1];
rz(-0.65139604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0754335) q[0];
sx q[0];
rz(-1.4071583) q[0];
sx q[0];
rz(-3.1343824) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3878787) q[2];
sx q[2];
rz(-1.1929034) q[2];
sx q[2];
rz(-1.479508) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0682023) q[1];
sx q[1];
rz(-1.7531839) q[1];
sx q[1];
rz(0.50478151) q[1];
rz(-0.15768361) q[3];
sx q[3];
rz(-1.7947949) q[3];
sx q[3];
rz(-1.7832613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.035630781) q[2];
sx q[2];
rz(-2.0023846) q[2];
sx q[2];
rz(0.6967217) q[2];
rz(0.68909711) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(2.885163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7032787) q[0];
sx q[0];
rz(-0.2916446) q[0];
sx q[0];
rz(-2.1412204) q[0];
rz(1.4601624) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(1.9283074) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5942237) q[0];
sx q[0];
rz(-1.4495736) q[0];
sx q[0];
rz(1.8531043) q[0];
x q[1];
rz(2.6268509) q[2];
sx q[2];
rz(-0.79823433) q[2];
sx q[2];
rz(1.2003743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.225093) q[1];
sx q[1];
rz(-1.8986393) q[1];
sx q[1];
rz(1.4429379) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.090661006) q[3];
sx q[3];
rz(-2.1075776) q[3];
sx q[3];
rz(1.2556374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0393684) q[2];
sx q[2];
rz(-2.3185456) q[2];
sx q[2];
rz(-3.0774806) q[2];
rz(2.4168849) q[3];
sx q[3];
rz(-1.2084081) q[3];
sx q[3];
rz(-2.7418315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475912) q[0];
sx q[0];
rz(-1.2971017) q[0];
sx q[0];
rz(0.79175788) q[0];
rz(-1.1119615) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.3349104) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9324786) q[0];
sx q[0];
rz(-2.1585585) q[0];
sx q[0];
rz(-2.5665119) q[0];
rz(-pi) q[1];
rz(2.9826775) q[2];
sx q[2];
rz(-0.58664413) q[2];
sx q[2];
rz(1.9288837) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2446652) q[1];
sx q[1];
rz(-1.4003955) q[1];
sx q[1];
rz(3.0651613) q[1];
rz(-1.2615292) q[3];
sx q[3];
rz(-2.5430508) q[3];
sx q[3];
rz(-0.33613294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2935334) q[2];
sx q[2];
rz(-2.2424825) q[2];
sx q[2];
rz(1.1406356) q[2];
rz(-0.10661495) q[3];
sx q[3];
rz(-1.6116319) q[3];
sx q[3];
rz(0.30985668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3265729) q[0];
sx q[0];
rz(-0.44875479) q[0];
sx q[0];
rz(3.1383681) q[0];
rz(-3.0942753) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(-1.3501732) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91302204) q[0];
sx q[0];
rz(-1.6294894) q[0];
sx q[0];
rz(1.8236158) q[0];
rz(-2.2590738) q[2];
sx q[2];
rz(-2.3436758) q[2];
sx q[2];
rz(-1.5882815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.167693) q[1];
sx q[1];
rz(-1.0118139) q[1];
sx q[1];
rz(-2.2476303) q[1];
rz(-pi) q[2];
rz(3.111114) q[3];
sx q[3];
rz(-1.3502096) q[3];
sx q[3];
rz(0.36530803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99001592) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(-0.081175096) q[2];
rz(0.89933991) q[3];
sx q[3];
rz(-0.81493655) q[3];
sx q[3];
rz(-1.8544633) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4261632) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(0.50450605) q[0];
rz(-2.1465178) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(-0.02034932) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58344719) q[0];
sx q[0];
rz(-2.7648395) q[0];
sx q[0];
rz(-1.3135629) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3624914) q[2];
sx q[2];
rz(-1.6951188) q[2];
sx q[2];
rz(0.97036874) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6183747) q[1];
sx q[1];
rz(-1.3058387) q[1];
sx q[1];
rz(-0.19740236) q[1];
rz(-pi) q[2];
rz(-1.4844839) q[3];
sx q[3];
rz(-1.5222957) q[3];
sx q[3];
rz(2.8750471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4703579) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(1.3746369) q[2];
rz(1.6124407) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(3.0488739) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2328211) q[0];
sx q[0];
rz(-0.11255539) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(-2.1807561) q[1];
sx q[1];
rz(-1.4832393) q[1];
sx q[1];
rz(-0.94295162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1652007) q[0];
sx q[0];
rz(-1.2726674) q[0];
sx q[0];
rz(2.0701029) q[0];
rz(-2.9081557) q[2];
sx q[2];
rz(-0.95165157) q[2];
sx q[2];
rz(-0.34282986) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8186989) q[1];
sx q[1];
rz(-2.2320691) q[1];
sx q[1];
rz(-0.7612919) q[1];
rz(1.8887159) q[3];
sx q[3];
rz(-0.32264454) q[3];
sx q[3];
rz(-1.8436197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1711787) q[2];
sx q[2];
rz(-1.3033988) q[2];
sx q[2];
rz(-1.1478434) q[2];
rz(2.7796699) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(-1.743478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4107133) q[0];
sx q[0];
rz(-2.78237) q[0];
sx q[0];
rz(-0.10243375) q[0];
rz(-0.61406413) q[1];
sx q[1];
rz(-2.1153317) q[1];
sx q[1];
rz(-0.817743) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9303904) q[0];
sx q[0];
rz(-1.0387392) q[0];
sx q[0];
rz(-0.802687) q[0];
rz(-pi) q[1];
rz(-0.89985116) q[2];
sx q[2];
rz(-1.4799812) q[2];
sx q[2];
rz(-2.1911774) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7505111) q[1];
sx q[1];
rz(-2.4121248) q[1];
sx q[1];
rz(-1.5424278) q[1];
rz(-pi) q[2];
rz(2.2028186) q[3];
sx q[3];
rz(-1.1891439) q[3];
sx q[3];
rz(-1.3248688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3488591) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(-2.8279772) q[2];
rz(0.46401986) q[3];
sx q[3];
rz(-0.84657621) q[3];
sx q[3];
rz(-0.919842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.2713276) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(-0.23649293) q[0];
rz(2.927921) q[1];
sx q[1];
rz(-2.2387319) q[1];
sx q[1];
rz(2.9170091) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8471223) q[0];
sx q[0];
rz(-2.0531539) q[0];
sx q[0];
rz(0.76089528) q[0];
x q[1];
rz(2.6170391) q[2];
sx q[2];
rz(-2.4894425) q[2];
sx q[2];
rz(-0.99936501) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27222363) q[1];
sx q[1];
rz(-2.096513) q[1];
sx q[1];
rz(1.7030667) q[1];
x q[2];
rz(0.41992374) q[3];
sx q[3];
rz(-2.2444199) q[3];
sx q[3];
rz(2.6633584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97682041) q[2];
sx q[2];
rz(-0.87146622) q[2];
sx q[2];
rz(-2.9898047) q[2];
rz(3.1101036) q[3];
sx q[3];
rz(-1.3969235) q[3];
sx q[3];
rz(0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0781773) q[0];
sx q[0];
rz(-1.6178394) q[0];
sx q[0];
rz(2.7375426) q[0];
rz(-2.0582485) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(2.940098) q[2];
sx q[2];
rz(-1.5200184) q[2];
sx q[2];
rz(1.308174) q[2];
rz(2.8491889) q[3];
sx q[3];
rz(-0.3429827) q[3];
sx q[3];
rz(-0.61957785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
