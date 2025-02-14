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
rz(0.88275498) q[0];
sx q[0];
rz(-2.9967699) q[0];
sx q[0];
rz(0.75135279) q[0];
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
rz(1.0517901) q[0];
sx q[0];
rz(-1.8500916) q[0];
sx q[0];
rz(-2.0816878) q[0];
rz(-0.037161552) q[2];
sx q[2];
rz(-0.72326173) q[2];
sx q[2];
rz(0.0497555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9589267) q[1];
sx q[1];
rz(-1.2616871) q[1];
sx q[1];
rz(1.4953509) q[1];
x q[2];
rz(-3.0823067) q[3];
sx q[3];
rz(-2.2036457) q[3];
sx q[3];
rz(0.92801731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2396607) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(-0.61061668) q[2];
rz(-2.2527952) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(-2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4438542) q[0];
sx q[0];
rz(-2.839851) q[0];
sx q[0];
rz(1.3296211) q[0];
rz(1.4854206) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(0.86404538) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5764424) q[0];
sx q[0];
rz(-0.2731384) q[0];
sx q[0];
rz(-1.9594203) q[0];
rz(-pi) q[1];
rz(2.0338221) q[2];
sx q[2];
rz(-1.8512176) q[2];
sx q[2];
rz(-0.011649557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.043034779) q[1];
sx q[1];
rz(-1.7435257) q[1];
sx q[1];
rz(0.14202001) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1883499) q[3];
sx q[3];
rz(-2.2526178) q[3];
sx q[3];
rz(0.98495959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0835421) q[2];
sx q[2];
rz(-0.6821878) q[2];
sx q[2];
rz(-0.19134276) q[2];
rz(2.8355016) q[3];
sx q[3];
rz(-1.4602665) q[3];
sx q[3];
rz(2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8323583) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(-0.424463) q[0];
rz(1.8008495) q[1];
sx q[1];
rz(-2.2258591) q[1];
sx q[1];
rz(2.4901966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.021928) q[0];
sx q[0];
rz(-2.9777973) q[0];
sx q[0];
rz(-1.6144362) q[0];
rz(-pi) q[1];
rz(-2.7119262) q[2];
sx q[2];
rz(-2.7236807) q[2];
sx q[2];
rz(1.197345) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0682023) q[1];
sx q[1];
rz(-1.7531839) q[1];
sx q[1];
rz(2.6368111) q[1];
rz(-pi) q[2];
rz(0.96730729) q[3];
sx q[3];
rz(-2.8684232) q[3];
sx q[3];
rz(-0.73707132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1059619) q[2];
sx q[2];
rz(-1.1392081) q[2];
sx q[2];
rz(0.6967217) q[2];
rz(-0.68909711) q[3];
sx q[3];
rz(-1.9870116) q[3];
sx q[3];
rz(2.885163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7032787) q[0];
sx q[0];
rz(-0.2916446) q[0];
sx q[0];
rz(-1.0003723) q[0];
rz(1.6814303) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(1.2132852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4184237) q[0];
sx q[0];
rz(-2.8349986) q[0];
sx q[0];
rz(-1.9830389) q[0];
rz(-pi) q[1];
rz(2.0385267) q[2];
sx q[2];
rz(-0.89780318) q[2];
sx q[2];
rz(1.8813934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53663838) q[1];
sx q[1];
rz(-0.35105303) q[1];
sx q[1];
rz(2.7829079) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0509316) q[3];
sx q[3];
rz(-1.034015) q[3];
sx q[3];
rz(-1.2556374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0393684) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(-3.0774806) q[2];
rz(-2.4168849) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(-2.7418315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475912) q[0];
sx q[0];
rz(-1.2971017) q[0];
sx q[0];
rz(0.79175788) q[0];
rz(-2.0296312) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.8066822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34599538) q[0];
sx q[0];
rz(-0.79774081) q[0];
sx q[0];
rz(-2.255385) q[0];
rz(1.6755988) q[2];
sx q[2];
rz(-0.99250472) q[2];
sx q[2];
rz(2.1190018) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89692749) q[1];
sx q[1];
rz(-1.7411971) q[1];
sx q[1];
rz(3.0651613) q[1];
rz(-pi) q[2];
rz(-1.2615292) q[3];
sx q[3];
rz(-2.5430508) q[3];
sx q[3];
rz(2.8054597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84805924) q[2];
sx q[2];
rz(-2.2424825) q[2];
sx q[2];
rz(-2.000957) q[2];
rz(-0.10661495) q[3];
sx q[3];
rz(-1.6116319) q[3];
sx q[3];
rz(0.30985668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3265729) q[0];
sx q[0];
rz(-2.6928379) q[0];
sx q[0];
rz(-3.1383681) q[0];
rz(0.047317304) q[1];
sx q[1];
rz(-1.4553757) q[1];
sx q[1];
rz(-1.7914194) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4989717) q[0];
sx q[0];
rz(-1.823171) q[0];
sx q[0];
rz(-0.060615505) q[0];
rz(-pi) q[1];
rz(-0.88251884) q[2];
sx q[2];
rz(-2.3436758) q[2];
sx q[2];
rz(1.5882815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1941095) q[1];
sx q[1];
rz(-1.0110185) q[1];
sx q[1];
rz(-2.4653788) q[1];
rz(-pi) q[2];
rz(3.111114) q[3];
sx q[3];
rz(-1.7913831) q[3];
sx q[3];
rz(2.7762846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1515767) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(-0.081175096) q[2];
rz(2.2422527) q[3];
sx q[3];
rz(-0.81493655) q[3];
sx q[3];
rz(-1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4261632) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(0.50450605) q[0];
rz(0.99507487) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(3.1212433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85914579) q[0];
sx q[0];
rz(-1.2070281) q[0];
sx q[0];
rz(0.10031853) q[0];
x q[1];
rz(1.3624914) q[2];
sx q[2];
rz(-1.6951188) q[2];
sx q[2];
rz(2.1712239) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.099906057) q[1];
sx q[1];
rz(-1.3803687) q[1];
sx q[1];
rz(-1.8407525) q[1];
rz(-pi) q[2];
rz(1.4844839) q[3];
sx q[3];
rz(-1.619297) q[3];
sx q[3];
rz(-0.26654551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67123479) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(-1.3746369) q[2];
rz(-1.5291519) q[3];
sx q[3];
rz(-2.0917454) q[3];
sx q[3];
rz(-3.0488739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90877157) q[0];
sx q[0];
rz(-3.0290373) q[0];
sx q[0];
rz(-2.0594647) q[0];
rz(2.1807561) q[1];
sx q[1];
rz(-1.4832393) q[1];
sx q[1];
rz(0.94295162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0884224) q[0];
sx q[0];
rz(-0.57500792) q[0];
sx q[0];
rz(1.0002329) q[0];
x q[1];
rz(-2.9081557) q[2];
sx q[2];
rz(-0.95165157) q[2];
sx q[2];
rz(2.7987628) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77737227) q[1];
sx q[1];
rz(-2.1465214) q[1];
sx q[1];
rz(-0.74933021) q[1];
rz(1.2528768) q[3];
sx q[3];
rz(-2.8189481) q[3];
sx q[3];
rz(1.2979729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1711787) q[2];
sx q[2];
rz(-1.3033988) q[2];
sx q[2];
rz(-1.1478434) q[2];
rz(2.7796699) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(1.3981147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73087937) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(-3.0391589) q[0];
rz(0.61406413) q[1];
sx q[1];
rz(-1.026261) q[1];
sx q[1];
rz(2.3238497) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2376643) q[0];
sx q[0];
rz(-2.2127164) q[0];
sx q[0];
rz(-0.68591811) q[0];
rz(-0.89985116) q[2];
sx q[2];
rz(-1.4799812) q[2];
sx q[2];
rz(0.95041529) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9407258) q[1];
sx q[1];
rz(-1.5897017) q[1];
sx q[1];
rz(2.3000642) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2028186) q[3];
sx q[3];
rz(-1.9524487) q[3];
sx q[3];
rz(1.8167239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7927336) q[2];
sx q[2];
rz(-2.0642955) q[2];
sx q[2];
rz(0.31361541) q[2];
rz(-2.6775728) q[3];
sx q[3];
rz(-0.84657621) q[3];
sx q[3];
rz(-0.919842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2713276) q[0];
sx q[0];
rz(-0.79065228) q[0];
sx q[0];
rz(2.9050997) q[0];
rz(-0.21367167) q[1];
sx q[1];
rz(-0.90286076) q[1];
sx q[1];
rz(0.22458354) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29447039) q[0];
sx q[0];
rz(-2.0531539) q[0];
sx q[0];
rz(-0.76089528) q[0];
rz(-pi) q[1];
rz(1.2055254) q[2];
sx q[2];
rz(-1.0177311) q[2];
sx q[2];
rz(1.5129364) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.013063518) q[1];
sx q[1];
rz(-0.54058248) q[1];
sx q[1];
rz(-2.9180727) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41992374) q[3];
sx q[3];
rz(-0.89717275) q[3];
sx q[3];
rz(2.6633584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.97682041) q[2];
sx q[2];
rz(-0.87146622) q[2];
sx q[2];
rz(0.15178794) q[2];
rz(3.1101036) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(-0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0634154) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(2.0582485) q[1];
sx q[1];
rz(-1.5647519) q[1];
sx q[1];
rz(-1.5819989) q[1];
rz(1.622621) q[2];
sx q[2];
rz(-1.3695649) q[2];
sx q[2];
rz(2.8893378) q[2];
rz(-1.4682228) q[3];
sx q[3];
rz(-1.2429308) q[3];
sx q[3];
rz(2.8314005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
