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
rz(1.524628) q[1];
sx q[1];
rz(-2.1722062) q[1];
sx q[1];
rz(0.17925395) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061926024) q[0];
sx q[0];
rz(-2.5653337) q[0];
sx q[0];
rz(2.1012596) q[0];
rz(-pi) q[1];
rz(-0.037161552) q[2];
sx q[2];
rz(-2.4183309) q[2];
sx q[2];
rz(-0.0497555) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7764531) q[1];
sx q[1];
rz(-1.6426597) q[1];
sx q[1];
rz(-0.30993575) q[1];
rz(-pi) q[2];
rz(0.93710812) q[3];
sx q[3];
rz(-1.5230012) q[3];
sx q[3];
rz(-2.5339047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2396607) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(0.61061668) q[2];
rz(2.2527952) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(-0.73386598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69773847) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(-1.8119716) q[0];
rz(-1.656172) q[1];
sx q[1];
rz(-1.679436) q[1];
sx q[1];
rz(-0.86404538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5764424) q[0];
sx q[0];
rz(-2.8684542) q[0];
sx q[0];
rz(1.1821724) q[0];
rz(1.1077706) q[2];
sx q[2];
rz(-1.290375) q[2];
sx q[2];
rz(3.1299431) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73693528) q[1];
sx q[1];
rz(-0.22316775) q[1];
sx q[1];
rz(0.88921247) q[1];
rz(-pi) q[2];
rz(0.61965539) q[3];
sx q[3];
rz(-2.2562422) q[3];
sx q[3];
rz(1.311613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0580505) q[2];
sx q[2];
rz(-0.6821878) q[2];
sx q[2];
rz(2.9502499) q[2];
rz(0.30609104) q[3];
sx q[3];
rz(-1.4602665) q[3];
sx q[3];
rz(0.93650854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3092344) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(0.424463) q[0];
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
rz(-1.1196647) q[0];
sx q[0];
rz(-0.16379539) q[0];
sx q[0];
rz(-1.6144362) q[0];
rz(-pi) q[1];
rz(-2.7578951) q[2];
sx q[2];
rz(-1.7406775) q[2];
sx q[2];
rz(-3.1184514) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0682023) q[1];
sx q[1];
rz(-1.3884087) q[1];
sx q[1];
rz(0.50478151) q[1];
rz(-pi) q[2];
x q[2];
rz(1.797514) q[3];
sx q[3];
rz(-1.4170839) q[3];
sx q[3];
rz(0.24776974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.035630781) q[2];
sx q[2];
rz(-2.0023846) q[2];
sx q[2];
rz(-2.4448709) q[2];
rz(-0.68909711) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(0.2564297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-1.5337475) q[1];
sx q[1];
rz(-1.2132852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98836658) q[0];
sx q[0];
rz(-1.8509764) q[0];
sx q[0];
rz(3.0154254) q[0];
x q[1];
rz(0.51474173) q[2];
sx q[2];
rz(-0.79823433) q[2];
sx q[2];
rz(1.9412184) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53663838) q[1];
sx q[1];
rz(-0.35105303) q[1];
sx q[1];
rz(-2.7829079) q[1];
x q[2];
rz(0.090661006) q[3];
sx q[3];
rz(-1.034015) q[3];
sx q[3];
rz(1.2556374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1022243) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(-3.0774806) q[2];
rz(-0.72470775) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(-0.39976111) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475912) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(2.3498348) q[0];
rz(-2.0296312) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.8066822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2091141) q[0];
sx q[0];
rz(-2.1585585) q[0];
sx q[0];
rz(0.57508075) q[0];
rz(-pi) q[1];
rz(-0.58081268) q[2];
sx q[2];
rz(-1.6585095) q[2];
sx q[2];
rz(2.650819) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47274703) q[1];
sx q[1];
rz(-2.9549874) q[1];
sx q[1];
rz(-1.9884459) q[1];
rz(-pi) q[2];
rz(-0.99467268) q[3];
sx q[3];
rz(-1.398456) q[3];
sx q[3];
rz(-1.4927499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84805924) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(2.000957) q[2];
rz(0.10661495) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(-2.831736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3265729) q[0];
sx q[0];
rz(-0.44875479) q[0];
sx q[0];
rz(3.1383681) q[0];
rz(3.0942753) q[1];
sx q[1];
rz(-1.4553757) q[1];
sx q[1];
rz(1.7914194) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88104507) q[0];
sx q[0];
rz(-2.8821917) q[0];
sx q[0];
rz(-1.3400643) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5642942) q[2];
sx q[2];
rz(-2.1566628) q[2];
sx q[2];
rz(-0.68638869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1941095) q[1];
sx q[1];
rz(-2.1305741) q[1];
sx q[1];
rz(2.4653788) q[1];
rz(-pi) q[2];
rz(-1.7914823) q[3];
sx q[3];
rz(-1.5410564) q[3];
sx q[3];
rz(1.1988175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99001592) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4261632) q[0];
sx q[0];
rz(-1.8303215) q[0];
sx q[0];
rz(0.50450605) q[0];
rz(-2.1465178) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(3.1212433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5581455) q[0];
sx q[0];
rz(-2.7648395) q[0];
sx q[0];
rz(-1.3135629) q[0];
rz(1.7791012) q[2];
sx q[2];
rz(-1.6951188) q[2];
sx q[2];
rz(-2.1712239) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0416866) q[1];
sx q[1];
rz(-1.3803687) q[1];
sx q[1];
rz(1.3008402) q[1];
rz(-1.6571088) q[3];
sx q[3];
rz(-1.619297) q[3];
sx q[3];
rz(-0.26654551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67123479) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(-1.7669558) q[2];
rz(-1.5291519) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(-0.092718743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90877157) q[0];
sx q[0];
rz(-0.11255539) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(2.1807561) q[1];
sx q[1];
rz(-1.4832393) q[1];
sx q[1];
rz(-2.198641) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0884224) q[0];
sx q[0];
rz(-0.57500792) q[0];
sx q[0];
rz(-2.1413598) q[0];
rz(-1.2569179) q[2];
sx q[2];
rz(-2.4853443) q[2];
sx q[2];
rz(-0.73168025) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3642204) q[1];
sx q[1];
rz(-0.99507123) q[1];
sx q[1];
rz(-0.74933021) q[1];
x q[2];
rz(0.10412962) q[3];
sx q[3];
rz(-1.2648598) q[3];
sx q[3];
rz(-0.96398523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97041398) q[2];
sx q[2];
rz(-1.8381939) q[2];
sx q[2];
rz(1.1478434) q[2];
rz(0.36192274) q[3];
sx q[3];
rz(-1.3779093) q[3];
sx q[3];
rz(-1.743478) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4107133) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(-3.0391589) q[0];
rz(2.5275285) q[1];
sx q[1];
rz(-2.1153317) q[1];
sx q[1];
rz(-0.817743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9303904) q[0];
sx q[0];
rz(-1.0387392) q[0];
sx q[0];
rz(0.802687) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2417415) q[2];
sx q[2];
rz(-1.6616115) q[2];
sx q[2];
rz(0.95041529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7885544) q[1];
sx q[1];
rz(-0.8416881) q[1];
sx q[1];
rz(-0.025351449) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46157591) q[3];
sx q[3];
rz(-2.1510923) q[3];
sx q[3];
rz(-3.1212774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7927336) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(0.31361541) q[2];
rz(2.6775728) q[3];
sx q[3];
rz(-2.2950164) q[3];
sx q[3];
rz(2.2217506) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2713276) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(2.9050997) q[0];
rz(-0.21367167) q[1];
sx q[1];
rz(-2.2387319) q[1];
sx q[1];
rz(-0.22458354) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8471223) q[0];
sx q[0];
rz(-1.0884388) q[0];
sx q[0];
rz(-0.76089528) q[0];
x q[1];
rz(-0.58402365) q[2];
sx q[2];
rz(-1.2619602) q[2];
sx q[2];
rz(-0.14039224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.013063518) q[1];
sx q[1];
rz(-2.6010102) q[1];
sx q[1];
rz(0.22352) q[1];
rz(-pi) q[2];
rz(-2.7216689) q[3];
sx q[3];
rz(-0.89717275) q[3];
sx q[3];
rz(-2.6633584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.97682041) q[2];
sx q[2];
rz(-2.2701264) q[2];
sx q[2];
rz(-2.9898047) q[2];
rz(0.031489059) q[3];
sx q[3];
rz(-1.3969235) q[3];
sx q[3];
rz(-0.12846863) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0781773) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(1.0833441) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(2.892911) q[2];
sx q[2];
rz(-0.20771019) q[2];
sx q[2];
rz(-0.50616979) q[2];
rz(-0.32947551) q[3];
sx q[3];
rz(-1.4737045) q[3];
sx q[3];
rz(1.2274689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
