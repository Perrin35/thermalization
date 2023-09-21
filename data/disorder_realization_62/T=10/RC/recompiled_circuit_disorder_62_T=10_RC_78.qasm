OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463319) q[0];
sx q[0];
rz(-1.594461) q[0];
sx q[0];
rz(2.6803826) q[0];
rz(-pi) q[1];
rz(0.77779777) q[2];
sx q[2];
rz(-1.3133089) q[2];
sx q[2];
rz(-2.4553026) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9259778) q[1];
sx q[1];
rz(-1.3091334) q[1];
sx q[1];
rz(-1.91933) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1575559) q[3];
sx q[3];
rz(-0.26502702) q[3];
sx q[3];
rz(2.7812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87876451) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(-3.0905241) q[2];
rz(0.55705327) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(-0.59659514) q[0];
rz(-2.3157628) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.9155496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81235028) q[0];
sx q[0];
rz(-1.7698405) q[0];
sx q[0];
rz(-0.45254405) q[0];
rz(-pi) q[1];
rz(-2.366757) q[2];
sx q[2];
rz(-2.2171387) q[2];
sx q[2];
rz(-1.6006084) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1184428) q[1];
sx q[1];
rz(-0.21986248) q[1];
sx q[1];
rz(-1.5688194) q[1];
rz(-pi) q[2];
rz(-0.99818228) q[3];
sx q[3];
rz(-1.5628353) q[3];
sx q[3];
rz(-1.5919459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3423959) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-2.8105695) q[2];
rz(-2.3349169) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(1.249041) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(2.1121315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2217076) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(0.94421454) q[0];
rz(-pi) q[1];
rz(-2.4918409) q[2];
sx q[2];
rz(-1.3403112) q[2];
sx q[2];
rz(2.8284555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6368235) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(-1.7117281) q[1];
rz(-pi) q[2];
rz(1.449552) q[3];
sx q[3];
rz(-2.3070934) q[3];
sx q[3];
rz(2.8672225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.133698) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(-1.7525904) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(1.5456276) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(-1.625659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8221178) q[0];
sx q[0];
rz(-1.428077) q[0];
sx q[0];
rz(-0.6832173) q[0];
x q[1];
rz(0.085725351) q[2];
sx q[2];
rz(-2.1282196) q[2];
sx q[2];
rz(-2.2103708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8856636) q[1];
sx q[1];
rz(-0.71223488) q[1];
sx q[1];
rz(0.42339143) q[1];
rz(0.58873119) q[3];
sx q[3];
rz(-1.3874467) q[3];
sx q[3];
rz(-0.4962894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76379124) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.8466922) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(-1.5198583) q[0];
rz(-0.63201085) q[1];
sx q[1];
rz(-0.72223392) q[1];
sx q[1];
rz(2.246726) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3940112) q[0];
sx q[0];
rz(-2.3590439) q[0];
sx q[0];
rz(-2.5400019) q[0];
rz(-2.52389) q[2];
sx q[2];
rz(-0.82759826) q[2];
sx q[2];
rz(-0.4862116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.85256165) q[1];
sx q[1];
rz(-2.3173996) q[1];
sx q[1];
rz(2.186071) q[1];
rz(-pi) q[2];
rz(0.026167913) q[3];
sx q[3];
rz(-2.8550365) q[3];
sx q[3];
rz(-2.9122796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.616509) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.3775795) q[0];
rz(-2.6691943) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(0.46498743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156292) q[0];
sx q[0];
rz(-1.1383801) q[0];
sx q[0];
rz(-0.36884357) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9372413) q[2];
sx q[2];
rz(-2.7933279) q[2];
sx q[2];
rz(-1.5262926) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2496693) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(0.93530099) q[1];
rz(-2.898071) q[3];
sx q[3];
rz(-2.4176819) q[3];
sx q[3];
rz(0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.650699) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(0.4450376) q[2];
rz(2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0141107) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(1.7215464) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(-2.9856317) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0359081) q[0];
sx q[0];
rz(-1.7419635) q[0];
sx q[0];
rz(2.8842584) q[0];
rz(1.1440802) q[2];
sx q[2];
rz(-2.0442171) q[2];
sx q[2];
rz(-1.13525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1289039) q[1];
sx q[1];
rz(-1.3959937) q[1];
sx q[1];
rz(0.054794475) q[1];
rz(-1.0714674) q[3];
sx q[3];
rz(-1.0896177) q[3];
sx q[3];
rz(2.1035946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.069313958) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(-2.8159451) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9496562) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(-0.33777133) q[0];
rz(1.0900963) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(-2.8930194) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45390689) q[0];
sx q[0];
rz(-1.0605863) q[0];
sx q[0];
rz(-1.8574255) q[0];
rz(-0.82069223) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(-1.5936268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12411815) q[1];
sx q[1];
rz(-2.0548471) q[1];
sx q[1];
rz(-0.89213051) q[1];
rz(-pi) q[2];
x q[2];
rz(1.400984) q[3];
sx q[3];
rz(-1.4949993) q[3];
sx q[3];
rz(-2.4236987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(-2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(2.8588262) q[0];
rz(-2.4400318) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(-1.3185906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628172) q[0];
sx q[0];
rz(-1.6512198) q[0];
sx q[0];
rz(-1.1286939) q[0];
rz(-pi) q[1];
rz(2.1815427) q[2];
sx q[2];
rz(-1.0078444) q[2];
sx q[2];
rz(-0.27573953) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56837592) q[1];
sx q[1];
rz(-1.7163367) q[1];
sx q[1];
rz(-0.88370609) q[1];
rz(-pi) q[2];
rz(2.3913279) q[3];
sx q[3];
rz(-1.2816458) q[3];
sx q[3];
rz(2.9671448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-0.39917699) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(-1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76319641) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.8803966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1051837) q[0];
sx q[0];
rz(-1.1622218) q[0];
sx q[0];
rz(-1.4559792) q[0];
rz(-pi) q[1];
rz(3.10896) q[2];
sx q[2];
rz(-2.543078) q[2];
sx q[2];
rz(-2.0805217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2640472) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(1.643814) q[1];
x q[2];
rz(1.7301171) q[3];
sx q[3];
rz(-0.89309249) q[3];
sx q[3];
rz(-2.2417828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5579055) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(1.9231046) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(2.5789554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8284843) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(-2.5333511) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(-3.1153395) q[2];
sx q[2];
rz(-2.142341) q[2];
sx q[2];
rz(-0.061948902) q[2];
rz(2.9363587) q[3];
sx q[3];
rz(-0.60094613) q[3];
sx q[3];
rz(-0.080106674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];