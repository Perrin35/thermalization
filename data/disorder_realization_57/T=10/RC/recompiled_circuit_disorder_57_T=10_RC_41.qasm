OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(2.6101987) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(1.3637435) q[1];
sx q[1];
rz(10.627828) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8578313) q[0];
sx q[0];
rz(-0.28781048) q[0];
sx q[0];
rz(-2.3057111) q[0];
rz(-pi) q[1];
rz(1.3538023) q[2];
sx q[2];
rz(-1.3218244) q[2];
sx q[2];
rz(-2.4169902) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.10567084) q[1];
sx q[1];
rz(-2.8997313) q[1];
sx q[1];
rz(0.24005228) q[1];
rz(-0.12227998) q[3];
sx q[3];
rz(-2.8465726) q[3];
sx q[3];
rz(2.2515841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27665859) q[2];
sx q[2];
rz(-2.0099137) q[2];
sx q[2];
rz(-2.9425088) q[2];
rz(-1.52786) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(-0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99877015) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(0.81940991) q[0];
rz(2.8557414) q[1];
sx q[1];
rz(-0.91393036) q[1];
sx q[1];
rz(-1.8751289) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3730777) q[0];
sx q[0];
rz(-0.74155945) q[0];
sx q[0];
rz(-2.5624647) q[0];
rz(-pi) q[1];
rz(1.9752713) q[2];
sx q[2];
rz(-1.7911439) q[2];
sx q[2];
rz(1.599556) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.61050615) q[1];
sx q[1];
rz(-0.55263457) q[1];
sx q[1];
rz(2.3081739) q[1];
rz(-pi) q[2];
rz(-2.0224138) q[3];
sx q[3];
rz(-2.0997117) q[3];
sx q[3];
rz(3.0144453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1566029) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(1.161969) q[2];
rz(3.0544288) q[3];
sx q[3];
rz(-1.6379084) q[3];
sx q[3];
rz(-3.0676837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53482985) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(1.3820232) q[1];
sx q[1];
rz(-1.8654114) q[1];
sx q[1];
rz(-2.4940925) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377797) q[0];
sx q[0];
rz(-2.505629) q[0];
sx q[0];
rz(-3.0657835) q[0];
rz(-pi) q[1];
rz(2.1157672) q[2];
sx q[2];
rz(-0.93169824) q[2];
sx q[2];
rz(0.35312032) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4609047) q[1];
sx q[1];
rz(-2.540743) q[1];
sx q[1];
rz(-2.2520573) q[1];
rz(-pi) q[2];
rz(-2.2494499) q[3];
sx q[3];
rz(-0.8222848) q[3];
sx q[3];
rz(-3.126614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2993762) q[2];
sx q[2];
rz(-0.78263485) q[2];
sx q[2];
rz(-1.5608609) q[2];
rz(-0.98172274) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2340387) q[0];
sx q[0];
rz(-0.83775318) q[0];
sx q[0];
rz(-0.87189829) q[0];
rz(-2.7543228) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(-0.29104582) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1211348) q[0];
sx q[0];
rz(-0.99636787) q[0];
sx q[0];
rz(0.27340425) q[0];
rz(-pi) q[1];
rz(-0.83424763) q[2];
sx q[2];
rz(-0.4220037) q[2];
sx q[2];
rz(-3.1230502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2070771) q[1];
sx q[1];
rz(-1.0672489) q[1];
sx q[1];
rz(0.7657004) q[1];
rz(-2.980551) q[3];
sx q[3];
rz(-1.0003918) q[3];
sx q[3];
rz(0.17810861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7724093) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(2.5975361) q[2];
rz(0.50928515) q[3];
sx q[3];
rz(-2.0777168) q[3];
sx q[3];
rz(2.2756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0020224) q[0];
sx q[0];
rz(-2.1152571) q[0];
sx q[0];
rz(-0.64506662) q[0];
rz(0.6257261) q[1];
sx q[1];
rz(-1.2236243) q[1];
sx q[1];
rz(-0.63017875) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038534855) q[0];
sx q[0];
rz(-1.4487421) q[0];
sx q[0];
rz(-2.5621668) q[0];
rz(-1.9899881) q[2];
sx q[2];
rz(-0.79608166) q[2];
sx q[2];
rz(0.58005709) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2423357) q[1];
sx q[1];
rz(-1.7944272) q[1];
sx q[1];
rz(2.1171338) q[1];
rz(-2.9003521) q[3];
sx q[3];
rz(-0.57096982) q[3];
sx q[3];
rz(-2.8591448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3877635) q[2];
sx q[2];
rz(-1.8176273) q[2];
sx q[2];
rz(-1.7774263) q[2];
rz(-1.8917313) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(-2.2201339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0548148) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(0.71682799) q[0];
rz(-2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(0.81370083) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98052927) q[0];
sx q[0];
rz(-1.6335532) q[0];
sx q[0];
rz(1.784523) q[0];
rz(-pi) q[1];
rz(0.72197638) q[2];
sx q[2];
rz(-1.1140031) q[2];
sx q[2];
rz(-1.7320088) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.54436436) q[1];
sx q[1];
rz(-1.2372984) q[1];
sx q[1];
rz(2.3436105) q[1];
x q[2];
rz(-1.9244476) q[3];
sx q[3];
rz(-1.2815223) q[3];
sx q[3];
rz(-1.7464964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8852691) q[2];
sx q[2];
rz(-0.43748125) q[2];
sx q[2];
rz(1.1876594) q[2];
rz(-2.2096283) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(-0.37117547) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9880144) q[0];
sx q[0];
rz(-2.1126641) q[0];
sx q[0];
rz(2.4555092) q[0];
rz(-1.6185121) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(2.4783321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2880741) q[0];
sx q[0];
rz(-0.92068866) q[0];
sx q[0];
rz(0.030036146) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63753193) q[2];
sx q[2];
rz(-1.7469284) q[2];
sx q[2];
rz(-2.9498364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38695947) q[1];
sx q[1];
rz(-0.52801758) q[1];
sx q[1];
rz(1.1623043) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5797964) q[3];
sx q[3];
rz(-1.1601163) q[3];
sx q[3];
rz(-0.79515776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.474581) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(-3.1265124) q[2];
rz(1.0117426) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(2.570178) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886154) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(-3.0342039) q[0];
rz(-0.30474162) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(-2.6838141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37516884) q[0];
sx q[0];
rz(-0.94978118) q[0];
sx q[0];
rz(0.39874886) q[0];
rz(2.1019943) q[2];
sx q[2];
rz(-0.89745144) q[2];
sx q[2];
rz(1.4199867) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.01955186) q[1];
sx q[1];
rz(-1.8035839) q[1];
sx q[1];
rz(-1.9368534) q[1];
rz(-pi) q[2];
rz(-1.9321953) q[3];
sx q[3];
rz(-2.0933588) q[3];
sx q[3];
rz(2.2613329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7765939) q[2];
sx q[2];
rz(-1.4564617) q[2];
sx q[2];
rz(-1.7101074) q[2];
rz(1.4252023) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(-1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15437056) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(0.95348683) q[0];
rz(-1.3257239) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(-2.5568533) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1138211) q[0];
sx q[0];
rz(-0.50243938) q[0];
sx q[0];
rz(-0.2343982) q[0];
rz(-pi) q[1];
rz(-1.936475) q[2];
sx q[2];
rz(-1.056991) q[2];
sx q[2];
rz(-0.50525451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19496275) q[1];
sx q[1];
rz(-1.1319036) q[1];
sx q[1];
rz(1.4942601) q[1];
rz(-2.0546012) q[3];
sx q[3];
rz(-0.30011794) q[3];
sx q[3];
rz(2.8407612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8087625) q[2];
sx q[2];
rz(-2.3213883) q[2];
sx q[2];
rz(-0.27627036) q[2];
rz(0.55109465) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0693414) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(1.9845225) q[1];
sx q[1];
rz(-1.3815222) q[1];
sx q[1];
rz(1.5225333) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9690543) q[0];
sx q[0];
rz(-1.8803055) q[0];
sx q[0];
rz(-0.68185101) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40752878) q[2];
sx q[2];
rz(-0.73020836) q[2];
sx q[2];
rz(2.4804295) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9138227) q[1];
sx q[1];
rz(-1.4064944) q[1];
sx q[1];
rz(-0.63197244) q[1];
rz(2.7195815) q[3];
sx q[3];
rz(-1.2000298) q[3];
sx q[3];
rz(-1.5164204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9860501) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(-0.84214169) q[2];
rz(1.6048253) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-0.026528927) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2355272) q[0];
sx q[0];
rz(-1.9530095) q[0];
sx q[0];
rz(-0.45146913) q[0];
rz(0.72262598) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(0.026635219) q[2];
sx q[2];
rz(-1.3616818) q[2];
sx q[2];
rz(1.5875265) q[2];
rz(2.261653) q[3];
sx q[3];
rz(-1.6321245) q[3];
sx q[3];
rz(1.5499916) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
