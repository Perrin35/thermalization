OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(-0.92209446) q[0];
sx q[0];
rz(2.3715012) q[0];
rz(-2.4251179) q[1];
sx q[1];
rz(-0.94243503) q[1];
sx q[1];
rz(2.4997349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7866369) q[0];
sx q[0];
rz(-1.4572976) q[0];
sx q[0];
rz(-2.9391975) q[0];
rz(-0.0014512295) q[2];
sx q[2];
rz(-1.5696758) q[2];
sx q[2];
rz(-1.6482774) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.227046) q[1];
sx q[1];
rz(-1.7265336) q[1];
sx q[1];
rz(-2.7776633) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58310469) q[3];
sx q[3];
rz(-2.2197086) q[3];
sx q[3];
rz(0.16498868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98051071) q[2];
sx q[2];
rz(-1.0970205) q[2];
sx q[2];
rz(2.3316627) q[2];
rz(3.1071281) q[3];
sx q[3];
rz(-0.66027111) q[3];
sx q[3];
rz(0.1035498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63475364) q[0];
sx q[0];
rz(-2.6320808) q[0];
sx q[0];
rz(-0.65938812) q[0];
rz(-1.7022645) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(2.6606681) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.07671) q[0];
sx q[0];
rz(-2.0079029) q[0];
sx q[0];
rz(-0.22910288) q[0];
rz(2.2873053) q[2];
sx q[2];
rz(-0.97161869) q[2];
sx q[2];
rz(2.6024659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7742851) q[1];
sx q[1];
rz(-1.7789504) q[1];
sx q[1];
rz(-0.03117301) q[1];
x q[2];
rz(2.0388076) q[3];
sx q[3];
rz(-2.0986522) q[3];
sx q[3];
rz(-2.6754987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3769569) q[2];
sx q[2];
rz(-0.83676338) q[2];
sx q[2];
rz(2.5118206) q[2];
rz(-1.9624286) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(2.0298957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4513007) q[0];
sx q[0];
rz(-3.1307104) q[0];
sx q[0];
rz(-2.1731398) q[0];
rz(-0.14006607) q[1];
sx q[1];
rz(-1.7882971) q[1];
sx q[1];
rz(0.5828988) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4752858) q[0];
sx q[0];
rz(-0.17284849) q[0];
sx q[0];
rz(-0.26474196) q[0];
rz(-pi) q[1];
rz(-0.83816041) q[2];
sx q[2];
rz(-1.3930905) q[2];
sx q[2];
rz(3.1237683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0576833) q[1];
sx q[1];
rz(-2.1321802) q[1];
sx q[1];
rz(2.3777005) q[1];
rz(-pi) q[2];
rz(-1.8898592) q[3];
sx q[3];
rz(-1.3033861) q[3];
sx q[3];
rz(1.2280994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6110903) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(0.21406847) q[2];
rz(-0.30238447) q[3];
sx q[3];
rz(-0.74367911) q[3];
sx q[3];
rz(-2.7157057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532303) q[0];
sx q[0];
rz(-1.0091877) q[0];
sx q[0];
rz(-2.0971712) q[0];
rz(-2.2638679) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(-2.9728319) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2917738) q[0];
sx q[0];
rz(-0.55233228) q[0];
sx q[0];
rz(-2.0122347) q[0];
rz(-0.20511583) q[2];
sx q[2];
rz(-1.2502708) q[2];
sx q[2];
rz(-1.6728624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8873621) q[1];
sx q[1];
rz(-1.7866644) q[1];
sx q[1];
rz(2.5309023) q[1];
rz(-0.041958001) q[3];
sx q[3];
rz(-0.86554147) q[3];
sx q[3];
rz(-2.369997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7793444) q[2];
sx q[2];
rz(-1.2835953) q[2];
sx q[2];
rz(1.0603504) q[2];
rz(-0.29252163) q[3];
sx q[3];
rz(-2.4157603) q[3];
sx q[3];
rz(3.0574851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58920687) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(0.86828434) q[0];
rz(-2.5153416) q[1];
sx q[1];
rz(-1.7528844) q[1];
sx q[1];
rz(-1.7832322) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043644301) q[0];
sx q[0];
rz(-1.4820423) q[0];
sx q[0];
rz(2.8525145) q[0];
rz(0.50164671) q[2];
sx q[2];
rz(-2.7124462) q[2];
sx q[2];
rz(-2.8535064) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2430181) q[1];
sx q[1];
rz(-2.822091) q[1];
sx q[1];
rz(2.363149) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63774469) q[3];
sx q[3];
rz(-1.2079835) q[3];
sx q[3];
rz(-2.2408298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2061578) q[2];
sx q[2];
rz(-2.3249966) q[2];
sx q[2];
rz(-2.6182776) q[2];
rz(-0.37007904) q[3];
sx q[3];
rz(-0.78032929) q[3];
sx q[3];
rz(1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0321781) q[0];
sx q[0];
rz(-0.21571708) q[0];
sx q[0];
rz(-0.35414645) q[0];
rz(0.94193637) q[1];
sx q[1];
rz(-1.6862005) q[1];
sx q[1];
rz(1.3166434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8982074) q[0];
sx q[0];
rz(-1.4665717) q[0];
sx q[0];
rz(-0.51527713) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5923827) q[2];
sx q[2];
rz(-1.5592279) q[2];
sx q[2];
rz(2.9988097) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9256014) q[1];
sx q[1];
rz(-2.9570438) q[1];
sx q[1];
rz(-0.33949236) q[1];
x q[2];
rz(-1.8527669) q[3];
sx q[3];
rz(-0.94579711) q[3];
sx q[3];
rz(-0.19240141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3256623) q[2];
sx q[2];
rz(-2.7553835) q[2];
sx q[2];
rz(2.8032934) q[2];
rz(0.47652388) q[3];
sx q[3];
rz(-0.75348133) q[3];
sx q[3];
rz(-0.34059718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4395831) q[0];
sx q[0];
rz(-1.5583353) q[0];
sx q[0];
rz(-0.53771341) q[0];
rz(1.4078377) q[1];
sx q[1];
rz(-2.6815963) q[1];
sx q[1];
rz(0.62526155) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7571676) q[0];
sx q[0];
rz(-1.0119857) q[0];
sx q[0];
rz(1.1435056) q[0];
rz(-pi) q[1];
rz(0.12282108) q[2];
sx q[2];
rz(-2.5993957) q[2];
sx q[2];
rz(1.922883) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1792153) q[1];
sx q[1];
rz(-1.3562225) q[1];
sx q[1];
rz(-2.5291165) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2042562) q[3];
sx q[3];
rz(-0.82895422) q[3];
sx q[3];
rz(2.3693565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9739146) q[2];
sx q[2];
rz(-0.46678552) q[2];
sx q[2];
rz(-1.3803049) q[2];
rz(2.7050833) q[3];
sx q[3];
rz(-2.1197539) q[3];
sx q[3];
rz(-0.71353394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1768782) q[0];
sx q[0];
rz(-2.5202993) q[0];
sx q[0];
rz(2.6253413) q[0];
rz(2.3618354) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(2.0947184) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31016997) q[0];
sx q[0];
rz(-2.8084233) q[0];
sx q[0];
rz(-2.6875671) q[0];
rz(-0.70430906) q[2];
sx q[2];
rz(-2.4302135) q[2];
sx q[2];
rz(0.34397438) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.88491625) q[1];
sx q[1];
rz(-3.0312928) q[1];
sx q[1];
rz(-2.1259456) q[1];
x q[2];
rz(-0.67892142) q[3];
sx q[3];
rz(-0.58829868) q[3];
sx q[3];
rz(3.1117518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.20389916) q[2];
sx q[2];
rz(-0.1940618) q[2];
sx q[2];
rz(-0.95721179) q[2];
rz(-0.29414487) q[3];
sx q[3];
rz(-1.1295986) q[3];
sx q[3];
rz(2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1421563) q[0];
sx q[0];
rz(-2.9111828) q[0];
sx q[0];
rz(0.18705046) q[0];
rz(-0.43141463) q[1];
sx q[1];
rz(-0.43771935) q[1];
sx q[1];
rz(1.4923219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6607912) q[0];
sx q[0];
rz(-1.4993164) q[0];
sx q[0];
rz(3.0707703) q[0];
rz(-0.19047462) q[2];
sx q[2];
rz(-1.8996793) q[2];
sx q[2];
rz(2.637907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2577511) q[1];
sx q[1];
rz(-1.7353936) q[1];
sx q[1];
rz(1.1562892) q[1];
rz(2.1077376) q[3];
sx q[3];
rz(-2.9599422) q[3];
sx q[3];
rz(0.27508116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46818587) q[2];
sx q[2];
rz(-0.91172051) q[2];
sx q[2];
rz(2.1569596) q[2];
rz(0.51472384) q[3];
sx q[3];
rz(-0.52557164) q[3];
sx q[3];
rz(-1.908186) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24432261) q[0];
sx q[0];
rz(-1.4800625) q[0];
sx q[0];
rz(-2.3829714) q[0];
rz(-1.1933391) q[1];
sx q[1];
rz(-1.9494282) q[1];
sx q[1];
rz(1.6903445) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0139931) q[0];
sx q[0];
rz(-1.7993449) q[0];
sx q[0];
rz(-0.20968135) q[0];
rz(-pi) q[1];
rz(1.7733943) q[2];
sx q[2];
rz(-0.2066863) q[2];
sx q[2];
rz(-2.8967173) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32827863) q[1];
sx q[1];
rz(-2.5254888) q[1];
sx q[1];
rz(-1.5244984) q[1];
rz(-pi) q[2];
x q[2];
rz(1.162446) q[3];
sx q[3];
rz(-1.2733885) q[3];
sx q[3];
rz(1.4434467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54185581) q[2];
sx q[2];
rz(-2.2208322) q[2];
sx q[2];
rz(2.3403781) q[2];
rz(0.93252212) q[3];
sx q[3];
rz(-1.9263809) q[3];
sx q[3];
rz(-2.7264989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5781317) q[0];
sx q[0];
rz(-1.3787855) q[0];
sx q[0];
rz(-0.61510573) q[0];
rz(-0.32147944) q[1];
sx q[1];
rz(-0.97578661) q[1];
sx q[1];
rz(-1.6937561) q[1];
rz(0.16531113) q[2];
sx q[2];
rz(-1.4480235) q[2];
sx q[2];
rz(1.5245246) q[2];
rz(2.2533909) q[3];
sx q[3];
rz(-2.785688) q[3];
sx q[3];
rz(1.9733081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
