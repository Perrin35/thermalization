OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(0.93044257) q[0];
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(1.1073444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005172) q[0];
sx q[0];
rz(-2.374875) q[0];
sx q[0];
rz(0.47857743) q[0];
rz(-pi) q[1];
rz(2.6644457) q[2];
sx q[2];
rz(-2.2331182) q[2];
sx q[2];
rz(-1.9805816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90802898) q[1];
sx q[1];
rz(-1.4587914) q[1];
sx q[1];
rz(1.8506552) q[1];
rz(-pi) q[2];
rz(2.0290306) q[3];
sx q[3];
rz(-0.68800612) q[3];
sx q[3];
rz(-2.1306761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.063623108) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(-1.3280274) q[2];
rz(0.32087457) q[3];
sx q[3];
rz(-2.1556373) q[3];
sx q[3];
rz(-3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6593453) q[0];
sx q[0];
rz(-3.0292065) q[0];
sx q[0];
rz(2.2609718) q[0];
rz(-1.2940787) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(2.3243288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6200703) q[0];
sx q[0];
rz(-1.2261454) q[0];
sx q[0];
rz(0.14981139) q[0];
rz(-pi) q[1];
rz(1.9643289) q[2];
sx q[2];
rz(-0.52429188) q[2];
sx q[2];
rz(-2.2146068) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0751794) q[1];
sx q[1];
rz(-1.8261837) q[1];
sx q[1];
rz(0.65402072) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0975295) q[3];
sx q[3];
rz(-1.739199) q[3];
sx q[3];
rz(0.67697224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2237504) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(0.98637995) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(-1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.36689511) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(0.96631518) q[0];
rz(-2.9486588) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(-1.6945217) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23611785) q[0];
sx q[0];
rz(-0.21572278) q[0];
sx q[0];
rz(-3.059379) q[0];
rz(-pi) q[1];
rz(-1.2800531) q[2];
sx q[2];
rz(-1.650052) q[2];
sx q[2];
rz(0.62627625) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24070534) q[1];
sx q[1];
rz(-1.9449688) q[1];
sx q[1];
rz(0.64970533) q[1];
x q[2];
rz(-2.0855911) q[3];
sx q[3];
rz(-0.47009531) q[3];
sx q[3];
rz(-1.521829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43859279) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(1.1052216) q[2];
rz(-0.74622074) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(-2.1658649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352017) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.4659457) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-1.0703215) q[1];
sx q[1];
rz(-0.38898653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5327832) q[0];
sx q[0];
rz(-0.36768498) q[0];
sx q[0];
rz(-0.68433783) q[0];
x q[1];
rz(2.0573425) q[2];
sx q[2];
rz(-1.7440737) q[2];
sx q[2];
rz(1.372352) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3868689) q[1];
sx q[1];
rz(-1.9935973) q[1];
sx q[1];
rz(-2.5368284) q[1];
x q[2];
rz(1.7763406) q[3];
sx q[3];
rz(-2.3081429) q[3];
sx q[3];
rz(1.2780241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7230364) q[2];
sx q[2];
rz(-1.1636461) q[2];
sx q[2];
rz(0.45670613) q[2];
rz(1.6263973) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91963768) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(0.36002457) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(2.6470851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3213145) q[0];
sx q[0];
rz(-1.2914133) q[0];
sx q[0];
rz(1.4415635) q[0];
rz(-0.039426609) q[2];
sx q[2];
rz(-1.4354424) q[2];
sx q[2];
rz(0.53265041) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8620017) q[1];
sx q[1];
rz(-1.5108567) q[1];
sx q[1];
rz(-2.5111591) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6974929) q[3];
sx q[3];
rz(-1.5937623) q[3];
sx q[3];
rz(0.75957739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71022025) q[2];
sx q[2];
rz(-2.4857095) q[2];
sx q[2];
rz(2.8430856) q[2];
rz(-0.31202894) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(-0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9962149) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(-0.19590713) q[0];
rz(3.1205102) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(-1.235199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346749) q[0];
sx q[0];
rz(-1.2376681) q[0];
sx q[0];
rz(-1.8963277) q[0];
rz(-pi) q[1];
rz(1.2321129) q[2];
sx q[2];
rz(-2.2465696) q[2];
sx q[2];
rz(2.0896926) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1393309) q[1];
sx q[1];
rz(-1.1552703) q[1];
sx q[1];
rz(-2.5050487) q[1];
x q[2];
rz(-1.487021) q[3];
sx q[3];
rz(-1.8063074) q[3];
sx q[3];
rz(-1.868353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6482676) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(-0.43169272) q[2];
rz(1.4124983) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7504808) q[0];
sx q[0];
rz(-1.1688122) q[0];
sx q[0];
rz(-3.0294763) q[0];
rz(-0.21513367) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(-2.0281866) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33681413) q[0];
sx q[0];
rz(-1.830173) q[0];
sx q[0];
rz(0.39774261) q[0];
x q[1];
rz(2.8261975) q[2];
sx q[2];
rz(-1.3239667) q[2];
sx q[2];
rz(-2.18404) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.73247319) q[1];
sx q[1];
rz(-1.2670244) q[1];
sx q[1];
rz(-0.71883808) q[1];
rz(-pi) q[2];
rz(0.21344276) q[3];
sx q[3];
rz(-1.1323954) q[3];
sx q[3];
rz(-1.9813117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.001174288) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(0.12602885) q[2];
rz(2.0942988) q[3];
sx q[3];
rz(-1.8208561) q[3];
sx q[3];
rz(2.4333911) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66013181) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(-1.2449645) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(1.2089027) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395994) q[0];
sx q[0];
rz(-1.8327466) q[0];
sx q[0];
rz(0.082017935) q[0];
x q[1];
rz(0.65931321) q[2];
sx q[2];
rz(-0.94617832) q[2];
sx q[2];
rz(-2.5059932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0749082) q[1];
sx q[1];
rz(-1.2411331) q[1];
sx q[1];
rz(-2.3782905) q[1];
rz(-pi) q[2];
rz(-2.9221411) q[3];
sx q[3];
rz(-2.1900574) q[3];
sx q[3];
rz(1.0964364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(1.3195999) q[2];
rz(0.59213263) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(-0.035141703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-0.2329344) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(1.7804902) q[0];
rz(-1.9305485) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(2.7499054) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93145934) q[0];
sx q[0];
rz(-1.1362695) q[0];
sx q[0];
rz(-1.5772485) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3707156) q[2];
sx q[2];
rz(-2.0152425) q[2];
sx q[2];
rz(2.0796298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45410941) q[1];
sx q[1];
rz(-2.3309163) q[1];
sx q[1];
rz(-1.1078542) q[1];
rz(0.52569315) q[3];
sx q[3];
rz(-1.634053) q[3];
sx q[3];
rz(-0.69378187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52035511) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(1.8048145) q[2];
rz(-0.76198602) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(-2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8005463) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(-2.5706932) q[0];
rz(1.7123429) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(-2.9796519) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83308342) q[0];
sx q[0];
rz(-0.97531318) q[0];
sx q[0];
rz(2.7916629) q[0];
x q[1];
rz(-0.0091153092) q[2];
sx q[2];
rz(-1.7593804) q[2];
sx q[2];
rz(1.7878704) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0429749) q[1];
sx q[1];
rz(-2.9264755) q[1];
sx q[1];
rz(2.2153562) q[1];
rz(-pi) q[2];
rz(3.083769) q[3];
sx q[3];
rz(-2.6177546) q[3];
sx q[3];
rz(0.1334838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1841715) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(-1.4282248) q[2];
rz(-1.1994294) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(-1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7006871) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.7715001) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(2.2489108) q[2];
sx q[2];
rz(-0.18845367) q[2];
sx q[2];
rz(-1.0343196) q[2];
rz(0.45392848) q[3];
sx q[3];
rz(-0.47149999) q[3];
sx q[3];
rz(-0.75054689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];