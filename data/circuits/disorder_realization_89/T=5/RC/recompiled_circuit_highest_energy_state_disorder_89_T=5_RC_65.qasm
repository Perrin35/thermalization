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
rz(1.8472449) q[0];
sx q[0];
rz(-1.6142774) q[0];
sx q[0];
rz(2.9543028) q[0];
rz(-1.3297431) q[1];
sx q[1];
rz(-2.5951374) q[1];
sx q[1];
rz(1.3475013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2696323) q[0];
sx q[0];
rz(-0.56817164) q[0];
sx q[0];
rz(-1.7164036) q[0];
x q[1];
rz(2.7653469) q[2];
sx q[2];
rz(-1.0828138) q[2];
sx q[2];
rz(-0.03586344) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.076278585) q[1];
sx q[1];
rz(-2.4597628) q[1];
sx q[1];
rz(-1.5449468) q[1];
rz(-pi) q[2];
rz(2.4745117) q[3];
sx q[3];
rz(-0.32235185) q[3];
sx q[3];
rz(-1.9867112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61907855) q[2];
sx q[2];
rz(-1.1828902) q[2];
sx q[2];
rz(-2.8600233) q[2];
rz(1.36739) q[3];
sx q[3];
rz(-1.2308246) q[3];
sx q[3];
rz(-0.27845964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.53681579) q[0];
sx q[0];
rz(-0.29412687) q[0];
sx q[0];
rz(1.7915223) q[0];
rz(0.049909441) q[1];
sx q[1];
rz(-0.3388181) q[1];
sx q[1];
rz(1.2443589) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45634633) q[0];
sx q[0];
rz(-2.7875627) q[0];
sx q[0];
rz(0.1683321) q[0];
rz(-pi) q[1];
rz(-2.9551465) q[2];
sx q[2];
rz(-2.1121587) q[2];
sx q[2];
rz(0.10554927) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0383953) q[1];
sx q[1];
rz(-1.8382676) q[1];
sx q[1];
rz(0.20305509) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4723827) q[3];
sx q[3];
rz(-2.4871181) q[3];
sx q[3];
rz(1.7184642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0758948) q[2];
sx q[2];
rz(-1.9588797) q[2];
sx q[2];
rz(-3.0720224) q[2];
rz(2.4075497) q[3];
sx q[3];
rz(-2.3591122) q[3];
sx q[3];
rz(2.5231584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0760536) q[0];
sx q[0];
rz(-2.9422308) q[0];
sx q[0];
rz(-2.8850436) q[0];
rz(1.8007295) q[1];
sx q[1];
rz(-0.95756617) q[1];
sx q[1];
rz(1.0440913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1397823) q[0];
sx q[0];
rz(-0.35813791) q[0];
sx q[0];
rz(2.0689808) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3816517) q[2];
sx q[2];
rz(-0.30549288) q[2];
sx q[2];
rz(-1.0585275) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0716182) q[1];
sx q[1];
rz(-0.41980068) q[1];
sx q[1];
rz(-0.75666766) q[1];
x q[2];
rz(-0.14117853) q[3];
sx q[3];
rz(-2.1394844) q[3];
sx q[3];
rz(1.6077667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9112245) q[2];
sx q[2];
rz(-2.0679097) q[2];
sx q[2];
rz(2.4252841) q[2];
rz(2.6835942) q[3];
sx q[3];
rz(-1.674998) q[3];
sx q[3];
rz(-0.23076335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79391795) q[0];
sx q[0];
rz(-2.1383998) q[0];
sx q[0];
rz(-2.4804261) q[0];
rz(-1.2541153) q[1];
sx q[1];
rz(-1.5690119) q[1];
sx q[1];
rz(-1.5987781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1791627) q[0];
sx q[0];
rz(-1.822285) q[0];
sx q[0];
rz(-0.18119399) q[0];
x q[1];
rz(-2.6933367) q[2];
sx q[2];
rz(-2.087062) q[2];
sx q[2];
rz(-0.22656952) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99847163) q[1];
sx q[1];
rz(-2.4354092) q[1];
sx q[1];
rz(0.097340214) q[1];
x q[2];
rz(-2.1041252) q[3];
sx q[3];
rz(-1.2308443) q[3];
sx q[3];
rz(-1.8890017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8864112) q[2];
sx q[2];
rz(-0.35021439) q[2];
sx q[2];
rz(1.1345471) q[2];
rz(2.5875731) q[3];
sx q[3];
rz(-1.3628682) q[3];
sx q[3];
rz(-2.5567283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2968813) q[0];
sx q[0];
rz(-0.0021332707) q[0];
sx q[0];
rz(2.010349) q[0];
rz(3.0834815) q[1];
sx q[1];
rz(-1.2429712) q[1];
sx q[1];
rz(1.4505454) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03019985) q[0];
sx q[0];
rz(-0.95472017) q[0];
sx q[0];
rz(0.48171039) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8314231) q[2];
sx q[2];
rz(-0.7855987) q[2];
sx q[2];
rz(0.73332649) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1102396) q[1];
sx q[1];
rz(-1.3340557) q[1];
sx q[1];
rz(2.8048022) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2048619) q[3];
sx q[3];
rz(-1.2311683) q[3];
sx q[3];
rz(-0.32138164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45945534) q[2];
sx q[2];
rz(-2.5294999) q[2];
sx q[2];
rz(-1.7876145) q[2];
rz(0.55245429) q[3];
sx q[3];
rz(-2.2508636) q[3];
sx q[3];
rz(2.6433105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.246493) q[0];
sx q[0];
rz(-2.1857078) q[0];
sx q[0];
rz(1.0981052) q[0];
rz(1.8662628) q[1];
sx q[1];
rz(-1.2397436) q[1];
sx q[1];
rz(1.0414418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7324243) q[0];
sx q[0];
rz(-2.083626) q[0];
sx q[0];
rz(0.70968117) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8813754) q[2];
sx q[2];
rz(-0.12311664) q[2];
sx q[2];
rz(-2.3359063) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84717709) q[1];
sx q[1];
rz(-1.4694459) q[1];
sx q[1];
rz(2.0854998) q[1];
rz(-0.28686974) q[3];
sx q[3];
rz(-0.86599819) q[3];
sx q[3];
rz(-1.4554086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2991403) q[2];
sx q[2];
rz(-1.9879397) q[2];
sx q[2];
rz(3.0843206) q[2];
rz(0.20036571) q[3];
sx q[3];
rz(-0.55557957) q[3];
sx q[3];
rz(0.19181767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361236) q[0];
sx q[0];
rz(-2.496688) q[0];
sx q[0];
rz(-2.7591925) q[0];
rz(2.1740225) q[1];
sx q[1];
rz(-2.7303374) q[1];
sx q[1];
rz(1.8866906) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8604831) q[0];
sx q[0];
rz(-2.5279928) q[0];
sx q[0];
rz(-3.0439418) q[0];
x q[1];
rz(-0.91302769) q[2];
sx q[2];
rz(-2.0951197) q[2];
sx q[2];
rz(-1.9425478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69051504) q[1];
sx q[1];
rz(-1.9743063) q[1];
sx q[1];
rz(-2.7355641) q[1];
rz(-1.8462049) q[3];
sx q[3];
rz(-0.81116889) q[3];
sx q[3];
rz(-1.0397183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3249698) q[2];
sx q[2];
rz(-0.35494706) q[2];
sx q[2];
rz(-0.85186446) q[2];
rz(2.8999117) q[3];
sx q[3];
rz(-1.207374) q[3];
sx q[3];
rz(-2.2309301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8082751) q[0];
sx q[0];
rz(-2.9589544) q[0];
sx q[0];
rz(-2.0727378) q[0];
rz(-1.3775657) q[1];
sx q[1];
rz(-2.865538) q[1];
sx q[1];
rz(-0.68881234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6816872) q[0];
sx q[0];
rz(-1.9294159) q[0];
sx q[0];
rz(-1.9214517) q[0];
rz(-pi) q[1];
x q[1];
rz(2.942932) q[2];
sx q[2];
rz(-2.5294494) q[2];
sx q[2];
rz(2.8191301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.090052) q[1];
sx q[1];
rz(-1.2512815) q[1];
sx q[1];
rz(-0.40488147) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52940552) q[3];
sx q[3];
rz(-1.7200618) q[3];
sx q[3];
rz(-0.68219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0546725) q[2];
sx q[2];
rz(-0.72267795) q[2];
sx q[2];
rz(1.5642536) q[2];
rz(-0.65586048) q[3];
sx q[3];
rz(-1.4411996) q[3];
sx q[3];
rz(-0.76914579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5250788) q[0];
sx q[0];
rz(-1.0322796) q[0];
sx q[0];
rz(2.9922564) q[0];
rz(-2.0556045) q[1];
sx q[1];
rz(-2.1859152) q[1];
sx q[1];
rz(2.65926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457561) q[0];
sx q[0];
rz(-2.0424423) q[0];
sx q[0];
rz(-0.095998569) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6012548) q[2];
sx q[2];
rz(-0.23269146) q[2];
sx q[2];
rz(0.57664492) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0808951) q[1];
sx q[1];
rz(-1.2978862) q[1];
sx q[1];
rz(-0.70571438) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6685358) q[3];
sx q[3];
rz(-1.9888788) q[3];
sx q[3];
rz(3.0168777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.09717354) q[2];
sx q[2];
rz(-0.83690518) q[2];
sx q[2];
rz(-2.6824299) q[2];
rz(-0.098371355) q[3];
sx q[3];
rz(-1.8561074) q[3];
sx q[3];
rz(1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2761053) q[0];
sx q[0];
rz(-1.2609755) q[0];
sx q[0];
rz(-0.14044811) q[0];
rz(1.2175995) q[1];
sx q[1];
rz(-1.1414889) q[1];
sx q[1];
rz(-1.5906895) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63839943) q[0];
sx q[0];
rz(-2.2941883) q[0];
sx q[0];
rz(1.6795078) q[0];
x q[1];
rz(1.7459007) q[2];
sx q[2];
rz(-1.7410472) q[2];
sx q[2];
rz(-2.424602) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5441011) q[1];
sx q[1];
rz(-0.95817503) q[1];
sx q[1];
rz(-1.2493916) q[1];
rz(2.2416653) q[3];
sx q[3];
rz(-1.560671) q[3];
sx q[3];
rz(-0.72521082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5156252) q[2];
sx q[2];
rz(-1.0270065) q[2];
sx q[2];
rz(-1.1930126) q[2];
rz(2.6595645) q[3];
sx q[3];
rz(-2.7218282) q[3];
sx q[3];
rz(2.1217864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0057509) q[0];
sx q[0];
rz(-1.0938205) q[0];
sx q[0];
rz(-0.83458207) q[0];
rz(2.3781378) q[1];
sx q[1];
rz(-1.8370942) q[1];
sx q[1];
rz(-2.7948517) q[1];
rz(1.4540666) q[2];
sx q[2];
rz(-1.5925831) q[2];
sx q[2];
rz(1.5769373) q[2];
rz(-2.5478195) q[3];
sx q[3];
rz(-1.9304921) q[3];
sx q[3];
rz(-1.0747128) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
