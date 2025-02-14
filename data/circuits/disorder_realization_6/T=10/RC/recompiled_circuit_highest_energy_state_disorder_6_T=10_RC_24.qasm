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
rz(0.9917292) q[0];
sx q[0];
rz(-0.60757929) q[0];
sx q[0];
rz(0.99555957) q[0];
rz(-2.6140656) q[1];
sx q[1];
rz(-1.0310643) q[1];
sx q[1];
rz(2.6568025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6226396) q[0];
sx q[0];
rz(-2.6128809) q[0];
sx q[0];
rz(1.8334808) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44102168) q[2];
sx q[2];
rz(-0.82953757) q[2];
sx q[2];
rz(1.5462359) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0409702) q[1];
sx q[1];
rz(-0.61196152) q[1];
sx q[1];
rz(-1.4684567) q[1];
x q[2];
rz(1.6438585) q[3];
sx q[3];
rz(-1.8612075) q[3];
sx q[3];
rz(-0.452347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.63554865) q[2];
sx q[2];
rz(-0.94749331) q[2];
sx q[2];
rz(-0.49723899) q[2];
rz(2.5908568) q[3];
sx q[3];
rz(-0.68998718) q[3];
sx q[3];
rz(2.0042888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15744844) q[0];
sx q[0];
rz(-1.4009615) q[0];
sx q[0];
rz(-2.3781811) q[0];
rz(2.5570671) q[1];
sx q[1];
rz(-1.3817363) q[1];
sx q[1];
rz(-1.0113299) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9824163) q[0];
sx q[0];
rz(-1.455602) q[0];
sx q[0];
rz(-1.4388491) q[0];
rz(0.44707985) q[2];
sx q[2];
rz(-0.26105389) q[2];
sx q[2];
rz(2.3787468) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3608777) q[1];
sx q[1];
rz(-1.3643801) q[1];
sx q[1];
rz(1.92093) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2143651) q[3];
sx q[3];
rz(-0.6112186) q[3];
sx q[3];
rz(0.54655308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6611019) q[2];
sx q[2];
rz(-1.9384408) q[2];
sx q[2];
rz(-2.1413546) q[2];
rz(-2.5601322) q[3];
sx q[3];
rz(-0.73271078) q[3];
sx q[3];
rz(-1.2015517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8345399) q[0];
sx q[0];
rz(-2.5426799) q[0];
sx q[0];
rz(-0.55759984) q[0];
rz(-0.3071951) q[1];
sx q[1];
rz(-2.5418044) q[1];
sx q[1];
rz(-0.44816005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9014943) q[0];
sx q[0];
rz(-1.0567825) q[0];
sx q[0];
rz(-1.4185403) q[0];
rz(-pi) q[1];
rz(1.2101754) q[2];
sx q[2];
rz(-0.54979392) q[2];
sx q[2];
rz(0.87042337) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8413755) q[1];
sx q[1];
rz(-0.63134495) q[1];
sx q[1];
rz(-1.6168827) q[1];
x q[2];
rz(1.5968733) q[3];
sx q[3];
rz(-1.136565) q[3];
sx q[3];
rz(-0.48420478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8223411) q[2];
sx q[2];
rz(-2.4288869) q[2];
sx q[2];
rz(2.5437497) q[2];
rz(0.44761014) q[3];
sx q[3];
rz(-1.2667344) q[3];
sx q[3];
rz(-3.083057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7287801) q[0];
sx q[0];
rz(-0.030981177) q[0];
sx q[0];
rz(0.70984167) q[0];
rz(-2.1624883) q[1];
sx q[1];
rz(-0.85936463) q[1];
sx q[1];
rz(-1.3213347) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0036572) q[0];
sx q[0];
rz(-0.70323479) q[0];
sx q[0];
rz(0.08589311) q[0];
rz(-0.15726451) q[2];
sx q[2];
rz(-2.3549821) q[2];
sx q[2];
rz(0.32341126) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7622457) q[1];
sx q[1];
rz(-1.8384461) q[1];
sx q[1];
rz(2.2257811) q[1];
x q[2];
rz(-0.31339733) q[3];
sx q[3];
rz(-2.0923231) q[3];
sx q[3];
rz(1.0115305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9909624) q[2];
sx q[2];
rz(-1.4272828) q[2];
sx q[2];
rz(2.4347351) q[2];
rz(1.3458716) q[3];
sx q[3];
rz(-1.6790877) q[3];
sx q[3];
rz(-2.4367387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83955806) q[0];
sx q[0];
rz(-2.3607881) q[0];
sx q[0];
rz(1.6395521) q[0];
rz(-1.543965) q[1];
sx q[1];
rz(-0.85999703) q[1];
sx q[1];
rz(-1.647321) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65022138) q[0];
sx q[0];
rz(-1.5724535) q[0];
sx q[0];
rz(2.1412388) q[0];
rz(-pi) q[1];
rz(-1.5780894) q[2];
sx q[2];
rz(-0.97541684) q[2];
sx q[2];
rz(-0.0059758107) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8710821) q[1];
sx q[1];
rz(-1.8746261) q[1];
sx q[1];
rz(0.31131502) q[1];
rz(-pi) q[2];
rz(-0.022541209) q[3];
sx q[3];
rz(-0.95571858) q[3];
sx q[3];
rz(-2.817077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26022252) q[2];
sx q[2];
rz(-1.2726731) q[2];
sx q[2];
rz(-0.55848813) q[2];
rz(0.50885606) q[3];
sx q[3];
rz(-1.5285834) q[3];
sx q[3];
rz(-1.2671027) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59638554) q[0];
sx q[0];
rz(-1.8983497) q[0];
sx q[0];
rz(-0.16045706) q[0];
rz(-0.67333737) q[1];
sx q[1];
rz(-1.2168177) q[1];
sx q[1];
rz(1.1268667) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9551487) q[0];
sx q[0];
rz(-1.2683682) q[0];
sx q[0];
rz(-1.1919829) q[0];
rz(0.33529386) q[2];
sx q[2];
rz(-0.79293434) q[2];
sx q[2];
rz(1.060134) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8452478) q[1];
sx q[1];
rz(-2.5242477) q[1];
sx q[1];
rz(-1.9250573) q[1];
rz(-pi) q[2];
rz(-2.1127719) q[3];
sx q[3];
rz(-1.3888956) q[3];
sx q[3];
rz(-2.101909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1274073) q[2];
sx q[2];
rz(-1.5994453) q[2];
sx q[2];
rz(0.7737774) q[2];
rz(2.1524147) q[3];
sx q[3];
rz(-2.7343605) q[3];
sx q[3];
rz(2.0686843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2268552) q[0];
sx q[0];
rz(-1.4829153) q[0];
sx q[0];
rz(-0.7731272) q[0];
rz(-1.5559366) q[1];
sx q[1];
rz(-1.8525886) q[1];
sx q[1];
rz(1.9645436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0814514) q[0];
sx q[0];
rz(-0.8340652) q[0];
sx q[0];
rz(-3.0990776) q[0];
x q[1];
rz(-2.8297205) q[2];
sx q[2];
rz(-2.1215028) q[2];
sx q[2];
rz(-2.4518509) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6625632) q[1];
sx q[1];
rz(-1.8445043) q[1];
sx q[1];
rz(-2.231313) q[1];
rz(1.7103393) q[3];
sx q[3];
rz(-1.7467032) q[3];
sx q[3];
rz(-1.6493451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3561463) q[2];
sx q[2];
rz(-0.36474228) q[2];
sx q[2];
rz(-2.6981603) q[2];
rz(2.7555079) q[3];
sx q[3];
rz(-1.4177136) q[3];
sx q[3];
rz(-0.20080565) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76483738) q[0];
sx q[0];
rz(-1.5926462) q[0];
sx q[0];
rz(3.1078597) q[0];
rz(2.4434166) q[1];
sx q[1];
rz(-0.59711421) q[1];
sx q[1];
rz(-2.4765292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8830983) q[0];
sx q[0];
rz(-1.8332714) q[0];
sx q[0];
rz(1.1061944) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4735293) q[2];
sx q[2];
rz(-2.2396834) q[2];
sx q[2];
rz(0.39448002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4951952) q[1];
sx q[1];
rz(-0.84719275) q[1];
sx q[1];
rz(-0.88314914) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.063015032) q[3];
sx q[3];
rz(-0.32013963) q[3];
sx q[3];
rz(-3.0885765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.72208059) q[2];
sx q[2];
rz(-0.91700143) q[2];
sx q[2];
rz(-1.7783995) q[2];
rz(-2.6204387) q[3];
sx q[3];
rz(-1.3241973) q[3];
sx q[3];
rz(-2.6496729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1727961) q[0];
sx q[0];
rz(-1.3461312) q[0];
sx q[0];
rz(0.81934339) q[0];
rz(0.68483812) q[1];
sx q[1];
rz(-1.2872773) q[1];
sx q[1];
rz(0.12231621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4222833) q[0];
sx q[0];
rz(-0.78351142) q[0];
sx q[0];
rz(2.1859403) q[0];
rz(-pi) q[1];
rz(-2.7610131) q[2];
sx q[2];
rz(-2.1316445) q[2];
sx q[2];
rz(2.3083936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51309359) q[1];
sx q[1];
rz(-0.6969531) q[1];
sx q[1];
rz(2.6234187) q[1];
x q[2];
rz(-2.0360394) q[3];
sx q[3];
rz(-2.1085582) q[3];
sx q[3];
rz(0.28013438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.43716064) q[2];
sx q[2];
rz(-1.8081343) q[2];
sx q[2];
rz(-2.2740347) q[2];
rz(1.7433172) q[3];
sx q[3];
rz(-2.7531392) q[3];
sx q[3];
rz(-0.5790264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45359465) q[0];
sx q[0];
rz(-1.2247676) q[0];
sx q[0];
rz(1.0892185) q[0];
rz(-0.047529686) q[1];
sx q[1];
rz(-1.0462953) q[1];
sx q[1];
rz(2.4542123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6733765) q[0];
sx q[0];
rz(-0.49096732) q[0];
sx q[0];
rz(-1.6223905) q[0];
rz(-1.1401165) q[2];
sx q[2];
rz(-1.8269369) q[2];
sx q[2];
rz(1.3821951) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2694305) q[1];
sx q[1];
rz(-0.80842847) q[1];
sx q[1];
rz(2.4926659) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0484588) q[3];
sx q[3];
rz(-2.2054005) q[3];
sx q[3];
rz(2.840691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1462732) q[2];
sx q[2];
rz(-0.61890382) q[2];
sx q[2];
rz(-1.4825561) q[2];
rz(-0.65685529) q[3];
sx q[3];
rz(-1.1494136) q[3];
sx q[3];
rz(2.759554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.252608) q[0];
sx q[0];
rz(-1.1532619) q[0];
sx q[0];
rz(1.1070195) q[0];
rz(-2.5629015) q[1];
sx q[1];
rz(-2.3042669) q[1];
sx q[1];
rz(-0.28490983) q[1];
rz(-0.22004057) q[2];
sx q[2];
rz(-0.36424988) q[2];
sx q[2];
rz(3.0314432) q[2];
rz(-1.4887078) q[3];
sx q[3];
rz(-1.0552867) q[3];
sx q[3];
rz(-1.1724006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
