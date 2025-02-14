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
rz(-2.4515406) q[0];
sx q[0];
rz(-2.2909988) q[0];
sx q[0];
rz(2.9414862) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(5.7601647) q[1];
sx q[1];
rz(10.756607) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1855186) q[0];
sx q[0];
rz(-1.8022493) q[0];
sx q[0];
rz(1.1197907) q[0];
x q[1];
rz(2.5426504) q[2];
sx q[2];
rz(-1.9859196) q[2];
sx q[2];
rz(-2.6225363) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6864904) q[1];
sx q[1];
rz(-2.2481866) q[1];
sx q[1];
rz(-1.6325006) q[1];
rz(-0.08272127) q[3];
sx q[3];
rz(-2.045407) q[3];
sx q[3];
rz(0.079291346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.072210463) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(-1.3145831) q[2];
rz(2.2383111) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697934) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(2.2858009) q[0];
rz(-0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(0.78786293) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.60074) q[0];
sx q[0];
rz(-2.8898281) q[0];
sx q[0];
rz(1.4815848) q[0];
x q[1];
rz(1.5283952) q[2];
sx q[2];
rz(-1.2268492) q[2];
sx q[2];
rz(-0.49371546) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.012177906) q[1];
sx q[1];
rz(-1.9099292) q[1];
sx q[1];
rz(-2.0268834) q[1];
x q[2];
rz(0.66707261) q[3];
sx q[3];
rz(-1.4049585) q[3];
sx q[3];
rz(-1.1950243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.065980109) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(2.7454929) q[2];
rz(1.570545) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(1.3554696) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1876672) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(3.0809825) q[0];
rz(0.68471471) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(2.8025467) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44797541) q[0];
sx q[0];
rz(-2.1537499) q[0];
sx q[0];
rz(1.5459803) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57444797) q[2];
sx q[2];
rz(-0.22413218) q[2];
sx q[2];
rz(1.6696861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0319034) q[1];
sx q[1];
rz(-1.2983783) q[1];
sx q[1];
rz(-1.2485571) q[1];
rz(-pi) q[2];
rz(0.011914201) q[3];
sx q[3];
rz(-0.77687009) q[3];
sx q[3];
rz(2.7912031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48245779) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(2.6711312) q[2];
rz(0.86236924) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(-1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(-2.733574) q[0];
rz(1.3511924) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(-1.8203576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6141339) q[0];
sx q[0];
rz(-1.2686994) q[0];
sx q[0];
rz(-0.52665773) q[0];
x q[1];
rz(3.0555326) q[2];
sx q[2];
rz(-0.91740184) q[2];
sx q[2];
rz(-1.5039521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76391782) q[1];
sx q[1];
rz(-1.498068) q[1];
sx q[1];
rz(3.0200028) q[1];
x q[2];
rz(-1.1492689) q[3];
sx q[3];
rz(-1.1771923) q[3];
sx q[3];
rz(-0.43962653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(-0.72285405) q[2];
rz(0.84960788) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(-2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(2.5816259) q[0];
rz(-2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(-0.29108873) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2084853) q[0];
sx q[0];
rz(-1.4015084) q[0];
sx q[0];
rz(2.8693958) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6339602) q[2];
sx q[2];
rz(-0.77253714) q[2];
sx q[2];
rz(-0.26320266) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6711802) q[1];
sx q[1];
rz(-2.414783) q[1];
sx q[1];
rz(1.2695168) q[1];
x q[2];
rz(1.4104615) q[3];
sx q[3];
rz(-0.82538) q[3];
sx q[3];
rz(-1.3813409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9186972) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(-2.843294) q[2];
rz(1.1693303) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.8339612) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(-2.5391915) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(1.0317624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96249798) q[0];
sx q[0];
rz(-0.17498762) q[0];
sx q[0];
rz(-1.4797158) q[0];
rz(-2.736892) q[2];
sx q[2];
rz(-1.7621303) q[2];
sx q[2];
rz(1.9218685) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.885434) q[1];
sx q[1];
rz(-0.51894655) q[1];
sx q[1];
rz(0.38796723) q[1];
rz(-0.86182819) q[3];
sx q[3];
rz(-2.3579881) q[3];
sx q[3];
rz(2.6605822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3370257) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(1.5186914) q[2];
rz(-2.2931781) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(-0.039610473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1161716) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(-2.1436932) q[0];
rz(-0.2746703) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.3708699) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8664198) q[0];
sx q[0];
rz(-2.201797) q[0];
sx q[0];
rz(-1.5874528) q[0];
rz(-1.9657126) q[2];
sx q[2];
rz(-0.38287258) q[2];
sx q[2];
rz(-1.9667786) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6295885) q[1];
sx q[1];
rz(-1.8114212) q[1];
sx q[1];
rz(2.5179067) q[1];
rz(2.3207449) q[3];
sx q[3];
rz(-0.92536649) q[3];
sx q[3];
rz(-0.57725932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85713282) q[2];
sx q[2];
rz(-0.94228116) q[2];
sx q[2];
rz(0.029646309) q[2];
rz(2.5263785) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(-0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089559473) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(-0.7830559) q[0];
rz(1.1198593) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(1.7436183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92358855) q[0];
sx q[0];
rz(-2.5177885) q[0];
sx q[0];
rz(0.82360928) q[0];
rz(-1.8552165) q[2];
sx q[2];
rz(-1.3462974) q[2];
sx q[2];
rz(2.7530991) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82247558) q[1];
sx q[1];
rz(-0.90293316) q[1];
sx q[1];
rz(-3.130359) q[1];
rz(0.78432958) q[3];
sx q[3];
rz(-1.722986) q[3];
sx q[3];
rz(1.0206211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.3207377) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(2.1871908) q[0];
rz(-1.154254) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(2.0571713) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53756489) q[0];
sx q[0];
rz(-1.5286235) q[0];
sx q[0];
rz(-0.86215638) q[0];
rz(1.2899288) q[2];
sx q[2];
rz(-1.8164338) q[2];
sx q[2];
rz(-2.7410067) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22218765) q[1];
sx q[1];
rz(-2.0424941) q[1];
sx q[1];
rz(-3.0604048) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7550657) q[3];
sx q[3];
rz(-2.3109155) q[3];
sx q[3];
rz(-2.4746462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.083100975) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(-2.7044738) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(-3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(-1.7310671) q[0];
rz(-0.22077416) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-0.9332307) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30770424) q[0];
sx q[0];
rz(-1.5571463) q[0];
sx q[0];
rz(-2.3962014) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6548806) q[2];
sx q[2];
rz(-0.67811869) q[2];
sx q[2];
rz(2.202452) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9737967) q[1];
sx q[1];
rz(-2.3575511) q[1];
sx q[1];
rz(1.7586437) q[1];
rz(-pi) q[2];
rz(-0.93818061) q[3];
sx q[3];
rz(-1.6611735) q[3];
sx q[3];
rz(0.5394906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.71020469) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(-0.63146511) q[2];
rz(-0.7190052) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(-0.24148153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5030293) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(-1.491811) q[2];
sx q[2];
rz(-0.38182237) q[2];
sx q[2];
rz(0.67425722) q[2];
rz(1.6735703) q[3];
sx q[3];
rz(-2.1379875) q[3];
sx q[3];
rz(0.47982346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
