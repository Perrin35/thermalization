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
rz(0.69005203) q[0];
sx q[0];
rz(-0.85059387) q[0];
sx q[0];
rz(-2.9414862) q[0];
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
rz(-2.645675) q[0];
sx q[0];
rz(-2.0089275) q[0];
sx q[0];
rz(-0.25610957) q[0];
rz(-pi) q[1];
rz(2.5426504) q[2];
sx q[2];
rz(-1.9859196) q[2];
sx q[2];
rz(0.51905635) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15439776) q[1];
sx q[1];
rz(-1.5227277) q[1];
sx q[1];
rz(0.67832077) q[1];
rz(-1.7302527) q[3];
sx q[3];
rz(-0.48122367) q[3];
sx q[3];
rz(-0.10018292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.072210463) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(-1.3145831) q[2];
rz(2.2383111) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(-2.7956853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2717993) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(-2.2858009) q[0];
rz(2.3675512) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(2.3537297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44875328) q[0];
sx q[0];
rz(-1.3200545) q[0];
sx q[0];
rz(3.1186799) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11779479) q[2];
sx q[2];
rz(-2.795145) q[2];
sx q[2];
rz(2.7730377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1545002) q[1];
sx q[1];
rz(-2.5804511) q[1];
sx q[1];
rz(2.2461056) q[1];
x q[2];
rz(-2.47452) q[3];
sx q[3];
rz(-1.4049585) q[3];
sx q[3];
rz(-1.1950243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.065980109) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(2.7454929) q[2];
rz(1.570545) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(-1.7861231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1876672) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(-3.0809825) q[0];
rz(2.4568779) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(-0.33904591) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44797541) q[0];
sx q[0];
rz(-2.1537499) q[0];
sx q[0];
rz(-1.5459803) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6940368) q[2];
sx q[2];
rz(-1.3831105) q[2];
sx q[2];
rz(2.0579684) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10968929) q[1];
sx q[1];
rz(-1.8432143) q[1];
sx q[1];
rz(-1.2485571) q[1];
rz(-0.011914201) q[3];
sx q[3];
rz(-2.3647226) q[3];
sx q[3];
rz(-0.35038951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48245779) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(0.47046146) q[2];
rz(0.86236924) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(-2.733574) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(-1.8203576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119638) q[0];
sx q[0];
rz(-2.5416059) q[0];
sx q[0];
rz(0.55498755) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91561134) q[2];
sx q[2];
rz(-1.5024937) q[2];
sx q[2];
rz(-3.0223522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8157573) q[1];
sx q[1];
rz(-1.6920631) q[1];
sx q[1];
rz(-1.6440637) q[1];
x q[2];
rz(-0.42709728) q[3];
sx q[3];
rz(-1.9582886) q[3];
sx q[3];
rz(2.1807293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2230175) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(2.4187386) q[2];
rz(2.2919848) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(-1.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1866813) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(0.55996672) q[0];
rz(0.2050744) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(2.8505039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2084853) q[0];
sx q[0];
rz(-1.7400842) q[0];
sx q[0];
rz(-0.27219682) q[0];
rz(1.5076324) q[2];
sx q[2];
rz(-2.3690555) q[2];
sx q[2];
rz(2.87839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4704125) q[1];
sx q[1];
rz(-2.414783) q[1];
sx q[1];
rz(1.8720759) q[1];
rz(-pi) q[2];
rz(1.4104615) q[3];
sx q[3];
rz(-0.82538) q[3];
sx q[3];
rz(1.7602518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2228955) q[2];
sx q[2];
rz(-1.6632068) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8339612) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(0.6024012) q[0];
rz(0.3715474) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(-2.1098302) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87001291) q[0];
sx q[0];
rz(-1.3965415) q[0];
sx q[0];
rz(-0.016079024) q[0];
rz(-pi) q[1];
rz(-1.7784836) q[2];
sx q[2];
rz(-1.1738994) q[2];
sx q[2];
rz(2.7092421) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25615869) q[1];
sx q[1];
rz(-2.6226461) q[1];
sx q[1];
rz(0.38796723) q[1];
x q[2];
rz(2.2183275) q[3];
sx q[3];
rz(-2.0482691) q[3];
sx q[3];
rz(-1.5058194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80456698) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(1.5186914) q[2];
rz(-2.2931781) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(-3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0254211) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(-0.99789944) q[0];
rz(-0.2746703) q[1];
sx q[1];
rz(-2.2614567) q[1];
sx q[1];
rz(-1.3708699) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8664198) q[0];
sx q[0];
rz(-2.201797) q[0];
sx q[0];
rz(1.5541398) q[0];
x q[1];
rz(2.9878658) q[2];
sx q[2];
rz(-1.9228553) q[2];
sx q[2];
rz(1.597126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22860195) q[1];
sx q[1];
rz(-0.96768846) q[1];
sx q[1];
rz(-1.8643537) q[1];
rz(-pi) q[2];
rz(-0.73565817) q[3];
sx q[3];
rz(-0.94649411) q[3];
sx q[3];
rz(-1.5746547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2844598) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(3.1119463) q[2];
rz(2.5263785) q[3];
sx q[3];
rz(-1.3968202) q[3];
sx q[3];
rz(-2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089559473) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(0.7830559) q[0];
rz(-1.1198593) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(1.7436183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92358855) q[0];
sx q[0];
rz(-2.5177885) q[0];
sx q[0];
rz(-0.82360928) q[0];
rz(0.2335642) q[2];
sx q[2];
rz(-1.2937045) q[2];
sx q[2];
rz(-1.8943) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74136342) q[1];
sx q[1];
rz(-1.5796164) q[1];
sx q[1];
rz(-0.90290248) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3574202) q[3];
sx q[3];
rz(-2.3436574) q[3];
sx q[3];
rz(2.7415581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(1.3207377) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(-1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9418697) q[0];
sx q[0];
rz(-1.0933192) q[0];
sx q[0];
rz(-2.1871908) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(2.0571713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0592151) q[0];
sx q[0];
rz(-0.70967662) q[0];
sx q[0];
rz(-1.6355455) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2899288) q[2];
sx q[2];
rz(-1.3251588) q[2];
sx q[2];
rz(-0.40058595) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.045001205) q[1];
sx q[1];
rz(-0.47811478) q[1];
sx q[1];
rz(-1.4131312) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7550657) q[3];
sx q[3];
rz(-2.3109155) q[3];
sx q[3];
rz(2.4746462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0584917) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(-2.7044738) q[2];
rz(1.3433749) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(1.7310671) q[0];
rz(2.9208185) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-2.208362) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8338884) q[0];
sx q[0];
rz(-1.5571463) q[0];
sx q[0];
rz(-2.3962014) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5228517) q[2];
sx q[2];
rz(-1.2729984) q[2];
sx q[2];
rz(2.1190475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9737967) q[1];
sx q[1];
rz(-2.3575511) q[1];
sx q[1];
rz(1.7586437) q[1];
rz(-1.4187063) q[3];
sx q[3];
rz(-2.503431) q[3];
sx q[3];
rz(-1.1537976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.71020469) q[2];
sx q[2];
rz(-2.7993671) q[2];
sx q[2];
rz(-0.63146511) q[2];
rz(2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(-2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5030293) q[0];
sx q[0];
rz(-1.6521492) q[0];
sx q[0];
rz(-1.8304992) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(3.1099215) q[2];
sx q[2];
rz(-1.1902255) q[2];
sx q[2];
rz(-2.5524216) q[2];
rz(1.4680223) q[3];
sx q[3];
rz(-1.0036052) q[3];
sx q[3];
rz(-2.6617692) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
