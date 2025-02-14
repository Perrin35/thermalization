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
rz(-2.8430804) q[0];
sx q[0];
rz(-1.2895583) q[0];
sx q[0];
rz(3.1182752) q[0];
rz(0.98195568) q[1];
sx q[1];
rz(-0.73749956) q[1];
sx q[1];
rz(1.4758543) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6029089) q[0];
sx q[0];
rz(-1.6404466) q[0];
sx q[0];
rz(3.0979741) q[0];
rz(-pi) q[1];
rz(3.1013203) q[2];
sx q[2];
rz(-1.9007517) q[2];
sx q[2];
rz(1.9391935) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75778537) q[1];
sx q[1];
rz(-1.1769036) q[1];
sx q[1];
rz(-0.76683137) q[1];
x q[2];
rz(-0.89453434) q[3];
sx q[3];
rz(-1.9893551) q[3];
sx q[3];
rz(2.7194617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4618571) q[2];
sx q[2];
rz(-1.4718141) q[2];
sx q[2];
rz(-2.8409345) q[2];
rz(0.65827185) q[3];
sx q[3];
rz(-2.7418147) q[3];
sx q[3];
rz(2.5532653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7268426) q[0];
sx q[0];
rz(-2.2756133) q[0];
sx q[0];
rz(-0.17587371) q[0];
rz(0.56100065) q[1];
sx q[1];
rz(-2.439552) q[1];
sx q[1];
rz(-0.19634136) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25717673) q[0];
sx q[0];
rz(-1.1286583) q[0];
sx q[0];
rz(0.37143645) q[0];
rz(-pi) q[1];
rz(-0.11949338) q[2];
sx q[2];
rz(-1.6867078) q[2];
sx q[2];
rz(1.3136065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2159816) q[1];
sx q[1];
rz(-0.51541939) q[1];
sx q[1];
rz(-1.5748596) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34336932) q[3];
sx q[3];
rz(-1.6968622) q[3];
sx q[3];
rz(2.8186563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.089513) q[2];
sx q[2];
rz(-0.61605805) q[2];
sx q[2];
rz(1.6717795) q[2];
rz(0.52516627) q[3];
sx q[3];
rz(-1.7142121) q[3];
sx q[3];
rz(0.66450459) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0199652) q[0];
sx q[0];
rz(-2.2522734) q[0];
sx q[0];
rz(-2.5110733) q[0];
rz(1.9969253) q[1];
sx q[1];
rz(-1.2455218) q[1];
sx q[1];
rz(-0.05680457) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6358179) q[0];
sx q[0];
rz(-0.56418252) q[0];
sx q[0];
rz(0.83585541) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6654451) q[2];
sx q[2];
rz(-2.2000929) q[2];
sx q[2];
rz(-0.72222939) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60384974) q[1];
sx q[1];
rz(-1.7122388) q[1];
sx q[1];
rz(-0.46700041) q[1];
x q[2];
rz(0.2127368) q[3];
sx q[3];
rz(-1.1520583) q[3];
sx q[3];
rz(-1.1513082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28222617) q[2];
sx q[2];
rz(-1.4713918) q[2];
sx q[2];
rz(-2.7024506) q[2];
rz(1.3965083) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(-0.5784353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0466995) q[0];
sx q[0];
rz(-2.4201604) q[0];
sx q[0];
rz(-1.0149581) q[0];
rz(3.0789442) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(0.28422022) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361523) q[0];
sx q[0];
rz(-1.9198737) q[0];
sx q[0];
rz(2.6830868) q[0];
x q[1];
rz(-2.6437628) q[2];
sx q[2];
rz(-1.4477028) q[2];
sx q[2];
rz(-0.13409889) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9584459) q[1];
sx q[1];
rz(-1.5827521) q[1];
sx q[1];
rz(-0.35838106) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.013945097) q[3];
sx q[3];
rz(-1.317465) q[3];
sx q[3];
rz(-1.0903181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8404428) q[2];
sx q[2];
rz(-2.132685) q[2];
sx q[2];
rz(-0.76048771) q[2];
rz(-2.8216951) q[3];
sx q[3];
rz(-1.1287929) q[3];
sx q[3];
rz(-1.0285876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29192057) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(2.736295) q[0];
rz(-2.6536476) q[1];
sx q[1];
rz(-1.1064203) q[1];
sx q[1];
rz(-2.778756) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80984801) q[0];
sx q[0];
rz(-0.54143006) q[0];
sx q[0];
rz(-0.62744139) q[0];
x q[1];
rz(2.6020701) q[2];
sx q[2];
rz(-2.4837257) q[2];
sx q[2];
rz(0.69275975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.664714) q[1];
sx q[1];
rz(-1.6194317) q[1];
sx q[1];
rz(0.22108476) q[1];
rz(-pi) q[2];
rz(-1.8764349) q[3];
sx q[3];
rz(-2.2430217) q[3];
sx q[3];
rz(1.6068939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.001658) q[2];
sx q[2];
rz(-1.9826865) q[2];
sx q[2];
rz(-0.72251594) q[2];
rz(-2.6288988) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(1.9374013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2672511) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(2.8668218) q[0];
rz(-2.784506) q[1];
sx q[1];
rz(-1.5029224) q[1];
sx q[1];
rz(1.9836609) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0561075) q[0];
sx q[0];
rz(-2.4153313) q[0];
sx q[0];
rz(1.1994491) q[0];
x q[1];
rz(0.59717565) q[2];
sx q[2];
rz(-2.1417326) q[2];
sx q[2];
rz(-1.939029) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8573945) q[1];
sx q[1];
rz(-1.3669087) q[1];
sx q[1];
rz(0.60867761) q[1];
rz(-0.88648667) q[3];
sx q[3];
rz(-2.2456256) q[3];
sx q[3];
rz(1.0533028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.57924119) q[2];
sx q[2];
rz(-1.3902105) q[2];
sx q[2];
rz(0.58970279) q[2];
rz(-1.8288745) q[3];
sx q[3];
rz(-1.4785942) q[3];
sx q[3];
rz(0.5635128) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8967459) q[0];
sx q[0];
rz(-2.1488996) q[0];
sx q[0];
rz(-2.8712811) q[0];
rz(-0.42701328) q[1];
sx q[1];
rz(-1.7662363) q[1];
sx q[1];
rz(0.40831533) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1751375) q[0];
sx q[0];
rz(-1.2027367) q[0];
sx q[0];
rz(1.6987263) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2683581) q[2];
sx q[2];
rz(-1.6484652) q[2];
sx q[2];
rz(0.49688646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9198614) q[1];
sx q[1];
rz(-1.9834922) q[1];
sx q[1];
rz(-2.2151715) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72260489) q[3];
sx q[3];
rz(-1.4935023) q[3];
sx q[3];
rz(-2.3376436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.938505) q[2];
sx q[2];
rz(-1.0838584) q[2];
sx q[2];
rz(-2.2173524) q[2];
rz(2.3067394) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(2.1448962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860109) q[0];
sx q[0];
rz(-0.38014933) q[0];
sx q[0];
rz(-0.36744776) q[0];
rz(0.5683178) q[1];
sx q[1];
rz(-1.808993) q[1];
sx q[1];
rz(0.41499358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5552705) q[0];
sx q[0];
rz(-1.30997) q[0];
sx q[0];
rz(-1.1010948) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6122323) q[2];
sx q[2];
rz(-2.0545914) q[2];
sx q[2];
rz(-0.31365487) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4492682) q[1];
sx q[1];
rz(-1.0639166) q[1];
sx q[1];
rz(-0.70517069) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9165809) q[3];
sx q[3];
rz(-1.5756851) q[3];
sx q[3];
rz(-1.9859615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35855287) q[2];
sx q[2];
rz(-0.8873322) q[2];
sx q[2];
rz(-0.83317327) q[2];
rz(0.35946515) q[3];
sx q[3];
rz(-1.0901674) q[3];
sx q[3];
rz(1.5045213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1737162) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(-0.014884431) q[0];
rz(2.0170276) q[1];
sx q[1];
rz(-2.896307) q[1];
sx q[1];
rz(-3.0344149) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8421321) q[0];
sx q[0];
rz(-2.3981574) q[0];
sx q[0];
rz(1.3495665) q[0];
rz(-pi) q[1];
rz(0.82294269) q[2];
sx q[2];
rz(-1.2185893) q[2];
sx q[2];
rz(2.4885881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1746782) q[1];
sx q[1];
rz(-2.0474829) q[1];
sx q[1];
rz(-1.7382453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.393287) q[3];
sx q[3];
rz(-1.8281947) q[3];
sx q[3];
rz(-2.7696115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5643481) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(-0.85961071) q[2];
rz(-0.25935069) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(-0.41659659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-0.49941007) q[0];
sx q[0];
rz(-3.0152617) q[0];
sx q[0];
rz(-1.9835749) q[0];
rz(2.6462789) q[1];
sx q[1];
rz(-1.9751578) q[1];
sx q[1];
rz(-2.8445629) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.934108) q[0];
sx q[0];
rz(-1.5667324) q[0];
sx q[0];
rz(1.8898176) q[0];
x q[1];
rz(-2.8412965) q[2];
sx q[2];
rz(-1.5321333) q[2];
sx q[2];
rz(3.0186332) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9969284) q[1];
sx q[1];
rz(-1.7283617) q[1];
sx q[1];
rz(-2.8466606) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2333651) q[3];
sx q[3];
rz(-1.55791) q[3];
sx q[3];
rz(0.21256178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.872252) q[2];
sx q[2];
rz(-2.7597235) q[2];
sx q[2];
rz(-2.2807109) q[2];
rz(2.6779029) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(2.004682) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1562445) q[0];
sx q[0];
rz(-1.5679659) q[0];
sx q[0];
rz(1.1563942) q[0];
rz(-1.3337749) q[1];
sx q[1];
rz(-1.540624) q[1];
sx q[1];
rz(-2.435138) q[1];
rz(-2.479066) q[2];
sx q[2];
rz(-2.0338197) q[2];
sx q[2];
rz(2.8070569) q[2];
rz(-2.4931277) q[3];
sx q[3];
rz(-0.59352661) q[3];
sx q[3];
rz(1.2980579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
