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
rz(1.0951618) q[0];
sx q[0];
rz(0.14552966) q[0];
sx q[0];
rz(9.9280823) q[0];
rz(1.1777999) q[1];
sx q[1];
rz(-1.5947394) q[1];
sx q[1];
rz(-1.6961179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58848042) q[0];
sx q[0];
rz(-2.7403826) q[0];
sx q[0];
rz(0.11102872) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4769449) q[2];
sx q[2];
rz(-2.2713813) q[2];
sx q[2];
rz(1.0657016) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9564285) q[1];
sx q[1];
rz(-1.899029) q[1];
sx q[1];
rz(2.4976455) q[1];
rz(-pi) q[2];
rz(1.5551994) q[3];
sx q[3];
rz(-1.6199779) q[3];
sx q[3];
rz(1.552207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4989) q[2];
sx q[2];
rz(-1.2314726) q[2];
sx q[2];
rz(-1.5042245) q[2];
rz(0.73689342) q[3];
sx q[3];
rz(-1.6156018) q[3];
sx q[3];
rz(2.1604497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7351643) q[0];
sx q[0];
rz(-0.46877113) q[0];
sx q[0];
rz(2.6306187) q[0];
rz(1.3276395) q[1];
sx q[1];
rz(-1.7402486) q[1];
sx q[1];
rz(-2.8151292) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8998979) q[0];
sx q[0];
rz(-2.7621671) q[0];
sx q[0];
rz(2.3510758) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9316177) q[2];
sx q[2];
rz(-0.89824235) q[2];
sx q[2];
rz(-0.58174101) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67124704) q[1];
sx q[1];
rz(-1.874028) q[1];
sx q[1];
rz(2.1821351) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94989461) q[3];
sx q[3];
rz(-2.6767459) q[3];
sx q[3];
rz(-0.45528938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9713126) q[2];
sx q[2];
rz(-2.9176517) q[2];
sx q[2];
rz(-1.2962606) q[2];
rz(-0.41804677) q[3];
sx q[3];
rz(-0.97801912) q[3];
sx q[3];
rz(-0.70438284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.068785) q[0];
sx q[0];
rz(-0.3827706) q[0];
sx q[0];
rz(2.5991154) q[0];
rz(-1.3756649) q[1];
sx q[1];
rz(-0.41789564) q[1];
sx q[1];
rz(-3.1105522) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9144273) q[0];
sx q[0];
rz(-1.7591159) q[0];
sx q[0];
rz(-0.31115599) q[0];
rz(-1.8607742) q[2];
sx q[2];
rz(-0.6781247) q[2];
sx q[2];
rz(-1.6999619) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.52306226) q[1];
sx q[1];
rz(-1.3253322) q[1];
sx q[1];
rz(2.1851468) q[1];
x q[2];
rz(1.5872845) q[3];
sx q[3];
rz(-2.1257536) q[3];
sx q[3];
rz(2.6944427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62640181) q[2];
sx q[2];
rz(-0.73651892) q[2];
sx q[2];
rz(-2.314563) q[2];
rz(-2.6229897) q[3];
sx q[3];
rz(-2.0293472) q[3];
sx q[3];
rz(1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14722918) q[0];
sx q[0];
rz(-1.8356859) q[0];
sx q[0];
rz(1.9299141) q[0];
rz(1.0481102) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(-1.3406219) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3163534) q[0];
sx q[0];
rz(-1.8779426) q[0];
sx q[0];
rz(0.65748837) q[0];
x q[1];
rz(-0.2366613) q[2];
sx q[2];
rz(-1.7156895) q[2];
sx q[2];
rz(1.7828816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2681899) q[1];
sx q[1];
rz(-0.65618578) q[1];
sx q[1];
rz(-0.49439685) q[1];
x q[2];
rz(1.9817438) q[3];
sx q[3];
rz(-1.8576429) q[3];
sx q[3];
rz(0.2528068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94053215) q[2];
sx q[2];
rz(-0.17720711) q[2];
sx q[2];
rz(2.000467) q[2];
rz(-1.5308258) q[3];
sx q[3];
rz(-0.96735668) q[3];
sx q[3];
rz(-2.358986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56213266) q[0];
sx q[0];
rz(-1.9714404) q[0];
sx q[0];
rz(-0.60254565) q[0];
rz(-2.5613979) q[1];
sx q[1];
rz(-1.5126901) q[1];
sx q[1];
rz(-2.3604732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21602042) q[0];
sx q[0];
rz(-1.4013774) q[0];
sx q[0];
rz(1.2792619) q[0];
rz(-pi) q[1];
rz(2.4915472) q[2];
sx q[2];
rz(-2.263265) q[2];
sx q[2];
rz(-1.4743069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8210817) q[1];
sx q[1];
rz(-2.2901504) q[1];
sx q[1];
rz(-3.0266552) q[1];
x q[2];
rz(0.1991462) q[3];
sx q[3];
rz(-2.4059787) q[3];
sx q[3];
rz(0.74935952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.96131229) q[2];
sx q[2];
rz(-1.140241) q[2];
sx q[2];
rz(-2.9126634) q[2];
rz(1.3198352) q[3];
sx q[3];
rz(-1.6596551) q[3];
sx q[3];
rz(2.2975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.22885403) q[0];
sx q[0];
rz(-1.5541394) q[0];
sx q[0];
rz(-1.9629021) q[0];
rz(-1.9121869) q[1];
sx q[1];
rz(-0.39437672) q[1];
sx q[1];
rz(1.2961402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79149198) q[0];
sx q[0];
rz(-1.434636) q[0];
sx q[0];
rz(-1.2873935) q[0];
rz(-pi) q[1];
rz(2.3713114) q[2];
sx q[2];
rz(-0.89777032) q[2];
sx q[2];
rz(1.9865303) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7414822) q[1];
sx q[1];
rz(-2.0026836) q[1];
sx q[1];
rz(-0.53439616) q[1];
rz(-pi) q[2];
rz(2.9518106) q[3];
sx q[3];
rz(-0.49473195) q[3];
sx q[3];
rz(3.1056946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65427762) q[2];
sx q[2];
rz(-0.9026022) q[2];
sx q[2];
rz(-2.8793092) q[2];
rz(2.8413963) q[3];
sx q[3];
rz(-1.2223949) q[3];
sx q[3];
rz(-2.9330971) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4406141) q[0];
sx q[0];
rz(-0.3955667) q[0];
sx q[0];
rz(-1.8390919) q[0];
rz(-2.473096) q[1];
sx q[1];
rz(-1.3553456) q[1];
sx q[1];
rz(1.0350593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.221091) q[0];
sx q[0];
rz(-1.971444) q[0];
sx q[0];
rz(-2.9095838) q[0];
rz(-1.9403753) q[2];
sx q[2];
rz(-1.3962922) q[2];
sx q[2];
rz(-2.7480887) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0332843) q[1];
sx q[1];
rz(-0.35221812) q[1];
sx q[1];
rz(-1.7914821) q[1];
rz(-pi) q[2];
rz(-0.74372282) q[3];
sx q[3];
rz(-1.111314) q[3];
sx q[3];
rz(2.2032025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0237026) q[2];
sx q[2];
rz(-1.002545) q[2];
sx q[2];
rz(-1.6560076) q[2];
rz(0.28247908) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(0.43454596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.78541237) q[0];
sx q[0];
rz(-1.0419351) q[0];
sx q[0];
rz(-2.8670512) q[0];
rz(1.2046332) q[1];
sx q[1];
rz(-0.58012539) q[1];
sx q[1];
rz(2.5801632) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0208541) q[0];
sx q[0];
rz(-1.0444114) q[0];
sx q[0];
rz(2.3096183) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37158575) q[2];
sx q[2];
rz(-2.4306261) q[2];
sx q[2];
rz(-1.154605) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16374396) q[1];
sx q[1];
rz(-1.1625966) q[1];
sx q[1];
rz(2.8099437) q[1];
x q[2];
rz(0.26430126) q[3];
sx q[3];
rz(-1.9076548) q[3];
sx q[3];
rz(1.565153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5795035) q[2];
sx q[2];
rz(-1.1056489) q[2];
sx q[2];
rz(1.7378463) q[2];
rz(1.3459407) q[3];
sx q[3];
rz(-0.74508777) q[3];
sx q[3];
rz(-2.6534206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76029921) q[0];
sx q[0];
rz(-2.5385222) q[0];
sx q[0];
rz(2.2316933) q[0];
rz(-1.1970041) q[1];
sx q[1];
rz(-2.0884114) q[1];
sx q[1];
rz(-0.88814703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67535454) q[0];
sx q[0];
rz(-2.9976237) q[0];
sx q[0];
rz(0.86021249) q[0];
rz(0.65983628) q[2];
sx q[2];
rz(-1.6183231) q[2];
sx q[2];
rz(-2.8328676) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17287317) q[1];
sx q[1];
rz(-0.76593562) q[1];
sx q[1];
rz(-0.57491803) q[1];
rz(-pi) q[2];
rz(-3.1040583) q[3];
sx q[3];
rz(-1.8662631) q[3];
sx q[3];
rz(0.42643828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1278648) q[2];
sx q[2];
rz(-1.4415386) q[2];
sx q[2];
rz(1.0825276) q[2];
rz(-0.23076375) q[3];
sx q[3];
rz(-1.8133546) q[3];
sx q[3];
rz(2.6575991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3103264) q[0];
sx q[0];
rz(-2.3895097) q[0];
sx q[0];
rz(2.5600774) q[0];
rz(-0.61559081) q[1];
sx q[1];
rz(-2.1938727) q[1];
sx q[1];
rz(1.8448578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49823495) q[0];
sx q[0];
rz(-1.5811992) q[0];
sx q[0];
rz(1.5889945) q[0];
x q[1];
rz(-1.5803976) q[2];
sx q[2];
rz(-1.2266955) q[2];
sx q[2];
rz(2.8167958) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2263068) q[1];
sx q[1];
rz(-1.9939594) q[1];
sx q[1];
rz(2.8598815) q[1];
rz(-pi) q[2];
rz(-2.4695685) q[3];
sx q[3];
rz(-1.250953) q[3];
sx q[3];
rz(-0.32957382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41946188) q[2];
sx q[2];
rz(-3.0007134) q[2];
sx q[2];
rz(2.6922743) q[2];
rz(0.32014534) q[3];
sx q[3];
rz(-1.82205) q[3];
sx q[3];
rz(-1.8951353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37353361) q[0];
sx q[0];
rz(-0.67714416) q[0];
sx q[0];
rz(-1.9376391) q[0];
rz(-1.7569348) q[1];
sx q[1];
rz(-0.23778267) q[1];
sx q[1];
rz(2.3429088) q[1];
rz(-1.9282707) q[2];
sx q[2];
rz(-2.4926179) q[2];
sx q[2];
rz(0.31821584) q[2];
rz(-2.2438335) q[3];
sx q[3];
rz(-2.5644292) q[3];
sx q[3];
rz(0.226365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
