OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(-2.7859935) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0922001) q[0];
sx q[0];
rz(-2.85673) q[0];
sx q[0];
rz(1.0546791) q[0];
rz(-pi) q[1];
rz(1.7945292) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(-1.7681233) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8259748) q[1];
sx q[1];
rz(-0.74233913) q[1];
sx q[1];
rz(-0.62645341) q[1];
x q[2];
rz(1.7765331) q[3];
sx q[3];
rz(-1.7997026) q[3];
sx q[3];
rz(2.8347901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(0.65650666) q[2];
rz(2.4025829) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332557) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(2.546229) q[0];
rz(-3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-0.48746902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90855234) q[0];
sx q[0];
rz(-1.4899583) q[0];
sx q[0];
rz(3.0768865) q[0];
rz(-pi) q[1];
rz(2.7198118) q[2];
sx q[2];
rz(-1.1676719) q[2];
sx q[2];
rz(-3.1096854) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6007874) q[1];
sx q[1];
rz(-1.4122737) q[1];
sx q[1];
rz(-1.7608789) q[1];
x q[2];
rz(0.08667605) q[3];
sx q[3];
rz(-1.8868586) q[3];
sx q[3];
rz(-2.812127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(2.9188459) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4091464) q[0];
sx q[0];
rz(-1.1799066) q[0];
sx q[0];
rz(1.976165) q[0];
rz(-pi) q[1];
rz(1.9420625) q[2];
sx q[2];
rz(-1.7811333) q[2];
sx q[2];
rz(0.51031761) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5388515) q[1];
sx q[1];
rz(-2.1826084) q[1];
sx q[1];
rz(-1.6536504) q[1];
x q[2];
rz(0.24554159) q[3];
sx q[3];
rz(-1.2430827) q[3];
sx q[3];
rz(-1.4471734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(1.1941236) q[2];
rz(2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38917437) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-2.2213675) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3738149) q[0];
sx q[0];
rz(-0.91705634) q[0];
sx q[0];
rz(2.6184404) q[0];
rz(-pi) q[1];
rz(0.36107365) q[2];
sx q[2];
rz(-1.4898584) q[2];
sx q[2];
rz(1.0635478) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8060382) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(2.4469417) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7761049) q[3];
sx q[3];
rz(-1.9234386) q[3];
sx q[3];
rz(0.40311381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(-1.8614004) q[2];
rz(0.81930339) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-0.4531025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4015409) q[0];
sx q[0];
rz(-1.4972367) q[0];
sx q[0];
rz(3.1042276) q[0];
rz(-3.1066936) q[2];
sx q[2];
rz(-1.4279832) q[2];
sx q[2];
rz(0.81894433) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15935005) q[1];
sx q[1];
rz(-2.9061926) q[1];
sx q[1];
rz(1.2118641) q[1];
x q[2];
rz(0.94282486) q[3];
sx q[3];
rz(-2.4403265) q[3];
sx q[3];
rz(-0.065746106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033427514) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(0.23813716) q[0];
rz(-0.37995964) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(1.628081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1171922) q[0];
sx q[0];
rz(-1.2478561) q[0];
sx q[0];
rz(-0.49844235) q[0];
rz(-pi) q[1];
rz(-1.6983301) q[2];
sx q[2];
rz(-2.6886534) q[2];
sx q[2];
rz(2.4085101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39735079) q[1];
sx q[1];
rz(-2.2480818) q[1];
sx q[1];
rz(-2.1703297) q[1];
rz(-pi) q[2];
rz(0.38576491) q[3];
sx q[3];
rz(-1.1499975) q[3];
sx q[3];
rz(-2.5667131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.05904077) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(-2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1259574) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(2.7440199) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.1869173) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884739) q[0];
sx q[0];
rz(-1.0736199) q[0];
sx q[0];
rz(2.8209646) q[0];
rz(-pi) q[1];
rz(1.5415807) q[2];
sx q[2];
rz(-1.511682) q[2];
sx q[2];
rz(0.37994775) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6656875) q[1];
sx q[1];
rz(-0.74531065) q[1];
sx q[1];
rz(1.3728632) q[1];
rz(2.7508221) q[3];
sx q[3];
rz(-2.4739389) q[3];
sx q[3];
rz(0.30927502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.879803) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.7857893) q[2];
rz(-1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-2.1112679) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54661575) q[0];
sx q[0];
rz(-2.3175276) q[0];
sx q[0];
rz(-1.0975518) q[0];
x q[1];
rz(2.0975935) q[2];
sx q[2];
rz(-1.5287405) q[2];
sx q[2];
rz(-2.2823357) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6701263) q[1];
sx q[1];
rz(-1.2738436) q[1];
sx q[1];
rz(2.9832277) q[1];
rz(0.52328531) q[3];
sx q[3];
rz(-1.0417632) q[3];
sx q[3];
rz(2.7494591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8490303) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(0.070177468) q[2];
rz(-2.3146546) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-0.58644811) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6475911) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-0.31627396) q[0];
rz(-2.0896185) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(-2.0064328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0907744) q[0];
sx q[0];
rz(-0.43681112) q[0];
sx q[0];
rz(0.45849053) q[0];
rz(-pi) q[1];
rz(-3.1124434) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(2.0707891) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75505776) q[1];
sx q[1];
rz(-2.6549576) q[1];
sx q[1];
rz(-1.4450032) q[1];
rz(-pi) q[2];
x q[2];
rz(0.023852392) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(-2.9305487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(-2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233376) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(-2.7897575) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(2.9454254) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10712121) q[0];
sx q[0];
rz(-2.2959318) q[0];
sx q[0];
rz(-2.3496773) q[0];
x q[1];
rz(0.050373366) q[2];
sx q[2];
rz(-1.435624) q[2];
sx q[2];
rz(2.2809575) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0386104) q[1];
sx q[1];
rz(-2.3048232) q[1];
sx q[1];
rz(3.1208913) q[1];
rz(-pi) q[2];
rz(-2.2979126) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(-2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(3.0604559) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(1.0749764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022973013) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(-2.9076004) q[2];
sx q[2];
rz(-1.6178314) q[2];
sx q[2];
rz(-2.9754054) q[2];
rz(-3.0951981) q[3];
sx q[3];
rz(-1.8002602) q[3];
sx q[3];
rz(0.044063448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
