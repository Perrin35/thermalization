OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(-2.1043632) q[0];
sx q[0];
rz(-0.35559911) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(-1.27966) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58334938) q[0];
sx q[0];
rz(-1.8177176) q[0];
sx q[0];
rz(2.9980744) q[0];
rz(-0.46977629) q[2];
sx q[2];
rz(-2.6758286) q[2];
sx q[2];
rz(-0.85675203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76525926) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(0.63912649) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7198556) q[3];
sx q[3];
rz(-0.30656439) q[3];
sx q[3];
rz(-2.7048064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0063643) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(2.485086) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0083369) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(2.6541236) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23232846) q[0];
sx q[0];
rz(-3.0380913) q[0];
sx q[0];
rz(-0.89719015) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42178085) q[2];
sx q[2];
rz(-1.1676719) q[2];
sx q[2];
rz(-3.1096854) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.65709719) q[1];
sx q[1];
rz(-2.8946981) q[1];
sx q[1];
rz(0.8685649) q[1];
rz(0.08667605) q[3];
sx q[3];
rz(-1.254734) q[3];
sx q[3];
rz(2.812127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(-2.7462192) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.19668002) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(2.9188459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4653141) q[0];
sx q[0];
rz(-1.197581) q[0];
sx q[0];
rz(0.42155427) q[0];
x q[1];
rz(-1.1995302) q[2];
sx q[2];
rz(-1.3604593) q[2];
sx q[2];
rz(2.631275) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5388515) q[1];
sx q[1];
rz(-0.95898421) q[1];
sx q[1];
rz(-1.4879423) q[1];
rz(-pi) q[2];
rz(-1.907903) q[3];
sx q[3];
rz(-1.8030231) q[3];
sx q[3];
rz(3.0984578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2591851) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(1.9474691) q[2];
rz(-2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38917437) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46566761) q[0];
sx q[0];
rz(-1.9786069) q[0];
sx q[0];
rz(-0.84665914) q[0];
rz(-pi) q[1];
rz(1.4843066) q[2];
sx q[2];
rz(-1.2109586) q[2];
sx q[2];
rz(-0.53777018) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3355545) q[1];
sx q[1];
rz(-2.028855) q[1];
sx q[1];
rz(2.4469417) q[1];
rz(2.782015) q[3];
sx q[3];
rz(-1.3782856) q[3];
sx q[3];
rz(-2.0457091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(-1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(1.4047594) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-0.4531025) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3081007) q[0];
sx q[0];
rz(-1.5335324) q[0];
sx q[0];
rz(1.4971855) q[0];
x q[1];
rz(1.808851) q[2];
sx q[2];
rz(-0.14698725) q[2];
sx q[2];
rz(2.5631995) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52757712) q[1];
sx q[1];
rz(-1.7909389) q[1];
sx q[1];
rz(0.084053587) q[1];
rz(0.971332) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(0.99861162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0129464) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.5117234) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033427514) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-0.23813716) q[0];
rz(2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.5135117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1171922) q[0];
sx q[0];
rz(-1.8937366) q[0];
sx q[0];
rz(-0.49844235) q[0];
rz(-2.020535) q[2];
sx q[2];
rz(-1.626484) q[2];
sx q[2];
rz(-2.4186717) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39735079) q[1];
sx q[1];
rz(-2.2480818) q[1];
sx q[1];
rz(0.97126295) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38576491) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(2.5667131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(-3.011076) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(2.9706764) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1259574) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(2.7440199) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.9546753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884739) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(-2.8209646) q[0];
rz(-pi) q[1];
rz(3.0824532) q[2];
sx q[2];
rz(-1.5999609) q[2];
sx q[2];
rz(-1.189122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6656875) q[1];
sx q[1];
rz(-2.396282) q[1];
sx q[1];
rz(1.7687294) q[1];
rz(-2.7508221) q[3];
sx q[3];
rz(-0.6676538) q[3];
sx q[3];
rz(-2.8323176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26178965) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(-1.6342182) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-2.2858455) q[0];
rz(0.021727173) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-1.0303248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1925826) q[0];
sx q[0];
rz(-0.858925) q[0];
sx q[0];
rz(-2.683995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4872929) q[2];
sx q[2];
rz(-2.6132772) q[2];
sx q[2];
rz(-2.5022142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1460261) q[1];
sx q[1];
rz(-1.7221754) q[1];
sx q[1];
rz(-1.8712908) q[1];
rz(-pi) q[2];
rz(0.97708948) q[3];
sx q[3];
rz(-1.1247375) q[3];
sx q[3];
rz(-2.2462728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29256233) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(-0.070177468) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(2.8253187) q[0];
rz(-2.0896185) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(-2.0064328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508182) q[0];
sx q[0];
rz(-0.43681112) q[0];
sx q[0];
rz(0.45849053) q[0];
rz(-pi) q[1];
rz(-0.029149292) q[2];
sx q[2];
rz(-1.8910732) q[2];
sx q[2];
rz(2.0707891) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75505776) q[1];
sx q[1];
rz(-0.48663501) q[1];
sx q[1];
rz(-1.6965894) q[1];
x q[2];
rz(1.6039861) q[3];
sx q[3];
rz(-0.6233218) q[3];
sx q[3];
rz(-2.8896796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(1.5448145) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(-2.9454254) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0553186) q[0];
sx q[0];
rz(-1.0090758) q[0];
sx q[0];
rz(2.4713211) q[0];
x q[1];
rz(1.9253795) q[2];
sx q[2];
rz(-2.9973929) q[2];
sx q[2];
rz(1.9229638) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6876467) q[1];
sx q[1];
rz(-1.5554264) q[1];
sx q[1];
rz(2.3049298) q[1];
rz(0.65532834) q[3];
sx q[3];
rz(-0.95643759) q[3];
sx q[3];
rz(2.4650246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(-1.0749764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(1.5224456) q[2];
sx q[2];
rz(-1.804525) q[2];
sx q[2];
rz(-1.3934025) q[2];
rz(-1.7668016) q[3];
sx q[3];
rz(-2.9075665) q[3];
sx q[3];
rz(-2.8961765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
