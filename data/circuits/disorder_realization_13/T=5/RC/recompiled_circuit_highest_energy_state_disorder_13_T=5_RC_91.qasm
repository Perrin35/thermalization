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
rz(-2.9450671) q[0];
sx q[0];
rz(-2.3152469) q[0];
sx q[0];
rz(2.2235121) q[0];
rz(3.0300568) q[1];
sx q[1];
rz(-1.7938951) q[1];
sx q[1];
rz(-1.5674051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.471431) q[0];
sx q[0];
rz(-1.5613371) q[0];
sx q[0];
rz(-0.0034250101) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5882173) q[2];
sx q[2];
rz(-0.63734431) q[2];
sx q[2];
rz(3.1269249) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5985377) q[1];
sx q[1];
rz(-0.32078136) q[1];
sx q[1];
rz(-1.9108652) q[1];
rz(2.918675) q[3];
sx q[3];
rz(-0.50361982) q[3];
sx q[3];
rz(-0.38583392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7142746) q[2];
sx q[2];
rz(-1.6552507) q[2];
sx q[2];
rz(1.2662668) q[2];
rz(1.0424987) q[3];
sx q[3];
rz(-1.6008585) q[3];
sx q[3];
rz(1.7460167) q[3];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.105044) q[0];
sx q[0];
rz(-0.62750134) q[0];
sx q[0];
rz(-3.0446766) q[0];
rz(2.3333683) q[1];
sx q[1];
rz(-0.40882912) q[1];
sx q[1];
rz(1.2867297) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6280178) q[0];
sx q[0];
rz(-0.45125181) q[0];
sx q[0];
rz(1.4892764) q[0];
rz(0.98640826) q[2];
sx q[2];
rz(-2.2041568) q[2];
sx q[2];
rz(1.4937166) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5401973) q[1];
sx q[1];
rz(-2.2608065) q[1];
sx q[1];
rz(-0.50443919) q[1];
x q[2];
rz(0.80974354) q[3];
sx q[3];
rz(-2.6698723) q[3];
sx q[3];
rz(-0.37744409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6404932) q[2];
sx q[2];
rz(-0.88397908) q[2];
sx q[2];
rz(-0.38197771) q[2];
rz(-2.6169418) q[3];
sx q[3];
rz(-1.7347521) q[3];
sx q[3];
rz(0.14373556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76661888) q[0];
sx q[0];
rz(-1.5559649) q[0];
sx q[0];
rz(2.4554456) q[0];
rz(-2.0740017) q[1];
sx q[1];
rz(-2.144404) q[1];
sx q[1];
rz(0.080726191) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19182675) q[0];
sx q[0];
rz(-1.6097665) q[0];
sx q[0];
rz(-1.3471589) q[0];
x q[1];
rz(-2.2171634) q[2];
sx q[2];
rz(-1.0253128) q[2];
sx q[2];
rz(1.9540862) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1898274) q[1];
sx q[1];
rz(-1.3406585) q[1];
sx q[1];
rz(1.5946424) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3537797) q[3];
sx q[3];
rz(-1.9160877) q[3];
sx q[3];
rz(2.6097176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2637691) q[2];
sx q[2];
rz(-2.4407385) q[2];
sx q[2];
rz(2.7027255) q[2];
rz(1.8435439) q[3];
sx q[3];
rz(-0.38709199) q[3];
sx q[3];
rz(-0.41518655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91442672) q[0];
sx q[0];
rz(-0.54410797) q[0];
sx q[0];
rz(0.67657226) q[0];
rz(-0.1380955) q[1];
sx q[1];
rz(-1.0905677) q[1];
sx q[1];
rz(0.20060435) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56708589) q[0];
sx q[0];
rz(-0.73154035) q[0];
sx q[0];
rz(-1.0163496) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3820962) q[2];
sx q[2];
rz(-0.9498792) q[2];
sx q[2];
rz(-0.046227235) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.323206) q[1];
sx q[1];
rz(-0.96508677) q[1];
sx q[1];
rz(1.4097286) q[1];
rz(-pi) q[2];
rz(-0.225876) q[3];
sx q[3];
rz(-2.568733) q[3];
sx q[3];
rz(-1.7879888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5033919) q[2];
sx q[2];
rz(-1.2752504) q[2];
sx q[2];
rz(1.4208043) q[2];
rz(3.0270789) q[3];
sx q[3];
rz(-1.4736466) q[3];
sx q[3];
rz(0.99944559) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81922174) q[0];
sx q[0];
rz(-2.350816) q[0];
sx q[0];
rz(-2.1694699) q[0];
rz(-0.32360336) q[1];
sx q[1];
rz(-1.4515896) q[1];
sx q[1];
rz(1.1824898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33859461) q[0];
sx q[0];
rz(-1.3756244) q[0];
sx q[0];
rz(0.040935733) q[0];
rz(-2.5713872) q[2];
sx q[2];
rz(-1.5758762) q[2];
sx q[2];
rz(-0.069381086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3283493) q[1];
sx q[1];
rz(-2.9651838) q[1];
sx q[1];
rz(1.117068) q[1];
rz(1.0504405) q[3];
sx q[3];
rz(-2.0791868) q[3];
sx q[3];
rz(2.6843478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16367308) q[2];
sx q[2];
rz(-1.7719496) q[2];
sx q[2];
rz(-2.6046806) q[2];
rz(0.72530693) q[3];
sx q[3];
rz(-3.091843) q[3];
sx q[3];
rz(2.3427826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2742915) q[0];
sx q[0];
rz(-2.0766356) q[0];
sx q[0];
rz(1.4428447) q[0];
rz(0.2105712) q[1];
sx q[1];
rz(-1.6981533) q[1];
sx q[1];
rz(-2.494536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4763586) q[0];
sx q[0];
rz(-1.826735) q[0];
sx q[0];
rz(0.38689918) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9672616) q[2];
sx q[2];
rz(-2.5781879) q[2];
sx q[2];
rz(-2.8772815) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0004859) q[1];
sx q[1];
rz(-0.55505156) q[1];
sx q[1];
rz(0.6935814) q[1];
rz(-1.296455) q[3];
sx q[3];
rz(-1.8963433) q[3];
sx q[3];
rz(-0.59301585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1232542) q[2];
sx q[2];
rz(-2.047057) q[2];
sx q[2];
rz(1.9617762) q[2];
rz(1.2849464) q[3];
sx q[3];
rz(-2.732087) q[3];
sx q[3];
rz(3.0783317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1061851) q[0];
sx q[0];
rz(-1.4466865) q[0];
sx q[0];
rz(-2.2440198) q[0];
rz(1.3342185) q[1];
sx q[1];
rz(-1.370627) q[1];
sx q[1];
rz(-0.99475494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.008381) q[0];
sx q[0];
rz(-1.8803839) q[0];
sx q[0];
rz(-2.3214843) q[0];
x q[1];
rz(-2.686196) q[2];
sx q[2];
rz(-1.0173305) q[2];
sx q[2];
rz(2.0023416) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2487476) q[1];
sx q[1];
rz(-1.1214646) q[1];
sx q[1];
rz(2.8161418) q[1];
x q[2];
rz(1.654387) q[3];
sx q[3];
rz(-0.28161848) q[3];
sx q[3];
rz(2.6476988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.083954088) q[2];
sx q[2];
rz(-1.8123764) q[2];
sx q[2];
rz(-0.29339054) q[2];
rz(3.0883664) q[3];
sx q[3];
rz(-2.2727727) q[3];
sx q[3];
rz(1.4330385) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5720125) q[0];
sx q[0];
rz(-1.7520289) q[0];
sx q[0];
rz(-2.8386175) q[0];
rz(1.4888034) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(-0.66910076) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55941033) q[0];
sx q[0];
rz(-0.56920496) q[0];
sx q[0];
rz(2.4962382) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2190518) q[2];
sx q[2];
rz(-1.3430077) q[2];
sx q[2];
rz(-2.5643189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2337906) q[1];
sx q[1];
rz(-1.0668584) q[1];
sx q[1];
rz(-2.6556117) q[1];
rz(-pi) q[2];
rz(-1.331393) q[3];
sx q[3];
rz(-1.1576011) q[3];
sx q[3];
rz(-1.5904782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1423433) q[2];
sx q[2];
rz(-1.3331022) q[2];
sx q[2];
rz(-0.86714253) q[2];
rz(-2.0182746) q[3];
sx q[3];
rz(-0.85223782) q[3];
sx q[3];
rz(1.999202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2808696) q[0];
sx q[0];
rz(-0.46361247) q[0];
sx q[0];
rz(0.11775693) q[0];
rz(0.87751687) q[1];
sx q[1];
rz(-2.120647) q[1];
sx q[1];
rz(2.1868475) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68761008) q[0];
sx q[0];
rz(-1.1907401) q[0];
sx q[0];
rz(-0.65681547) q[0];
rz(0.81127848) q[2];
sx q[2];
rz(-2.5015321) q[2];
sx q[2];
rz(1.0803573) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6462209) q[1];
sx q[1];
rz(-0.89447656) q[1];
sx q[1];
rz(0.27797525) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3225962) q[3];
sx q[3];
rz(-2.1124509) q[3];
sx q[3];
rz(1.1198695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8942326) q[2];
sx q[2];
rz(-0.33984137) q[2];
sx q[2];
rz(1.6575238) q[2];
rz(-1.6190716) q[3];
sx q[3];
rz(-1.7584453) q[3];
sx q[3];
rz(1.1744261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-2.5840983) q[0];
sx q[0];
rz(-1.431594) q[0];
sx q[0];
rz(2.3884921) q[0];
rz(-1.151471) q[1];
sx q[1];
rz(-0.52422062) q[1];
sx q[1];
rz(1.6023191) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70645151) q[0];
sx q[0];
rz(-1.1223842) q[0];
sx q[0];
rz(0.20332341) q[0];
x q[1];
rz(-0.87290092) q[2];
sx q[2];
rz(-1.97556) q[2];
sx q[2];
rz(-1.0412316) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11048296) q[1];
sx q[1];
rz(-1.3785161) q[1];
sx q[1];
rz(-0.39327217) q[1];
rz(-pi) q[2];
rz(-1.3063823) q[3];
sx q[3];
rz(-2.2087277) q[3];
sx q[3];
rz(-0.3141981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3573542) q[2];
sx q[2];
rz(-2.9604762) q[2];
sx q[2];
rz(2.6757346) q[2];
rz(1.0957796) q[3];
sx q[3];
rz(-2.1491137) q[3];
sx q[3];
rz(-2.5994658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096238484) q[0];
sx q[0];
rz(-1.5753373) q[0];
sx q[0];
rz(-1.5684431) q[0];
rz(-1.0678328) q[1];
sx q[1];
rz(-2.7012431) q[1];
sx q[1];
rz(2.3213097) q[1];
rz(-2.4196923) q[2];
sx q[2];
rz(-2.7963287) q[2];
sx q[2];
rz(-2.1523274) q[2];
rz(-0.54084861) q[3];
sx q[3];
rz(-1.114308) q[3];
sx q[3];
rz(-1.4642117) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
