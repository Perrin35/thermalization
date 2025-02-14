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
rz(0.19658495) q[0];
sx q[0];
rz(6.7253334) q[0];
sx q[0];
rz(10.286177) q[0];
rz(1.5098894) q[1];
sx q[1];
rz(-1.8169401) q[1];
sx q[1];
rz(3.1060001) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0033158) q[0];
sx q[0];
rz(-2.7502923) q[0];
sx q[0];
rz(-2.6075493) q[0];
rz(-0.30649779) q[2];
sx q[2];
rz(-1.4671637) q[2];
sx q[2];
rz(0.95970861) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0390748) q[1];
sx q[1];
rz(-2.2968484) q[1];
sx q[1];
rz(-1.5229358) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.54856) q[3];
sx q[3];
rz(-2.4614868) q[3];
sx q[3];
rz(2.3176127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55277905) q[2];
sx q[2];
rz(-1.8401044) q[2];
sx q[2];
rz(1.0944875) q[2];
rz(-1.6257446) q[3];
sx q[3];
rz(-0.87259126) q[3];
sx q[3];
rz(-0.14357963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.0482408) q[0];
sx q[0];
rz(-0.99026647) q[0];
sx q[0];
rz(-0.37407237) q[0];
rz(-0.85809842) q[1];
sx q[1];
rz(-2.2007807) q[1];
sx q[1];
rz(-2.0657952) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825622) q[0];
sx q[0];
rz(-0.1452862) q[0];
sx q[0];
rz(-0.7913211) q[0];
rz(-0.38041842) q[2];
sx q[2];
rz(-2.5606692) q[2];
sx q[2];
rz(2.7564633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3871284) q[1];
sx q[1];
rz(-2.7242047) q[1];
sx q[1];
rz(-0.2902969) q[1];
rz(-pi) q[2];
rz(2.2298563) q[3];
sx q[3];
rz(-1.0056595) q[3];
sx q[3];
rz(-2.679058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6263803) q[2];
sx q[2];
rz(-0.43312803) q[2];
sx q[2];
rz(1.0832146) q[2];
rz(1.2927239) q[3];
sx q[3];
rz(-2.0079398) q[3];
sx q[3];
rz(2.2129272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2556297) q[0];
sx q[0];
rz(-0.77115458) q[0];
sx q[0];
rz(0.5249002) q[0];
rz(1.4595754) q[1];
sx q[1];
rz(-0.73760215) q[1];
sx q[1];
rz(2.9958013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4519212) q[0];
sx q[0];
rz(-1.6761685) q[0];
sx q[0];
rz(1.4793878) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86051382) q[2];
sx q[2];
rz(-2.4053221) q[2];
sx q[2];
rz(-2.3585573) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1555712) q[1];
sx q[1];
rz(-2.0029481) q[1];
sx q[1];
rz(1.4714824) q[1];
rz(-0.84050982) q[3];
sx q[3];
rz(-2.0183992) q[3];
sx q[3];
rz(2.2882035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24474457) q[2];
sx q[2];
rz(-1.1620099) q[2];
sx q[2];
rz(2.9761918) q[2];
rz(2.2798955) q[3];
sx q[3];
rz(-2.4498037) q[3];
sx q[3];
rz(2.944788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8621314) q[0];
sx q[0];
rz(-0.52424163) q[0];
sx q[0];
rz(-0.72823802) q[0];
rz(-0.32314745) q[1];
sx q[1];
rz(-1.0341945) q[1];
sx q[1];
rz(-1.039215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4711888) q[0];
sx q[0];
rz(-1.146933) q[0];
sx q[0];
rz(-2.7039862) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2065998) q[2];
sx q[2];
rz(-1.3424917) q[2];
sx q[2];
rz(1.993597) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8700712) q[1];
sx q[1];
rz(-1.0876552) q[1];
sx q[1];
rz(1.6752233) q[1];
rz(-pi) q[2];
rz(0.69833243) q[3];
sx q[3];
rz(-2.192333) q[3];
sx q[3];
rz(0.96635287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6497961) q[2];
sx q[2];
rz(-1.0504664) q[2];
sx q[2];
rz(-0.51898471) q[2];
rz(-1.0745878) q[3];
sx q[3];
rz(-1.855775) q[3];
sx q[3];
rz(0.49526596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.062926) q[0];
sx q[0];
rz(-0.92200297) q[0];
sx q[0];
rz(2.9700188) q[0];
rz(2.0473139) q[1];
sx q[1];
rz(-1.8405874) q[1];
sx q[1];
rz(-2.9069854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0594992) q[0];
sx q[0];
rz(-1.4429868) q[0];
sx q[0];
rz(1.2907842) q[0];
rz(-2.634356) q[2];
sx q[2];
rz(-1.5248858) q[2];
sx q[2];
rz(3.1243589) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7947096) q[1];
sx q[1];
rz(-0.88635072) q[1];
sx q[1];
rz(0.80839949) q[1];
x q[2];
rz(0.42886491) q[3];
sx q[3];
rz(-2.3842616) q[3];
sx q[3];
rz(2.0067818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2285063) q[2];
sx q[2];
rz(-1.9501016) q[2];
sx q[2];
rz(-1.4027493) q[2];
rz(2.5168822) q[3];
sx q[3];
rz(-2.1731845) q[3];
sx q[3];
rz(1.3431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45873555) q[0];
sx q[0];
rz(-3.1338437) q[0];
sx q[0];
rz(-0.32753456) q[0];
rz(-0.028118357) q[1];
sx q[1];
rz(-1.1437806) q[1];
sx q[1];
rz(0.99172529) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0077181) q[0];
sx q[0];
rz(-1.0622866) q[0];
sx q[0];
rz(0.74061142) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4430094) q[2];
sx q[2];
rz(-0.24925079) q[2];
sx q[2];
rz(-0.64068782) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27665972) q[1];
sx q[1];
rz(-0.75705479) q[1];
sx q[1];
rz(-1.6318342) q[1];
rz(-pi) q[2];
rz(2.670861) q[3];
sx q[3];
rz(-1.5464029) q[3];
sx q[3];
rz(2.7324146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.907054) q[2];
sx q[2];
rz(-1.3776642) q[2];
sx q[2];
rz(1.482359) q[2];
rz(2.0756857) q[3];
sx q[3];
rz(-1.9612471) q[3];
sx q[3];
rz(-1.0970241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7557573) q[0];
sx q[0];
rz(-0.11180728) q[0];
sx q[0];
rz(2.8461611) q[0];
rz(2.9040728) q[1];
sx q[1];
rz(-0.93370456) q[1];
sx q[1];
rz(0.74737731) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4656671) q[0];
sx q[0];
rz(-2.3479778) q[0];
sx q[0];
rz(2.3966952) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5146444) q[2];
sx q[2];
rz(-2.0184506) q[2];
sx q[2];
rz(-0.7426151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57495368) q[1];
sx q[1];
rz(-0.63640928) q[1];
sx q[1];
rz(-0.91729911) q[1];
rz(-pi) q[2];
rz(-2.5980972) q[3];
sx q[3];
rz(-1.8955909) q[3];
sx q[3];
rz(2.3372758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6682917) q[2];
sx q[2];
rz(-1.95582) q[2];
sx q[2];
rz(2.8866923) q[2];
rz(-2.9292246) q[3];
sx q[3];
rz(-0.86447159) q[3];
sx q[3];
rz(-0.00096850639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.85106987) q[0];
sx q[0];
rz(-1.2238598) q[0];
sx q[0];
rz(3.0176924) q[0];
rz(-2.2663785) q[1];
sx q[1];
rz(-2.3085935) q[1];
sx q[1];
rz(-2.5828054) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2921857) q[0];
sx q[0];
rz(-1.8208139) q[0];
sx q[0];
rz(2.1110737) q[0];
rz(2.8577639) q[2];
sx q[2];
rz(-2.0528194) q[2];
sx q[2];
rz(-1.1555156) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3271823) q[1];
sx q[1];
rz(-0.60403999) q[1];
sx q[1];
rz(-2.6534897) q[1];
x q[2];
rz(2.111192) q[3];
sx q[3];
rz(-2.4321788) q[3];
sx q[3];
rz(-1.4340925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.045804068) q[2];
sx q[2];
rz(-2.6498821) q[2];
sx q[2];
rz(2.4208505) q[2];
rz(0.36302429) q[3];
sx q[3];
rz(-1.9178773) q[3];
sx q[3];
rz(-1.5117517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1849798) q[0];
sx q[0];
rz(-2.9473801) q[0];
sx q[0];
rz(-0.089056253) q[0];
rz(-2.7748499) q[1];
sx q[1];
rz(-2.2373503) q[1];
sx q[1];
rz(1.6820224) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7686667) q[0];
sx q[0];
rz(-2.4690876) q[0];
sx q[0];
rz(1.9422533) q[0];
rz(-pi) q[1];
rz(2.4602997) q[2];
sx q[2];
rz(-0.80901399) q[2];
sx q[2];
rz(1.1920649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1398581) q[1];
sx q[1];
rz(-2.325104) q[1];
sx q[1];
rz(1.0553318) q[1];
x q[2];
rz(-2.0085387) q[3];
sx q[3];
rz(-1.0171842) q[3];
sx q[3];
rz(2.5598524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9350932) q[2];
sx q[2];
rz(-1.3672028) q[2];
sx q[2];
rz(1.2500786) q[2];
rz(-1.9153204) q[3];
sx q[3];
rz(-0.63392249) q[3];
sx q[3];
rz(-0.63898501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117689) q[0];
sx q[0];
rz(-2.0084232) q[0];
sx q[0];
rz(2.0015707) q[0];
rz(2.5584768) q[1];
sx q[1];
rz(-1.4864328) q[1];
sx q[1];
rz(-0.66782943) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44396469) q[0];
sx q[0];
rz(-2.2831342) q[0];
sx q[0];
rz(1.1348073) q[0];
rz(2.2726353) q[2];
sx q[2];
rz(-1.7545106) q[2];
sx q[2];
rz(-2.3784203) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.30582593) q[1];
sx q[1];
rz(-1.139031) q[1];
sx q[1];
rz(2.9271896) q[1];
x q[2];
rz(-1.1125426) q[3];
sx q[3];
rz(-2.309121) q[3];
sx q[3];
rz(2.3979681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2976611) q[2];
sx q[2];
rz(-2.6466978) q[2];
sx q[2];
rz(0.75221357) q[2];
rz(-0.080605896) q[3];
sx q[3];
rz(-1.3496496) q[3];
sx q[3];
rz(2.5940671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3362296) q[0];
sx q[0];
rz(-1.5267876) q[0];
sx q[0];
rz(3.0841854) q[0];
rz(-2.2930131) q[1];
sx q[1];
rz(-0.72863693) q[1];
sx q[1];
rz(-1.4562664) q[1];
rz(0.91114974) q[2];
sx q[2];
rz(-1.0187889) q[2];
sx q[2];
rz(2.7762882) q[2];
rz(-0.46075321) q[3];
sx q[3];
rz(-1.2450153) q[3];
sx q[3];
rz(1.0871441) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
