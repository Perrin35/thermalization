OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(2.7913845) q[0];
sx q[0];
rz(12.933001) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(-1.6860513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3722056) q[0];
sx q[0];
rz(-1.9908449) q[0];
sx q[0];
rz(2.7647892) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8800814) q[2];
sx q[2];
rz(-2.4158084) q[2];
sx q[2];
rz(-0.45553614) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9283596) q[1];
sx q[1];
rz(-1.8395437) q[1];
sx q[1];
rz(-2.1409722) q[1];
rz(-pi) q[2];
rz(-0.22110181) q[3];
sx q[3];
rz(-3.0924774) q[3];
sx q[3];
rz(2.5887095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(2.74995) q[2];
rz(0.13970217) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(-0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249107) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.2044027) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5088168) q[0];
sx q[0];
rz(-2.5460089) q[0];
sx q[0];
rz(-2.3217208) q[0];
rz(-pi) q[1];
rz(-0.18560974) q[2];
sx q[2];
rz(-2.2679272) q[2];
sx q[2];
rz(-0.44167659) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7960297) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(1.4644535) q[1];
rz(1.9235839) q[3];
sx q[3];
rz(-0.70397607) q[3];
sx q[3];
rz(2.9670027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5045972) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(-1.9223928) q[2];
rz(-2.7870264) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59042674) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(-3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(1.191167) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5844849) q[0];
sx q[0];
rz(-0.35042052) q[0];
sx q[0];
rz(2.0273897) q[0];
rz(-0.28352719) q[2];
sx q[2];
rz(-1.3361738) q[2];
sx q[2];
rz(-1.4153751) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35347963) q[1];
sx q[1];
rz(-2.8444926) q[1];
sx q[1];
rz(2.3658845) q[1];
rz(-2.5997945) q[3];
sx q[3];
rz(-1.9199315) q[3];
sx q[3];
rz(-0.19402129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1742192) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(-2.1276786) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(1.0323662) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8106666) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(-2.9372835) q[0];
rz(-1.3775685) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(0.92299443) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5343691) q[0];
sx q[0];
rz(-1.9019706) q[0];
sx q[0];
rz(1.2481199) q[0];
rz(-pi) q[1];
rz(1.5590018) q[2];
sx q[2];
rz(-1.6147385) q[2];
sx q[2];
rz(-0.61461385) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2642306) q[1];
sx q[1];
rz(-1.1738395) q[1];
sx q[1];
rz(-3.0287659) q[1];
rz(-pi) q[2];
rz(1.5444078) q[3];
sx q[3];
rz(-0.70453405) q[3];
sx q[3];
rz(-1.1772616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(2.6320809) q[2];
rz(-0.64309684) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5287857) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(0.91941419) q[0];
rz(-2.1557504) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(-1.3607508) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1795792) q[0];
sx q[0];
rz(-2.8498581) q[0];
sx q[0];
rz(-0.71266642) q[0];
rz(-pi) q[1];
rz(-2.1941357) q[2];
sx q[2];
rz(-0.91295487) q[2];
sx q[2];
rz(2.5156977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67140019) q[1];
sx q[1];
rz(-0.50506401) q[1];
sx q[1];
rz(-0.43882521) q[1];
rz(1.4578044) q[3];
sx q[3];
rz(-0.60243536) q[3];
sx q[3];
rz(-0.065954176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3551066) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(-0.11486593) q[2];
rz(-2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961261) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(-1.8967569) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(0.27522603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305785) q[0];
sx q[0];
rz(-2.0654581) q[0];
sx q[0];
rz(0.37822322) q[0];
rz(-2.2116824) q[2];
sx q[2];
rz(-1.9800131) q[2];
sx q[2];
rz(0.90604679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.59906193) q[1];
sx q[1];
rz(-1.9179357) q[1];
sx q[1];
rz(-2.6229726) q[1];
rz(-pi) q[2];
rz(-0.082645881) q[3];
sx q[3];
rz(-0.9517037) q[3];
sx q[3];
rz(-1.9649399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.747921) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(-2.6829524) q[2];
rz(0.7115055) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(-2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28073072) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(1.3635427) q[0];
rz(1.7215615) q[1];
sx q[1];
rz(-2.6224711) q[1];
sx q[1];
rz(0.71969676) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62455432) q[0];
sx q[0];
rz(-1.8554243) q[0];
sx q[0];
rz(-2.4863003) q[0];
rz(-2.4160956) q[2];
sx q[2];
rz(-0.91425397) q[2];
sx q[2];
rz(0.072349116) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.91832671) q[1];
sx q[1];
rz(-1.0263138) q[1];
sx q[1];
rz(1.5189927) q[1];
rz(-1.8411438) q[3];
sx q[3];
rz(-2.2336604) q[3];
sx q[3];
rz(0.28966749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.482243) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-0.68816319) q[2];
rz(-3.0996389) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(2.8614614) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(-2.9550609) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(1.6465181) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4743487) q[0];
sx q[0];
rz(-1.7843887) q[0];
sx q[0];
rz(2.7937907) q[0];
x q[1];
rz(1.375884) q[2];
sx q[2];
rz(-2.2576828) q[2];
sx q[2];
rz(-1.9538823) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9547792) q[1];
sx q[1];
rz(-0.28007945) q[1];
sx q[1];
rz(2.7571452) q[1];
x q[2];
rz(-2.2349615) q[3];
sx q[3];
rz(-2.2615221) q[3];
sx q[3];
rz(-3.0674792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8994393) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(3.126826) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.8321962) q[3];
sx q[3];
rz(1.7285255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(2.263608) q[0];
rz(0.19628482) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(1.7810129) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6920647) q[0];
sx q[0];
rz(-0.99943107) q[0];
sx q[0];
rz(-0.044762386) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94366818) q[2];
sx q[2];
rz(-1.5329648) q[2];
sx q[2];
rz(0.19247069) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.18371836) q[1];
sx q[1];
rz(-2.6596333) q[1];
sx q[1];
rz(2.3493489) q[1];
x q[2];
rz(-1.1514879) q[3];
sx q[3];
rz(-1.7125704) q[3];
sx q[3];
rz(-1.7820953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62375162) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(0.8927792) q[2];
rz(-3.1023074) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(-0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2980625) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(0.046534006) q[0];
rz(-0.90905601) q[1];
sx q[1];
rz(-0.87288705) q[1];
sx q[1];
rz(-2.4216901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35346183) q[0];
sx q[0];
rz(-2.8849368) q[0];
sx q[0];
rz(-1.8887595) q[0];
rz(-3.1245329) q[2];
sx q[2];
rz(-2.8327201) q[2];
sx q[2];
rz(2.942254) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9626999) q[1];
sx q[1];
rz(-1.3239412) q[1];
sx q[1];
rz(-0.015951338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4107237) q[3];
sx q[3];
rz(-2.1248397) q[3];
sx q[3];
rz(2.4027367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-3.1151248) q[2];
rz(1.3039533) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(-1.8855689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532886) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-2.9251255) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(2.4344865) q[2];
sx q[2];
rz(-1.8137992) q[2];
sx q[2];
rz(0.9546311) q[2];
rz(-1.6197694) q[3];
sx q[3];
rz(-2.5544142) q[3];
sx q[3];
rz(2.1094473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
