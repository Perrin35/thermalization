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
rz(-1.8347972) q[0];
sx q[0];
rz(-2.2936294) q[0];
sx q[0];
rz(-3.1041978) q[0];
rz(2.1283863) q[1];
sx q[1];
rz(-2.6599045) q[1];
sx q[1];
rz(-1.4451292) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73611063) q[0];
sx q[0];
rz(-0.47022143) q[0];
sx q[0];
rz(0.41023631) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6736534) q[2];
sx q[2];
rz(-1.6409931) q[2];
sx q[2];
rz(1.6818893) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7788275) q[1];
sx q[1];
rz(-2.3794075) q[1];
sx q[1];
rz(2.306421) q[1];
rz(-pi) q[2];
rz(2.1358498) q[3];
sx q[3];
rz(-1.6673207) q[3];
sx q[3];
rz(0.18894503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41363132) q[2];
sx q[2];
rz(-0.093158826) q[2];
sx q[2];
rz(1.6348582) q[2];
rz(0.19351752) q[3];
sx q[3];
rz(-2.3573124) q[3];
sx q[3];
rz(0.94582742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042260878) q[0];
sx q[0];
rz(-1.7758545) q[0];
sx q[0];
rz(2.9625764) q[0];
rz(-1.5366813) q[1];
sx q[1];
rz(-2.8004526) q[1];
sx q[1];
rz(-1.5515597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1354438) q[0];
sx q[0];
rz(-1.8601226) q[0];
sx q[0];
rz(3.0399714) q[0];
rz(-pi) q[1];
rz(0.29921542) q[2];
sx q[2];
rz(-0.73198527) q[2];
sx q[2];
rz(-0.74886403) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4578456) q[1];
sx q[1];
rz(-0.48626712) q[1];
sx q[1];
rz(-2.1621428) q[1];
rz(-pi) q[2];
rz(-2.670774) q[3];
sx q[3];
rz(-2.0662466) q[3];
sx q[3];
rz(-1.214787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.67148525) q[2];
sx q[2];
rz(-1.6652197) q[2];
sx q[2];
rz(0.022484953) q[2];
rz(0.89618987) q[3];
sx q[3];
rz(-0.78138566) q[3];
sx q[3];
rz(1.6270858) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3317868) q[0];
sx q[0];
rz(-2.9349194) q[0];
sx q[0];
rz(-1.5326387) q[0];
rz(-1.1398075) q[1];
sx q[1];
rz(-2.8228357) q[1];
sx q[1];
rz(-0.87361139) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667359) q[0];
sx q[0];
rz(-0.92467659) q[0];
sx q[0];
rz(0.91888756) q[0];
rz(-pi) q[1];
rz(0.68364667) q[2];
sx q[2];
rz(-2.3693759) q[2];
sx q[2];
rz(-1.5600086) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4406179) q[1];
sx q[1];
rz(-1.7243885) q[1];
sx q[1];
rz(-0.48753341) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8779712) q[3];
sx q[3];
rz(-0.28214165) q[3];
sx q[3];
rz(0.96708114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56796271) q[2];
sx q[2];
rz(-0.45845389) q[2];
sx q[2];
rz(-0.47631329) q[2];
rz(-1.025398) q[3];
sx q[3];
rz(-1.3566596) q[3];
sx q[3];
rz(2.8477113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3697701) q[0];
sx q[0];
rz(-0.85404587) q[0];
sx q[0];
rz(2.5452132) q[0];
rz(-0.75309938) q[1];
sx q[1];
rz(-1.3558847) q[1];
sx q[1];
rz(0.3515884) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0344926) q[0];
sx q[0];
rz(-1.4395243) q[0];
sx q[0];
rz(2.6126325) q[0];
x q[1];
rz(0.88732052) q[2];
sx q[2];
rz(-1.1850238) q[2];
sx q[2];
rz(1.224347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56573136) q[1];
sx q[1];
rz(-1.8003776) q[1];
sx q[1];
rz(-1.9750392) q[1];
rz(-pi) q[2];
rz(-2.5078689) q[3];
sx q[3];
rz(-2.124064) q[3];
sx q[3];
rz(0.80851698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72254649) q[2];
sx q[2];
rz(-1.6435577) q[2];
sx q[2];
rz(0.92918116) q[2];
rz(-1.9123214) q[3];
sx q[3];
rz(-3.1049187) q[3];
sx q[3];
rz(1.8721972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96255985) q[0];
sx q[0];
rz(-2.5293009) q[0];
sx q[0];
rz(0.52781934) q[0];
rz(-2.5714696) q[1];
sx q[1];
rz(-2.5716883) q[1];
sx q[1];
rz(-0.39400563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6878789) q[0];
sx q[0];
rz(-1.7648932) q[0];
sx q[0];
rz(0.48783036) q[0];
x q[1];
rz(0.5929596) q[2];
sx q[2];
rz(-1.8793794) q[2];
sx q[2];
rz(1.2013931) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2886823) q[1];
sx q[1];
rz(-0.25215071) q[1];
sx q[1];
rz(2.0618477) q[1];
rz(-pi) q[2];
rz(-0.72368716) q[3];
sx q[3];
rz(-2.1877648) q[3];
sx q[3];
rz(-2.8599515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3923308) q[2];
sx q[2];
rz(-1.8529842) q[2];
sx q[2];
rz(1.1345081) q[2];
rz(-3.115861) q[3];
sx q[3];
rz(-1.0545571) q[3];
sx q[3];
rz(2.499685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57539097) q[0];
sx q[0];
rz(-1.8113149) q[0];
sx q[0];
rz(0.45912418) q[0];
rz(0.8210012) q[1];
sx q[1];
rz(-0.88290015) q[1];
sx q[1];
rz(-2.6301036) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5465062) q[0];
sx q[0];
rz(-1.0985785) q[0];
sx q[0];
rz(-2.8123463) q[0];
rz(-pi) q[1];
rz(1.3364001) q[2];
sx q[2];
rz(-0.92051855) q[2];
sx q[2];
rz(-2.5087207) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74384159) q[1];
sx q[1];
rz(-1.4201179) q[1];
sx q[1];
rz(0.57445261) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7396163) q[3];
sx q[3];
rz(-2.4509015) q[3];
sx q[3];
rz(-2.9258779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56200999) q[2];
sx q[2];
rz(-1.1082114) q[2];
sx q[2];
rz(2.3902334) q[2];
rz(-0.56685081) q[3];
sx q[3];
rz(-1.6040498) q[3];
sx q[3];
rz(1.3273299) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71996561) q[0];
sx q[0];
rz(-0.1816853) q[0];
sx q[0];
rz(0.0075465329) q[0];
rz(-2.201572) q[1];
sx q[1];
rz(-0.66497856) q[1];
sx q[1];
rz(0.13253458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7203001) q[0];
sx q[0];
rz(-0.99908913) q[0];
sx q[0];
rz(-0.49457834) q[0];
rz(0.33905115) q[2];
sx q[2];
rz(-0.38290747) q[2];
sx q[2];
rz(-2.6660181) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30250994) q[1];
sx q[1];
rz(-1.1801475) q[1];
sx q[1];
rz(1.7724228) q[1];
rz(-pi) q[2];
rz(-1.2279195) q[3];
sx q[3];
rz(-1.7352805) q[3];
sx q[3];
rz(3.0066688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14564766) q[2];
sx q[2];
rz(-1.6463248) q[2];
sx q[2];
rz(-3.0403467) q[2];
rz(0.64949399) q[3];
sx q[3];
rz(-0.5972623) q[3];
sx q[3];
rz(-0.35972843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7262909) q[0];
sx q[0];
rz(-1.8446209) q[0];
sx q[0];
rz(-1.2671965) q[0];
rz(-1.2665117) q[1];
sx q[1];
rz(-0.15788618) q[1];
sx q[1];
rz(-2.3983541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38108724) q[0];
sx q[0];
rz(-0.60956565) q[0];
sx q[0];
rz(-2.7838092) q[0];
rz(-0.27979346) q[2];
sx q[2];
rz(-2.330707) q[2];
sx q[2];
rz(1.5992129) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2253902) q[1];
sx q[1];
rz(-1.4135201) q[1];
sx q[1];
rz(-1.3316403) q[1];
x q[2];
rz(-2.4489347) q[3];
sx q[3];
rz(-2.1941278) q[3];
sx q[3];
rz(2.3828363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1242975) q[2];
sx q[2];
rz(-2.9399019) q[2];
sx q[2];
rz(-1.4264433) q[2];
rz(-2.1258449) q[3];
sx q[3];
rz(-1.4371212) q[3];
sx q[3];
rz(2.3118741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5468686) q[0];
sx q[0];
rz(-2.4877553) q[0];
sx q[0];
rz(0.015901707) q[0];
rz(-1.6390027) q[1];
sx q[1];
rz(-0.67633164) q[1];
sx q[1];
rz(2.1471088) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976006) q[0];
sx q[0];
rz(-0.74403896) q[0];
sx q[0];
rz(0.97980325) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15764641) q[2];
sx q[2];
rz(-0.38262832) q[2];
sx q[2];
rz(-1.9668818) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2331824) q[1];
sx q[1];
rz(-0.841978) q[1];
sx q[1];
rz(0.89836095) q[1];
rz(-pi) q[2];
rz(1.8344457) q[3];
sx q[3];
rz(-0.52894512) q[3];
sx q[3];
rz(-2.6890713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9859621) q[2];
sx q[2];
rz(-1.4302284) q[2];
sx q[2];
rz(-2.963781) q[2];
rz(-1.6832247) q[3];
sx q[3];
rz(-0.42176133) q[3];
sx q[3];
rz(1.2678857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630702) q[0];
sx q[0];
rz(-1.2489742) q[0];
sx q[0];
rz(1.857969) q[0];
rz(-2.3578857) q[1];
sx q[1];
rz(-2.3678534) q[1];
sx q[1];
rz(-1.3073889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4625774) q[0];
sx q[0];
rz(-1.1333915) q[0];
sx q[0];
rz(1.6721729) q[0];
rz(1.6450591) q[2];
sx q[2];
rz(-0.38887244) q[2];
sx q[2];
rz(-2.996794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2210113) q[1];
sx q[1];
rz(-0.1836818) q[1];
sx q[1];
rz(1.8964975) q[1];
rz(0.41475716) q[3];
sx q[3];
rz(-2.0817753) q[3];
sx q[3];
rz(1.3172883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54138294) q[2];
sx q[2];
rz(-1.7317438) q[2];
sx q[2];
rz(0.48995885) q[2];
rz(1.1146891) q[3];
sx q[3];
rz(-1.1763108) q[3];
sx q[3];
rz(-2.8118242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212696) q[0];
sx q[0];
rz(-1.9135495) q[0];
sx q[0];
rz(1.3187153) q[0];
rz(-1.8439138) q[1];
sx q[1];
rz(-2.4083125) q[1];
sx q[1];
rz(-0.42793035) q[1];
rz(1.8280468) q[2];
sx q[2];
rz(-1.3197109) q[2];
sx q[2];
rz(-0.33725658) q[2];
rz(1.3379723) q[3];
sx q[3];
rz(-2.4216087) q[3];
sx q[3];
rz(-2.5523228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
