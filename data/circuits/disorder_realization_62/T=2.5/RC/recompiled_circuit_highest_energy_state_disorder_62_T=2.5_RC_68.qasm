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
rz(-0.38464889) q[0];
sx q[0];
rz(-2.3031213) q[0];
sx q[0];
rz(0.60929259) q[0];
rz(0.72120178) q[1];
sx q[1];
rz(-0.92858044) q[1];
sx q[1];
rz(-2.0452926) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168092) q[0];
sx q[0];
rz(-1.6832628) q[0];
sx q[0];
rz(2.2676668) q[0];
rz(-pi) q[1];
rz(-0.47639552) q[2];
sx q[2];
rz(-2.1557249) q[2];
sx q[2];
rz(-2.3830519) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9043353) q[1];
sx q[1];
rz(-1.0544027) q[1];
sx q[1];
rz(1.8947435) q[1];
x q[2];
rz(-2.9864086) q[3];
sx q[3];
rz(-1.6801356) q[3];
sx q[3];
rz(-2.9269699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91813749) q[2];
sx q[2];
rz(-1.8310941) q[2];
sx q[2];
rz(1.9523331) q[2];
rz(2.7506645) q[3];
sx q[3];
rz(-1.9783741) q[3];
sx q[3];
rz(2.0093567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696058) q[0];
sx q[0];
rz(-1.7867418) q[0];
sx q[0];
rz(-1.6433486) q[0];
rz(0.6779201) q[1];
sx q[1];
rz(-1.3005715) q[1];
sx q[1];
rz(0.71151412) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2133117) q[0];
sx q[0];
rz(-0.60325256) q[0];
sx q[0];
rz(0.9071945) q[0];
x q[1];
rz(0.46320148) q[2];
sx q[2];
rz(-1.3579662) q[2];
sx q[2];
rz(2.6821399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1700279) q[1];
sx q[1];
rz(-1.6508647) q[1];
sx q[1];
rz(2.7951334) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1101338) q[3];
sx q[3];
rz(-2.194187) q[3];
sx q[3];
rz(-1.7707847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4352162) q[2];
sx q[2];
rz(-1.7855568) q[2];
sx q[2];
rz(-2.0785275) q[2];
rz(1.4937909) q[3];
sx q[3];
rz(-2.6726275) q[3];
sx q[3];
rz(2.4013605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73055926) q[0];
sx q[0];
rz(-0.39091245) q[0];
sx q[0];
rz(2.539769) q[0];
rz(-0.14566323) q[1];
sx q[1];
rz(-2.115695) q[1];
sx q[1];
rz(-0.62785968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8191166) q[0];
sx q[0];
rz(-1.5014582) q[0];
sx q[0];
rz(3.022911) q[0];
rz(-pi) q[1];
rz(2.7854087) q[2];
sx q[2];
rz(-2.0938452) q[2];
sx q[2];
rz(-1.5580391) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9722189) q[1];
sx q[1];
rz(-1.3429321) q[1];
sx q[1];
rz(1.2078753) q[1];
rz(0.95730036) q[3];
sx q[3];
rz(-1.9578551) q[3];
sx q[3];
rz(-1.330803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7566028) q[2];
sx q[2];
rz(-1.0641229) q[2];
sx q[2];
rz(2.9260054) q[2];
rz(-2.0032517) q[3];
sx q[3];
rz(-1.3201069) q[3];
sx q[3];
rz(2.5707572) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1348006) q[0];
sx q[0];
rz(-1.4147867) q[0];
sx q[0];
rz(-3.0227645) q[0];
rz(-1.4745332) q[1];
sx q[1];
rz(-2.2524565) q[1];
sx q[1];
rz(0.80783358) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87619239) q[0];
sx q[0];
rz(-1.0256488) q[0];
sx q[0];
rz(2.4643174) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.045305552) q[2];
sx q[2];
rz(-1.8938365) q[2];
sx q[2];
rz(-1.3458061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3701757) q[1];
sx q[1];
rz(-0.96730212) q[1];
sx q[1];
rz(-2.9379528) q[1];
x q[2];
rz(-1.5553357) q[3];
sx q[3];
rz(-2.2490152) q[3];
sx q[3];
rz(-0.34567269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64041758) q[2];
sx q[2];
rz(-1.4150323) q[2];
sx q[2];
rz(-0.51587063) q[2];
rz(0.094203146) q[3];
sx q[3];
rz(-0.7754063) q[3];
sx q[3];
rz(-1.5532106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3696988) q[0];
sx q[0];
rz(-2.0409245) q[0];
sx q[0];
rz(-1.9544253) q[0];
rz(-0.59133235) q[1];
sx q[1];
rz(-1.2966803) q[1];
sx q[1];
rz(-2.3147413) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39385228) q[0];
sx q[0];
rz(-1.6918139) q[0];
sx q[0];
rz(-0.38328247) q[0];
rz(0.081080699) q[2];
sx q[2];
rz(-2.1274008) q[2];
sx q[2];
rz(-0.48559819) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0765269) q[1];
sx q[1];
rz(-2.4453116) q[1];
sx q[1];
rz(-1.3672921) q[1];
rz(-pi) q[2];
rz(-0.50179568) q[3];
sx q[3];
rz(-1.0579946) q[3];
sx q[3];
rz(2.7570908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7190711) q[2];
sx q[2];
rz(-0.25455385) q[2];
sx q[2];
rz(-2.1264326) q[2];
rz(-0.90967956) q[3];
sx q[3];
rz(-1.2404975) q[3];
sx q[3];
rz(-0.66143405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5231617) q[0];
sx q[0];
rz(-1.9310512) q[0];
sx q[0];
rz(0.30995187) q[0];
rz(0.018772086) q[1];
sx q[1];
rz(-1.680178) q[1];
sx q[1];
rz(-3.054256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3173125) q[0];
sx q[0];
rz(-2.1433733) q[0];
sx q[0];
rz(1.5322112) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83648615) q[2];
sx q[2];
rz(-0.84794551) q[2];
sx q[2];
rz(-1.370887) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0951335) q[1];
sx q[1];
rz(-1.9659871) q[1];
sx q[1];
rz(1.0254775) q[1];
rz(-pi) q[2];
rz(1.8363073) q[3];
sx q[3];
rz(-2.5975321) q[3];
sx q[3];
rz(1.8484704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7901018) q[2];
sx q[2];
rz(-2.6678706) q[2];
sx q[2];
rz(0.52004415) q[2];
rz(-0.13488787) q[3];
sx q[3];
rz(-1.8205732) q[3];
sx q[3];
rz(-1.7322056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.0745558) q[0];
sx q[0];
rz(-2.0679857) q[0];
sx q[0];
rz(-2.7556162) q[0];
rz(2.0416253) q[1];
sx q[1];
rz(-0.67953449) q[1];
sx q[1];
rz(-2.3540672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3271017) q[0];
sx q[0];
rz(-1.2145894) q[0];
sx q[0];
rz(-0.30762958) q[0];
rz(-pi) q[1];
rz(1.196127) q[2];
sx q[2];
rz(-0.9866937) q[2];
sx q[2];
rz(1.0874334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0213607) q[1];
sx q[1];
rz(-1.225705) q[1];
sx q[1];
rz(-0.80806513) q[1];
rz(-1.8058895) q[3];
sx q[3];
rz(-1.7361779) q[3];
sx q[3];
rz(-0.38136417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0218574) q[2];
sx q[2];
rz(-1.8334374) q[2];
sx q[2];
rz(1.6278527) q[2];
rz(1.9205903) q[3];
sx q[3];
rz(-1.9767714) q[3];
sx q[3];
rz(-0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30962238) q[0];
sx q[0];
rz(-1.4291052) q[0];
sx q[0];
rz(-0.23304644) q[0];
rz(-1.4876935) q[1];
sx q[1];
rz(-1.033604) q[1];
sx q[1];
rz(1.0071365) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171948) q[0];
sx q[0];
rz(-1.9097415) q[0];
sx q[0];
rz(-0.54817731) q[0];
x q[1];
rz(1.8677668) q[2];
sx q[2];
rz(-2.179702) q[2];
sx q[2];
rz(2.0021653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5598232) q[1];
sx q[1];
rz(-1.1743465) q[1];
sx q[1];
rz(0.15869402) q[1];
x q[2];
rz(1.5910025) q[3];
sx q[3];
rz(-1.3267702) q[3];
sx q[3];
rz(-0.013210162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86203376) q[2];
sx q[2];
rz(-1.1125914) q[2];
sx q[2];
rz(-2.2588008) q[2];
rz(2.8629996) q[3];
sx q[3];
rz(-1.9735347) q[3];
sx q[3];
rz(-0.45894233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19926628) q[0];
sx q[0];
rz(-0.47573221) q[0];
sx q[0];
rz(-1.5717773) q[0];
rz(0.78872952) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(2.4615361) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59916436) q[0];
sx q[0];
rz(-0.76272805) q[0];
sx q[0];
rz(0.086407089) q[0];
x q[1];
rz(3.0228258) q[2];
sx q[2];
rz(-1.4450018) q[2];
sx q[2];
rz(2.072352) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1627106) q[1];
sx q[1];
rz(-1.0274402) q[1];
sx q[1];
rz(-2.0175949) q[1];
rz(-pi) q[2];
rz(-2.1722776) q[3];
sx q[3];
rz(-2.318104) q[3];
sx q[3];
rz(1.4692509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(-1.1663158) q[2];
rz(-1.9370646) q[3];
sx q[3];
rz(-1.8185936) q[3];
sx q[3];
rz(0.62674633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6176497) q[0];
sx q[0];
rz(-0.88215041) q[0];
sx q[0];
rz(2.7045265) q[0];
rz(-0.79824671) q[1];
sx q[1];
rz(-0.50004807) q[1];
sx q[1];
rz(0.22013586) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4675389) q[0];
sx q[0];
rz(-2.3706145) q[0];
sx q[0];
rz(-1.5037554) q[0];
x q[1];
rz(-1.0041084) q[2];
sx q[2];
rz(-2.7729642) q[2];
sx q[2];
rz(-1.6246375) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1194296) q[1];
sx q[1];
rz(-2.3009752) q[1];
sx q[1];
rz(2.7709211) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.043702262) q[3];
sx q[3];
rz(-0.80207295) q[3];
sx q[3];
rz(1.5278096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.87307125) q[2];
sx q[2];
rz(-1.9067418) q[2];
sx q[2];
rz(-0.29042563) q[2];
rz(-2.5051266) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(1.9071473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3270522) q[0];
sx q[0];
rz(-0.70801133) q[0];
sx q[0];
rz(1.1109362) q[0];
rz(-0.064432714) q[1];
sx q[1];
rz(-1.671052) q[1];
sx q[1];
rz(-0.64238092) q[1];
rz(1.6173132) q[2];
sx q[2];
rz(-2.5462493) q[2];
sx q[2];
rz(0.52296849) q[2];
rz(-1.3051582) q[3];
sx q[3];
rz(-2.0943741) q[3];
sx q[3];
rz(1.5423273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
