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
rz(2.7569438) q[0];
sx q[0];
rz(-0.83847133) q[0];
sx q[0];
rz(-0.60929259) q[0];
rz(0.72120178) q[1];
sx q[1];
rz(5.3546049) q[1];
sx q[1];
rz(7.3794853) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4247835) q[0];
sx q[0];
rz(-1.6832628) q[0];
sx q[0];
rz(0.87392585) q[0];
rz(-0.9651411) q[2];
sx q[2];
rz(-0.7363626) q[2];
sx q[2];
rz(3.1346653) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.83447245) q[1];
sx q[1];
rz(-2.5398919) q[1];
sx q[1];
rz(-0.51096075) q[1];
x q[2];
rz(-1.6814548) q[3];
sx q[3];
rz(-1.7250463) q[3];
sx q[3];
rz(-1.339104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91813749) q[2];
sx q[2];
rz(-1.8310941) q[2];
sx q[2];
rz(-1.9523331) q[2];
rz(-0.39092815) q[3];
sx q[3];
rz(-1.9783741) q[3];
sx q[3];
rz(2.0093567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44553462) q[0];
sx q[0];
rz(-1.3548509) q[0];
sx q[0];
rz(1.6433486) q[0];
rz(-2.4636726) q[1];
sx q[1];
rz(-1.8410212) q[1];
sx q[1];
rz(2.4300785) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2133117) q[0];
sx q[0];
rz(-0.60325256) q[0];
sx q[0];
rz(0.9071945) q[0];
x q[1];
rz(2.6783912) q[2];
sx q[2];
rz(-1.3579662) q[2];
sx q[2];
rz(0.45945278) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7603858) q[1];
sx q[1];
rz(-2.7863656) q[1];
sx q[1];
rz(-0.23204259) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-2.4352162) q[2];
sx q[2];
rz(-1.7855568) q[2];
sx q[2];
rz(2.0785275) q[2];
rz(-1.6478018) q[3];
sx q[3];
rz(-2.6726275) q[3];
sx q[3];
rz(-0.74023214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4110334) q[0];
sx q[0];
rz(-0.39091245) q[0];
sx q[0];
rz(-2.539769) q[0];
rz(0.14566323) q[1];
sx q[1];
rz(-2.115695) q[1];
sx q[1];
rz(-2.513733) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2780298) q[0];
sx q[0];
rz(-3.0042227) q[0];
sx q[0];
rz(-2.6111215) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7854087) q[2];
sx q[2];
rz(-1.0477475) q[2];
sx q[2];
rz(-1.5835536) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2035297) q[1];
sx q[1];
rz(-0.4258241) q[1];
sx q[1];
rz(2.1494206) q[1];
rz(-pi) q[2];
rz(-2.6791191) q[3];
sx q[3];
rz(-2.1331027) q[3];
sx q[3];
rz(3.1218046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7566028) q[2];
sx q[2];
rz(-1.0641229) q[2];
sx q[2];
rz(0.21558726) q[2];
rz(1.138341) q[3];
sx q[3];
rz(-1.3201069) q[3];
sx q[3];
rz(-0.57083541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0067921) q[0];
sx q[0];
rz(-1.726806) q[0];
sx q[0];
rz(0.11882812) q[0];
rz(1.4745332) q[1];
sx q[1];
rz(-0.88913616) q[1];
sx q[1];
rz(-2.3337591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12187639) q[0];
sx q[0];
rz(-2.3000679) q[0];
sx q[0];
rz(-2.3725933) q[0];
x q[1];
rz(-1.7052682) q[2];
sx q[2];
rz(-0.32609144) q[2];
sx q[2];
rz(1.6539314) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68271676) q[1];
sx q[1];
rz(-1.7380875) q[1];
sx q[1];
rz(-0.95750995) q[1];
x q[2];
rz(-1.586257) q[3];
sx q[3];
rz(-0.8925775) q[3];
sx q[3];
rz(2.79592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5011751) q[2];
sx q[2];
rz(-1.7265604) q[2];
sx q[2];
rz(-2.625722) q[2];
rz(-3.0473895) q[3];
sx q[3];
rz(-2.3661864) q[3];
sx q[3];
rz(1.5532106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696988) q[0];
sx q[0];
rz(-2.0409245) q[0];
sx q[0];
rz(1.9544253) q[0];
rz(0.59133235) q[1];
sx q[1];
rz(-1.2966803) q[1];
sx q[1];
rz(-0.82685131) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2555485) q[0];
sx q[0];
rz(-0.40103087) q[0];
sx q[0];
rz(0.31440763) q[0];
x q[1];
rz(-0.081080699) q[2];
sx q[2];
rz(-1.0141918) q[2];
sx q[2];
rz(-0.48559819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8137774) q[1];
sx q[1];
rz(-2.2499488) q[1];
sx q[1];
rz(-0.16736729) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99991344) q[3];
sx q[3];
rz(-2.0032845) q[3];
sx q[3];
rz(1.692358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4225215) q[2];
sx q[2];
rz(-2.8870388) q[2];
sx q[2];
rz(-1.0151601) q[2];
rz(-2.2319131) q[3];
sx q[3];
rz(-1.2404975) q[3];
sx q[3];
rz(-2.4801586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61843094) q[0];
sx q[0];
rz(-1.9310512) q[0];
sx q[0];
rz(0.30995187) q[0];
rz(3.1228206) q[1];
sx q[1];
rz(-1.4614146) q[1];
sx q[1];
rz(-3.054256) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7531484) q[0];
sx q[0];
rz(-2.5678621) q[0];
sx q[0];
rz(0.059772003) q[0];
rz(-2.3051065) q[2];
sx q[2];
rz(-2.2936471) q[2];
sx q[2];
rz(-1.7707056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.046459196) q[1];
sx q[1];
rz(-1.1756056) q[1];
sx q[1];
rz(2.1161152) q[1];
rz(-1.0424006) q[3];
sx q[3];
rz(-1.7070407) q[3];
sx q[3];
rz(-2.6353177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35149082) q[2];
sx q[2];
rz(-2.6678706) q[2];
sx q[2];
rz(-0.52004415) q[2];
rz(-0.13488787) q[3];
sx q[3];
rz(-1.8205732) q[3];
sx q[3];
rz(-1.7322056) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745558) q[0];
sx q[0];
rz(-2.0679857) q[0];
sx q[0];
rz(-2.7556162) q[0];
rz(-2.0416253) q[1];
sx q[1];
rz(-0.67953449) q[1];
sx q[1];
rz(-0.78752548) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3271017) q[0];
sx q[0];
rz(-1.9270032) q[0];
sx q[0];
rz(2.8339631) q[0];
rz(-pi) q[1];
rz(2.5239713) q[2];
sx q[2];
rz(-1.8810399) q[2];
sx q[2];
rz(-0.26981416) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3508151) q[1];
sx q[1];
rz(-2.3190089) q[1];
sx q[1];
rz(-1.0910396) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1925364) q[3];
sx q[3];
rz(-0.28655401) q[3];
sx q[3];
rz(-1.3499945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0218574) q[2];
sx q[2];
rz(-1.8334374) q[2];
sx q[2];
rz(-1.6278527) q[2];
rz(1.9205903) q[3];
sx q[3];
rz(-1.1648213) q[3];
sx q[3];
rz(0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30962238) q[0];
sx q[0];
rz(-1.7124875) q[0];
sx q[0];
rz(-2.9085462) q[0];
rz(-1.4876935) q[1];
sx q[1];
rz(-1.033604) q[1];
sx q[1];
rz(-2.1344562) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6876935) q[0];
sx q[0];
rz(-1.0570044) q[0];
sx q[0];
rz(-1.1790685) q[0];
rz(0.39733823) q[2];
sx q[2];
rz(-2.4724744) q[2];
sx q[2];
rz(0.64815176) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.050747063) q[1];
sx q[1];
rz(-1.7170893) q[1];
sx q[1];
rz(-1.1698224) q[1];
rz(1.5910025) q[3];
sx q[3];
rz(-1.8148225) q[3];
sx q[3];
rz(-3.1283825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2795589) q[2];
sx q[2];
rz(-2.0290012) q[2];
sx q[2];
rz(0.88279185) q[2];
rz(0.27859303) q[3];
sx q[3];
rz(-1.9735347) q[3];
sx q[3];
rz(-2.6826503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9423264) q[0];
sx q[0];
rz(-0.47573221) q[0];
sx q[0];
rz(1.5717773) q[0];
rz(2.3528631) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(-2.4615361) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2325033) q[0];
sx q[0];
rz(-1.6304558) q[0];
sx q[0];
rz(-0.76086126) q[0];
rz(2.3236507) q[2];
sx q[2];
rz(-2.9688058) q[2];
sx q[2];
rz(-2.8326952) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9763261) q[1];
sx q[1];
rz(-1.9496456) q[1];
sx q[1];
rz(-2.5514609) q[1];
rz(-0.84362883) q[3];
sx q[3];
rz(-1.9988201) q[3];
sx q[3];
rz(0.53800636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(-1.1663158) q[2];
rz(1.9370646) q[3];
sx q[3];
rz(-1.8185936) q[3];
sx q[3];
rz(-0.62674633) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.523943) q[0];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6740538) q[0];
sx q[0];
rz(-0.77097818) q[0];
sx q[0];
rz(1.5037554) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.885845) q[2];
sx q[2];
rz(-1.3761259) q[2];
sx q[2];
rz(0.48182975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5497351) q[1];
sx q[1];
rz(-0.80313528) q[1];
sx q[1];
rz(1.1863043) q[1];
rz(-pi) q[2];
rz(-0.043702262) q[3];
sx q[3];
rz(-0.80207295) q[3];
sx q[3];
rz(-1.613783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87307125) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(-0.29042563) q[2];
rz(0.63646603) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(1.9071473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.3270522) q[0];
sx q[0];
rz(-2.4335813) q[0];
sx q[0];
rz(-2.0306564) q[0];
rz(-3.0771599) q[1];
sx q[1];
rz(-1.4705407) q[1];
sx q[1];
rz(2.4992117) q[1];
rz(3.1101075) q[2];
sx q[2];
rz(-2.1654072) q[2];
sx q[2];
rz(-2.5624599) q[2];
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
