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
rz(-0.77684075) q[0];
sx q[0];
rz(2.2651894) q[0];
sx q[0];
rz(9.6776008) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(4.0568772) q[1];
sx q[1];
rz(8.6203909) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641805) q[0];
sx q[0];
rz(-1.2578811) q[0];
sx q[0];
rz(1.1699647) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8041274) q[2];
sx q[2];
rz(-1.9905645) q[2];
sx q[2];
rz(2.6884318) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.731718) q[1];
sx q[1];
rz(-1.2267641) q[1];
sx q[1];
rz(0.39459497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8283903) q[3];
sx q[3];
rz(-2.7814908) q[3];
sx q[3];
rz(-1.2822156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.09314166) q[2];
sx q[2];
rz(-0.96258771) q[2];
sx q[2];
rz(-2.2691881) q[2];
rz(0.33556542) q[3];
sx q[3];
rz(-1.7458956) q[3];
sx q[3];
rz(0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63808477) q[0];
sx q[0];
rz(-0.26569772) q[0];
sx q[0];
rz(-0.22915325) q[0];
rz(2.9511662) q[1];
sx q[1];
rz(-1.8016022) q[1];
sx q[1];
rz(0.71358877) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23799831) q[0];
sx q[0];
rz(-2.2269214) q[0];
sx q[0];
rz(0.34060069) q[0];
x q[1];
rz(0.83723703) q[2];
sx q[2];
rz(-1.9847437) q[2];
sx q[2];
rz(2.4606103) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5428764) q[1];
sx q[1];
rz(-2.2440254) q[1];
sx q[1];
rz(1.9510516) q[1];
rz(-pi) q[2];
rz(2.3768164) q[3];
sx q[3];
rz(-1.6891317) q[3];
sx q[3];
rz(1.8272994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68636346) q[2];
sx q[2];
rz(-2.8296622) q[2];
sx q[2];
rz(-0.4757821) q[2];
rz(2.0717715) q[3];
sx q[3];
rz(-2.1333623) q[3];
sx q[3];
rz(2.3750335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0692724) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(1.9092165) q[0];
rz(2.6909289) q[1];
sx q[1];
rz(-1.181299) q[1];
sx q[1];
rz(-0.88952363) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7018676) q[0];
sx q[0];
rz(-0.98977463) q[0];
sx q[0];
rz(2.0608611) q[0];
rz(-pi) q[1];
rz(2.6085147) q[2];
sx q[2];
rz(-2.0642082) q[2];
sx q[2];
rz(0.62668884) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8107618) q[1];
sx q[1];
rz(-2.2137801) q[1];
sx q[1];
rz(0.26442702) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37556383) q[3];
sx q[3];
rz(-1.2455539) q[3];
sx q[3];
rz(-0.90876814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4543317) q[2];
sx q[2];
rz(-2.7165518) q[2];
sx q[2];
rz(-1.2588151) q[2];
rz(1.2339633) q[3];
sx q[3];
rz(-2.0089269) q[3];
sx q[3];
rz(-3.1335355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312129) q[0];
sx q[0];
rz(-0.90573913) q[0];
sx q[0];
rz(1.4696962) q[0];
rz(-1.9944893) q[1];
sx q[1];
rz(-1.0018145) q[1];
sx q[1];
rz(-0.9009487) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1484447) q[0];
sx q[0];
rz(-1.6560418) q[0];
sx q[0];
rz(2.8277525) q[0];
x q[1];
rz(0.78184266) q[2];
sx q[2];
rz(-0.78357176) q[2];
sx q[2];
rz(1.9279566) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21877737) q[1];
sx q[1];
rz(-2.4115857) q[1];
sx q[1];
rz(-0.14132146) q[1];
x q[2];
rz(-0.17662405) q[3];
sx q[3];
rz(-0.7407397) q[3];
sx q[3];
rz(2.3185042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49742302) q[2];
sx q[2];
rz(-1.3966509) q[2];
sx q[2];
rz(2.4141342) q[2];
rz(2.4265031) q[3];
sx q[3];
rz(-0.46052027) q[3];
sx q[3];
rz(0.49811825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46399507) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(0.736262) q[0];
rz(-2.5212506) q[1];
sx q[1];
rz(-2.5963929) q[1];
sx q[1];
rz(1.6166519) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0382181) q[0];
sx q[0];
rz(-0.056110121) q[0];
sx q[0];
rz(-2.9655064) q[0];
rz(0.97648804) q[2];
sx q[2];
rz(-0.39820489) q[2];
sx q[2];
rz(-1.9081685) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1641008) q[1];
sx q[1];
rz(-2.3445446) q[1];
sx q[1];
rz(2.1063095) q[1];
rz(-pi) q[2];
rz(-0.81186632) q[3];
sx q[3];
rz(-0.97211876) q[3];
sx q[3];
rz(-0.52385274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5357431) q[2];
sx q[2];
rz(-0.77631408) q[2];
sx q[2];
rz(-2.1898451) q[2];
rz(0.4176628) q[3];
sx q[3];
rz(-2.8916736) q[3];
sx q[3];
rz(1.1550268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34053892) q[0];
sx q[0];
rz(-1.8908353) q[0];
sx q[0];
rz(0.75403768) q[0];
rz(0.89318371) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(2.975614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8360031) q[0];
sx q[0];
rz(-2.2512332) q[0];
sx q[0];
rz(-1.7050118) q[0];
x q[1];
rz(2.6352623) q[2];
sx q[2];
rz(-1.4131318) q[2];
sx q[2];
rz(-0.28749301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7268035) q[1];
sx q[1];
rz(-0.94894743) q[1];
sx q[1];
rz(-1.0110823) q[1];
x q[2];
rz(0.13940553) q[3];
sx q[3];
rz(-1.8434815) q[3];
sx q[3];
rz(0.26095873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37504998) q[2];
sx q[2];
rz(-1.1376209) q[2];
sx q[2];
rz(-3.1362015) q[2];
rz(-3.052875) q[3];
sx q[3];
rz(-1.6088477) q[3];
sx q[3];
rz(2.3458792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8832815) q[0];
sx q[0];
rz(-1.2092051) q[0];
sx q[0];
rz(-0.2051556) q[0];
rz(-0.908665) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(3.0970792) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1703029) q[0];
sx q[0];
rz(-1.4409541) q[0];
sx q[0];
rz(-1.7287888) q[0];
rz(0.037864121) q[2];
sx q[2];
rz(-1.711039) q[2];
sx q[2];
rz(0.87272206) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8537112) q[1];
sx q[1];
rz(-1.3559623) q[1];
sx q[1];
rz(-0.98835215) q[1];
rz(-pi) q[2];
rz(-1.6904852) q[3];
sx q[3];
rz(-1.5482248) q[3];
sx q[3];
rz(-2.165497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32555106) q[2];
sx q[2];
rz(-2.6191235) q[2];
sx q[2];
rz(-2.6204056) q[2];
rz(-2.9442287) q[3];
sx q[3];
rz(-1.7557996) q[3];
sx q[3];
rz(2.638773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407783) q[0];
sx q[0];
rz(-0.1939119) q[0];
sx q[0];
rz(2.8863696) q[0];
rz(-2.9995645) q[1];
sx q[1];
rz(-1.7250926) q[1];
sx q[1];
rz(-2.7878888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3981374) q[0];
sx q[0];
rz(-2.6530735) q[0];
sx q[0];
rz(0.81870458) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92112215) q[2];
sx q[2];
rz(-2.4068458) q[2];
sx q[2];
rz(0.91780797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6981909) q[1];
sx q[1];
rz(-0.56487067) q[1];
sx q[1];
rz(2.2400212) q[1];
x q[2];
rz(-1.0068302) q[3];
sx q[3];
rz(-2.2368397) q[3];
sx q[3];
rz(-0.49426916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8346617) q[2];
sx q[2];
rz(-0.31535172) q[2];
sx q[2];
rz(-1.6174512) q[2];
rz(-0.8664242) q[3];
sx q[3];
rz(-1.6774991) q[3];
sx q[3];
rz(0.58969897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50188142) q[0];
sx q[0];
rz(-0.15637936) q[0];
sx q[0];
rz(1.7891275) q[0];
rz(-1.2347219) q[1];
sx q[1];
rz(-2.9094978) q[1];
sx q[1];
rz(0.51188767) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6366301) q[0];
sx q[0];
rz(-1.413336) q[0];
sx q[0];
rz(0.053863346) q[0];
x q[1];
rz(0.76546957) q[2];
sx q[2];
rz(-1.6846894) q[2];
sx q[2];
rz(1.423045) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8066387) q[1];
sx q[1];
rz(-2.6158164) q[1];
sx q[1];
rz(0.61383617) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7861106) q[3];
sx q[3];
rz(-2.1256697) q[3];
sx q[3];
rz(1.2303838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0074244) q[2];
sx q[2];
rz(-1.8210501) q[2];
sx q[2];
rz(-0.39719886) q[2];
rz(2.126179) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(-2.5889682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.7427202) q[0];
sx q[0];
rz(-1.9341368) q[0];
sx q[0];
rz(-2.7040828) q[0];
rz(-2.4028026) q[1];
sx q[1];
rz(-0.80016017) q[1];
sx q[1];
rz(-1.4791666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5115868) q[0];
sx q[0];
rz(-1.5151205) q[0];
sx q[0];
rz(-0.48752174) q[0];
rz(1.871256) q[2];
sx q[2];
rz(-2.766937) q[2];
sx q[2];
rz(-0.38196358) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65867546) q[1];
sx q[1];
rz(-1.632893) q[1];
sx q[1];
rz(3.0351522) q[1];
rz(0.78082943) q[3];
sx q[3];
rz(-1.025628) q[3];
sx q[3];
rz(3.0913251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9208357) q[2];
sx q[2];
rz(-1.3584542) q[2];
sx q[2];
rz(-2.9662568) q[2];
rz(-1.0026503) q[3];
sx q[3];
rz(-2.9084539) q[3];
sx q[3];
rz(1.9076777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6709082) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(-2.9910174) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(-0.19828037) q[2];
sx q[2];
rz(-1.5783327) q[2];
sx q[2];
rz(-0.21680149) q[2];
rz(-0.2507052) q[3];
sx q[3];
rz(-1.7964994) q[3];
sx q[3];
rz(-2.9356706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
