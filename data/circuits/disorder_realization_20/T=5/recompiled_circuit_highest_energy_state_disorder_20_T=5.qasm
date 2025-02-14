OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.196584954857826) q[0];
sx q[0];
rz(6.72533336480195) q[0];
sx q[0];
rz(10.2861767172734) q[0];
rz(1.50988936424255) q[1];
sx q[1];
rz(1.32465258439118) q[1];
sx q[1];
rz(9.46037046834036) q[1];
cx q[1],q[0];
rz(-1.13827681541443) q[0];
sx q[0];
rz(3.5328929742151) q[0];
sx q[0];
rz(10.4615308999936) q[0];
rz(-1.90274727344513) q[2];
sx q[2];
rz(6.60620537598664) q[2];
sx q[2];
rz(8.49776039122745) q[2];
cx q[2],q[1];
rz(-1.50006771087646) q[1];
sx q[1];
rz(4.74817314942414) q[1];
sx q[1];
rz(10.2689533591191) q[1];
rz(-0.0222363285720348) q[3];
sx q[3];
rz(2.46148684819276) q[3];
sx q[3];
rz(10.2487578749578) q[3];
cx q[3],q[2];
rz(3.69437170028687) q[2];
sx q[2];
rz(1.84010437329347) q[2];
sx q[2];
rz(4.23608014582797) q[2];
rz(1.62574458122253) q[3];
sx q[3];
rz(4.01418391068513) q[3];
sx q[3];
rz(12.4227909803311) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.04824078083038) q[0];
sx q[0];
rz(5.29291883309419) q[0];
sx q[0];
rz(9.05070558785602) q[0];
rz(-2.28349423408508) q[1];
sx q[1];
rz(7.22399726708467) q[1];
sx q[1];
rz(13.6421680211942) q[1];
cx q[1],q[0];
rz(2.61585235595703) q[0];
sx q[0];
rz(4.81554821332032) q[0];
sx q[0];
rz(7.75149474143192) q[0];
rz(3.38069009780884) q[2];
sx q[2];
rz(7.31930557091767) q[2];
sx q[2];
rz(10.2561032533567) q[2];
cx q[2],q[1];
rz(1.08292162418365) q[1];
sx q[1];
rz(1.45450142224366) q[1];
sx q[1];
rz(10.5937975406568) q[1];
rz(2.33874368667603) q[3];
sx q[3];
rz(3.98134532769258) q[3];
sx q[3];
rz(8.92142353057071) q[3];
cx q[3],q[2];
rz(-5.76797294616699) q[2];
sx q[2];
rz(0.433128031092235) q[2];
sx q[2];
rz(5.19997069834873) q[2];
rz(1.29272389411926) q[3];
sx q[3];
rz(5.14953246911103) q[3];
sx q[3];
rz(13.4950361013333) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.25562965869904) q[0];
sx q[0];
rz(3.91274723609025) q[0];
sx q[0];
rz(8.89987776278659) q[0];
rz(1.45957541465759) q[1];
sx q[1];
rz(5.54558316071565) q[1];
sx q[1];
rz(6.13739392756625) q[1];
cx q[1],q[0];
rz(-0.735066711902618) q[0];
sx q[0];
rz(3.00220884581143) q[0];
sx q[0];
rz(10.2835435628812) q[0];
rz(5.24611568450928) q[2];
sx q[2];
rz(4.17820552189881) q[2];
sx q[2];
rz(10.923718905441) q[2];
cx q[2],q[1];
rz(3.8940486907959) q[1];
sx q[1];
rz(3.58430731494958) q[1];
sx q[1];
rz(8.06573030947849) q[1];
rz(2.41130614280701) q[3];
sx q[3];
rz(5.15999189217622) q[3];
sx q[3];
rz(5.42979619502231) q[3];
cx q[3],q[2];
rz(-0.244744569063187) q[2];
sx q[2];
rz(1.16200986702973) q[2];
sx q[2];
rz(6.44858620165988) q[2];
rz(-0.861697196960449) q[3];
sx q[3];
rz(2.44980368216569) q[3];
sx q[3];
rz(9.22797333299323) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.86213135719299) q[0];
sx q[0];
rz(3.66583427985246) q[0];
sx q[0];
rz(7.01142332553073) q[0];
rz(-2.81844520568848) q[1];
sx q[1];
rz(4.1757871230417) q[1];
sx q[1];
rz(11.5271556138913) q[1];
cx q[1],q[0];
rz(-0.820742189884186) q[0];
sx q[0];
rz(3.74109903176362) q[0];
sx q[0];
rz(8.67073187827274) q[0];
rz(0.372972726821899) q[2];
sx q[2];
rz(2.47142377694184) q[2];
sx q[2];
rz(12.6916172265927) q[2];
cx q[2],q[1];
rz(-7.33281564712524) q[1];
sx q[1];
rz(6.77660981019075) q[1];
sx q[1];
rz(11.191714501373) q[1];
rz(-0.731527805328369) q[3];
sx q[3];
rz(5.38454881508882) q[3];
sx q[3];
rz(10.6356780290525) q[3];
cx q[3],q[2];
rz(-0.491796582937241) q[2];
sx q[2];
rz(4.19205907185609) q[2];
sx q[2];
rz(12.0473859071653) q[2];
rz(1.07458782196045) q[3];
sx q[3];
rz(4.99736765225465) q[3];
sx q[3];
rz(9.92004392146274) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.07866668701172) q[0];
sx q[0];
rz(4.06359562476213) q[0];
sx q[0];
rz(9.25320412813827) q[0];
rz(2.04731392860413) q[1];
sx q[1];
rz(1.84058740933473) q[1];
sx q[1];
rz(9.65938524007007) q[1];
cx q[1],q[0];
rz(0.0714834704995155) q[0];
sx q[0];
rz(3.44869592984254) q[0];
sx q[0];
rz(13.0016259908597) q[0];
rz(-2.07803297042847) q[2];
sx q[2];
rz(4.66647842724855) q[2];
sx q[2];
rz(9.40754423513218) q[2];
cx q[2],q[1];
rz(-1.80869102478027) q[1];
sx q[1];
rz(4.11764100392396) q[1];
sx q[1];
rz(7.15160987376376) q[1];
rz(1.14193141460419) q[3];
sx q[3];
rz(5.52585426171357) q[3];
sx q[3];
rz(8.28996715544864) q[3];
cx q[3],q[2];
rz(2.22850632667542) q[2];
sx q[2];
rz(1.19149103959138) q[2];
sx q[2];
rz(8.02202866076633) q[2];
rz(-3.76630306243896) q[3];
sx q[3];
rz(4.11000082095201) q[3];
sx q[3];
rz(10.7679621934812) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.60032820701599) q[0];
sx q[0];
rz(-0.00774892966216001) q[0];
sx q[0];
rz(9.09724340438052) q[0];
rz(-0.0281183570623398) q[1];
sx q[1];
rz(5.13940468628938) q[1];
sx q[1];
rz(13.5580959081571) q[1];
cx q[1],q[0];
rz(0.926279962062836) q[0];
sx q[0];
rz(4.01179656584794) q[0];
sx q[0];
rz(10.3051185965459) q[0];
rz(-3.25027561187744) q[2];
sx q[2];
rz(1.79554882844026) q[2];
sx q[2];
rz(10.5207521676938) q[2];
cx q[2],q[1];
rz(3.50213408470154) q[1];
sx q[1];
rz(7.09868493874604) q[1];
sx q[1];
rz(14.1947397947232) q[1];
rz(-0.0273689515888691) q[3];
sx q[3];
rz(5.18296924431855) q[3];
sx q[3];
rz(8.25074801444217) q[3];
cx q[3],q[2];
rz(1.23453867435455) q[2];
sx q[2];
rz(4.90552106698091) q[2];
sx q[2];
rz(7.76554427146121) q[2];
rz(-2.07568573951721) q[3];
sx q[3];
rz(1.96124711831147) q[3];
sx q[3];
rz(10.5218020439069) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.75575733184814) q[0];
sx q[0];
rz(3.25339993287856) q[0];
sx q[0];
rz(12.2709390878598) q[0];
rz(0.237519860267639) q[1];
sx q[1];
rz(4.07529720862443) q[1];
sx q[1];
rz(7.03056261538669) q[1];
cx q[1],q[0];
rz(-0.454936534166336) q[0];
sx q[0];
rz(2.12244084675843) q[0];
sx q[0];
rz(10.0281706809919) q[0];
rz(0.53520268201828) q[2];
sx q[2];
rz(5.26965514023835) q[2];
sx q[2];
rz(13.6984205007474) q[2];
cx q[2],q[1];
rz(0.574953675270081) q[1];
sx q[1];
rz(3.77800193627412) q[1];
sx q[1];
rz(10.0782751798551) q[1];
rz(2.76678323745728) q[3];
sx q[3];
rz(5.22458592255647) q[3];
sx q[3];
rz(8.4678297996442) q[3];
cx q[3],q[2];
rz(-1.47330093383789) q[2];
sx q[2];
rz(4.3273653109842) q[2];
sx q[2];
rz(12.8212709188382) q[2];
rz(-0.212368026375771) q[3];
sx q[3];
rz(5.41871372063691) q[3];
sx q[3];
rz(9.42574646715449) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.851069867610931) q[0];
sx q[0];
rz(5.05932548840577) q[0];
sx q[0];
rz(9.54867823272153) q[0];
rz(0.875214159488678) q[1];
sx q[1];
rz(5.45018616517121) q[1];
sx q[1];
rz(8.86599067448779) q[1];
cx q[1],q[0];
rz(-0.669847190380096) q[0];
sx q[0];
rz(2.55152431328828) q[0];
sx q[0];
rz(9.88558111190006) q[0];
rz(2.64260411262512) q[2];
sx q[2];
rz(1.82154157956178) q[2];
sx q[2];
rz(12.285500740997) q[2];
cx q[2],q[1];
rz(-2.32718229293823) q[1];
sx q[1];
rz(6.88722530205781) q[1];
sx q[1];
rz(8.34208462237521) q[1];
rz(-0.540395677089691) q[3];
sx q[3];
rz(2.43217876751954) q[3];
sx q[3];
rz(10.8588704824369) q[3];
cx q[3],q[2];
rz(0.0458040684461594) q[2];
sx q[2];
rz(5.79147473176057) q[2];
sx q[2];
rz(11.845628476135) q[2];
rz(-0.363024294376373) q[3];
sx q[3];
rz(4.36530795891816) q[3];
sx q[3];
rz(10.9365296125333) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.956612884998322) q[0];
sx q[0];
rz(2.94738011260564) q[0];
sx q[0];
rz(9.51383421420261) q[0];
rz(-0.366742759943008) q[1];
sx q[1];
rz(4.04583504994447) q[1];
sx q[1];
rz(10.8843482494275) q[1];
cx q[1],q[0];
rz(0.372925907373428) q[0];
sx q[0];
rz(3.81409773428971) q[0];
sx q[0];
rz(9.79623491167232) q[0];
rz(-4.0310959815979) q[2];
sx q[2];
rz(5.47417131264741) q[2];
sx q[2];
rz(11.3743057012479) q[2];
cx q[2],q[1];
rz(-0.68946099281311) q[1];
sx q[1];
rz(5.39917531807954) q[1];
sx q[1];
rz(10.5124530553739) q[1];
rz(0.437742412090302) q[3];
sx q[3];
rz(5.26600113709504) q[3];
sx q[3];
rz(10.0065182208936) q[3];
cx q[3],q[2];
rz(4.3480920791626) q[2];
sx q[2];
rz(4.50879541237886) q[2];
sx q[2];
rz(7.53326389788791) q[2];
rz(1.91532039642334) q[3];
sx q[3];
rz(5.64926281769807) q[3];
sx q[3];
rz(13.2053556203763) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.5117689371109) q[0];
sx q[0];
rz(-1.1331695000357) q[0];
sx q[0];
rz(10.5647998809735) q[0];
rz(-5.70006942749023) q[1];
sx q[1];
rz(4.62802544434602) q[1];
sx q[1];
rz(11.8985411882321) q[1];
cx q[1],q[0];
rz(-0.443964689970016) q[0];
sx q[0];
rz(5.42472687562043) q[0];
sx q[0];
rz(8.98878896831676) q[0];
rz(6.00295972824097) q[2];
sx q[2];
rz(3.86310848792131) q[2];
sx q[2];
rz(7.30361816882297) q[2];
cx q[2],q[1];
rz(1.78575038909912) q[1];
sx q[1];
rz(4.51792720158631) q[1];
sx q[1];
rz(9.8654128074567) q[1];
rz(2.02328085899353) q[3];
sx q[3];
rz(7.12877479394013) q[3];
sx q[3];
rz(12.4552106618802) q[3];
cx q[3],q[2];
rz(-3.98552417755127) q[2];
sx q[2];
rz(5.78829041321809) q[2];
sx q[2];
rz(5.53097174166843) q[2];
rz(3.06098675727844) q[3];
sx q[3];
rz(4.9335357268625) q[3];
sx q[3];
rz(8.8772524356763) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.232038393616676) q[0];
sx q[0];
rz(4.76974061329896) q[0];
sx q[0];
rz(6.23910400866672) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.695601046085358) q[1];
sx q[1];
rz(5.43519416649873) q[1];
sx q[1];
rz(7.75234231948062) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.948902249336243) q[2];
sx q[2];
rz(0.549336107569285) q[2];
sx q[2];
rz(10.961330986015) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.06864523887634) q[3];
sx q[3];
rz(5.00718358357484) q[3];
sx q[3];
rz(12.1105193853299) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
