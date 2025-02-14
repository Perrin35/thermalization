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
rz(-2.0105536) q[0];
sx q[0];
rz(-2.5656963) q[0];
sx q[0];
rz(2.0146712) q[0];
rz(2.1908886) q[1];
sx q[1];
rz(-2.5403412) q[1];
sx q[1];
rz(0.33369219) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2017444) q[0];
sx q[0];
rz(-0.44838312) q[0];
sx q[0];
rz(1.2192307) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26716455) q[2];
sx q[2];
rz(-2.0120624) q[2];
sx q[2];
rz(1.522832) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0970983) q[1];
sx q[1];
rz(-1.3584029) q[1];
sx q[1];
rz(-2.7907759) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80218704) q[3];
sx q[3];
rz(-1.3787375) q[3];
sx q[3];
rz(3.1046228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2600962) q[2];
sx q[2];
rz(-2.0659955) q[2];
sx q[2];
rz(-2.0913701) q[2];
rz(2.8460734) q[3];
sx q[3];
rz(-0.76885709) q[3];
sx q[3];
rz(1.3692921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1693901) q[0];
sx q[0];
rz(-1.0193595) q[0];
sx q[0];
rz(0.66725677) q[0];
rz(2.9908906) q[1];
sx q[1];
rz(-1.0548016) q[1];
sx q[1];
rz(-2.8299455) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3903605) q[0];
sx q[0];
rz(-1.6412819) q[0];
sx q[0];
rz(-1.0485534) q[0];
x q[1];
rz(2.8337651) q[2];
sx q[2];
rz(-1.9914978) q[2];
sx q[2];
rz(1.4811279) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2296621) q[1];
sx q[1];
rz(-1.7102825) q[1];
sx q[1];
rz(0.18499495) q[1];
rz(-pi) q[2];
rz(-2.9695791) q[3];
sx q[3];
rz(-0.40296754) q[3];
sx q[3];
rz(2.1864219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47505891) q[2];
sx q[2];
rz(-0.89343137) q[2];
sx q[2];
rz(0.7676355) q[2];
rz(2.660699) q[3];
sx q[3];
rz(-0.77156639) q[3];
sx q[3];
rz(1.5546999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.84592205) q[0];
sx q[0];
rz(-1.6039811) q[0];
sx q[0];
rz(-2.3620102) q[0];
rz(-2.7698611) q[1];
sx q[1];
rz(-1.0075684) q[1];
sx q[1];
rz(-0.76098162) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5190353) q[0];
sx q[0];
rz(-1.8599959) q[0];
sx q[0];
rz(-1.2361121) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43698796) q[2];
sx q[2];
rz(-2.2088364) q[2];
sx q[2];
rz(0.47666046) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8345388) q[1];
sx q[1];
rz(-2.4014856) q[1];
sx q[1];
rz(1.3241052) q[1];
x q[2];
rz(0.85270057) q[3];
sx q[3];
rz(-1.5745244) q[3];
sx q[3];
rz(-2.0223597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3373105) q[2];
sx q[2];
rz(-2.2411942) q[2];
sx q[2];
rz(0.4551746) q[2];
rz(0.69194397) q[3];
sx q[3];
rz(-0.71440905) q[3];
sx q[3];
rz(-0.80879319) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96875018) q[0];
sx q[0];
rz(-1.3469561) q[0];
sx q[0];
rz(0.2562879) q[0];
rz(-1.7068656) q[1];
sx q[1];
rz(-1.7211569) q[1];
sx q[1];
rz(3.0228379) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1349604) q[0];
sx q[0];
rz(-1.4156355) q[0];
sx q[0];
rz(2.6895903) q[0];
rz(-1.2776472) q[2];
sx q[2];
rz(-1.6934694) q[2];
sx q[2];
rz(2.6495767) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6765457) q[1];
sx q[1];
rz(-0.80261723) q[1];
sx q[1];
rz(1.8957183) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9471274) q[3];
sx q[3];
rz(-2.2897165) q[3];
sx q[3];
rz(0.014499078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28748301) q[2];
sx q[2];
rz(-2.867925) q[2];
sx q[2];
rz(2.0859094) q[2];
rz(1.4500729) q[3];
sx q[3];
rz(-1.6943211) q[3];
sx q[3];
rz(2.292574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5797822) q[0];
sx q[0];
rz(-2.9636443) q[0];
sx q[0];
rz(1.5280888) q[0];
rz(-1.1771419) q[1];
sx q[1];
rz(-1.8714995) q[1];
sx q[1];
rz(-1.5783763) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10229853) q[0];
sx q[0];
rz(-1.5015409) q[0];
sx q[0];
rz(2.5152492) q[0];
rz(-pi) q[1];
rz(2.7308942) q[2];
sx q[2];
rz(-1.8735527) q[2];
sx q[2];
rz(1.2891642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39644966) q[1];
sx q[1];
rz(-0.94303326) q[1];
sx q[1];
rz(-1.3954074) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5734736) q[3];
sx q[3];
rz(-1.9802046) q[3];
sx q[3];
rz(-2.1684627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13581181) q[2];
sx q[2];
rz(-1.9848738) q[2];
sx q[2];
rz(-1.8394252) q[2];
rz(2.8016413) q[3];
sx q[3];
rz(-0.8067185) q[3];
sx q[3];
rz(-0.66238856) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.133701) q[0];
sx q[0];
rz(-1.5971203) q[0];
sx q[0];
rz(-1.1997461) q[0];
rz(1.7802736) q[1];
sx q[1];
rz(-0.97313762) q[1];
sx q[1];
rz(2.5637085) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91980714) q[0];
sx q[0];
rz(-0.01520201) q[0];
sx q[0];
rz(1.1342681) q[0];
x q[1];
rz(1.4125713) q[2];
sx q[2];
rz(-1.998868) q[2];
sx q[2];
rz(-0.23382631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99022229) q[1];
sx q[1];
rz(-1.434902) q[1];
sx q[1];
rz(2.6841867) q[1];
rz(-3.1158524) q[3];
sx q[3];
rz(-1.8696068) q[3];
sx q[3];
rz(2.4853277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.63152385) q[2];
sx q[2];
rz(-2.5670467) q[2];
sx q[2];
rz(-2.6681382) q[2];
rz(1.8691285) q[3];
sx q[3];
rz(-1.748964) q[3];
sx q[3];
rz(-0.22245358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94721395) q[0];
sx q[0];
rz(-0.65776238) q[0];
sx q[0];
rz(2.8084602) q[0];
rz(0.22854742) q[1];
sx q[1];
rz(-1.8014149) q[1];
sx q[1];
rz(-2.5659335) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7299907) q[0];
sx q[0];
rz(-1.6766929) q[0];
sx q[0];
rz(0.81231711) q[0];
x q[1];
rz(-0.47708738) q[2];
sx q[2];
rz(-1.8680352) q[2];
sx q[2];
rz(1.5971668) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54652366) q[1];
sx q[1];
rz(-1.2449588) q[1];
sx q[1];
rz(2.930546) q[1];
rz(-pi) q[2];
rz(-2.6603465) q[3];
sx q[3];
rz(-1.7982881) q[3];
sx q[3];
rz(2.9862822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8754742) q[2];
sx q[2];
rz(-0.55142752) q[2];
sx q[2];
rz(1.6654588) q[2];
rz(-2.0406145) q[3];
sx q[3];
rz(-2.0783547) q[3];
sx q[3];
rz(-0.5180009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5975033) q[0];
sx q[0];
rz(-0.88633716) q[0];
sx q[0];
rz(-0.30712095) q[0];
rz(0.33044526) q[1];
sx q[1];
rz(-1.5437061) q[1];
sx q[1];
rz(-1.7764567) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6520321) q[0];
sx q[0];
rz(-0.083261641) q[0];
sx q[0];
rz(-0.21647446) q[0];
x q[1];
rz(-2.6583985) q[2];
sx q[2];
rz(-0.38887923) q[2];
sx q[2];
rz(-0.70975297) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0338034) q[1];
sx q[1];
rz(-0.50869903) q[1];
sx q[1];
rz(2.1230039) q[1];
rz(1.4151814) q[3];
sx q[3];
rz(-2.3978516) q[3];
sx q[3];
rz(1.9187601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27552989) q[2];
sx q[2];
rz(-0.58774647) q[2];
sx q[2];
rz(2.8928939) q[2];
rz(-1.9281049) q[3];
sx q[3];
rz(-1.4444193) q[3];
sx q[3];
rz(0.061633751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1520749) q[0];
sx q[0];
rz(-0.9062506) q[0];
sx q[0];
rz(-0.686598) q[0];
rz(-1.1129414) q[1];
sx q[1];
rz(-0.80250347) q[1];
sx q[1];
rz(-0.31563219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0604117) q[0];
sx q[0];
rz(-0.32627772) q[0];
sx q[0];
rz(-1.7148036) q[0];
rz(-pi) q[1];
rz(1.8191843) q[2];
sx q[2];
rz(-0.82652107) q[2];
sx q[2];
rz(2.7374817) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3179847) q[1];
sx q[1];
rz(-0.66901842) q[1];
sx q[1];
rz(1.4186335) q[1];
x q[2];
rz(-2.4040495) q[3];
sx q[3];
rz(-1.5182523) q[3];
sx q[3];
rz(-0.74636783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.687872) q[2];
sx q[2];
rz(-0.56679711) q[2];
sx q[2];
rz(0.68457121) q[2];
rz(-1.941393) q[3];
sx q[3];
rz(-1.129351) q[3];
sx q[3];
rz(0.17957345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0612563) q[0];
sx q[0];
rz(-2.6604524) q[0];
sx q[0];
rz(-3.1367593) q[0];
rz(1.5176516) q[1];
sx q[1];
rz(-0.74394092) q[1];
sx q[1];
rz(-0.25371107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2462354) q[0];
sx q[0];
rz(-2.9140436) q[0];
sx q[0];
rz(1.1832817) q[0];
rz(-pi) q[1];
rz(0.51400916) q[2];
sx q[2];
rz(-0.80482641) q[2];
sx q[2];
rz(2.6857306) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39832662) q[1];
sx q[1];
rz(-2.4383845) q[1];
sx q[1];
rz(1.875967) q[1];
x q[2];
rz(2.5911035) q[3];
sx q[3];
rz(-1.7409678) q[3];
sx q[3];
rz(-0.44743928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1634875) q[2];
sx q[2];
rz(-1.2349962) q[2];
sx q[2];
rz(1.182349) q[2];
rz(-2.4649418) q[3];
sx q[3];
rz(-2.7550321) q[3];
sx q[3];
rz(-2.1277908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5273298) q[0];
sx q[0];
rz(-2.3557721) q[0];
sx q[0];
rz(-2.8489805) q[0];
rz(1.0489427) q[1];
sx q[1];
rz(-1.4194149) q[1];
sx q[1];
rz(0.70379757) q[1];
rz(-1.9090777) q[2];
sx q[2];
rz(-1.1223553) q[2];
sx q[2];
rz(-1.9434402) q[2];
rz(1.8234689) q[3];
sx q[3];
rz(-2.7636486) q[3];
sx q[3];
rz(-1.6259125) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
