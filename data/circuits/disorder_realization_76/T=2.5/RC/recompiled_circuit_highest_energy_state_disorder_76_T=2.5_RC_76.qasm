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
rz(0.87585706) q[0];
sx q[0];
rz(2.1729204) q[0];
sx q[0];
rz(7.2090413) q[0];
rz(-2.5246188) q[1];
sx q[1];
rz(-2.4763835) q[1];
sx q[1];
rz(1.8170504) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5228341) q[0];
sx q[0];
rz(-1.5069739) q[0];
sx q[0];
rz(2.5413949) q[0];
rz(-pi) q[1];
rz(-1.6683031) q[2];
sx q[2];
rz(-0.98819369) q[2];
sx q[2];
rz(-2.1511457) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.263773) q[1];
sx q[1];
rz(-0.78697453) q[1];
sx q[1];
rz(0.67194463) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4670228) q[3];
sx q[3];
rz(-0.88817443) q[3];
sx q[3];
rz(0.6011403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2609743) q[2];
sx q[2];
rz(-1.3857434) q[2];
sx q[2];
rz(-0.79208881) q[2];
rz(-0.875862) q[3];
sx q[3];
rz(-3.0328817) q[3];
sx q[3];
rz(-0.66332269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.51139128) q[0];
sx q[0];
rz(-1.8237317) q[0];
sx q[0];
rz(0.9414916) q[0];
rz(-2.1425653) q[1];
sx q[1];
rz(-2.2143054) q[1];
sx q[1];
rz(0.082854465) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64697826) q[0];
sx q[0];
rz(-2.4002693) q[0];
sx q[0];
rz(-2.9927918) q[0];
rz(2.8993334) q[2];
sx q[2];
rz(-1.9177116) q[2];
sx q[2];
rz(1.0936979) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1260316) q[1];
sx q[1];
rz(-2.0210316) q[1];
sx q[1];
rz(2.5712988) q[1];
rz(-pi) q[2];
rz(-2.1762455) q[3];
sx q[3];
rz(-1.0167398) q[3];
sx q[3];
rz(0.29747552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89722172) q[2];
sx q[2];
rz(-1.7873849) q[2];
sx q[2];
rz(-1.2574035) q[2];
rz(0.8463549) q[3];
sx q[3];
rz(-2.9823163) q[3];
sx q[3];
rz(0.67811051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1988679) q[0];
sx q[0];
rz(-1.4582448) q[0];
sx q[0];
rz(2.9123836) q[0];
rz(0.21408679) q[1];
sx q[1];
rz(-2.2702718) q[1];
sx q[1];
rz(0.3814989) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70951916) q[0];
sx q[0];
rz(-2.1901178) q[0];
sx q[0];
rz(-2.0053869) q[0];
rz(-pi) q[1];
rz(-0.83865954) q[2];
sx q[2];
rz(-1.4268954) q[2];
sx q[2];
rz(-0.14140192) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22893045) q[1];
sx q[1];
rz(-3.0813363) q[1];
sx q[1];
rz(1.0723537) q[1];
x q[2];
rz(-2.7945064) q[3];
sx q[3];
rz(-1.5659837) q[3];
sx q[3];
rz(-0.50932415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7233589) q[2];
sx q[2];
rz(-0.25622076) q[2];
sx q[2];
rz(2.5733433) q[2];
rz(-1.4189643) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(2.0645781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-0.028246183) q[0];
sx q[0];
rz(-2.009511) q[0];
sx q[0];
rz(2.8908253) q[0];
rz(0.084550683) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(-1.9409404) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6687209) q[0];
sx q[0];
rz(-0.78470147) q[0];
sx q[0];
rz(-2.0913823) q[0];
rz(0.24742207) q[2];
sx q[2];
rz(-1.1411925) q[2];
sx q[2];
rz(-1.0216433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29087092) q[1];
sx q[1];
rz(-0.71858908) q[1];
sx q[1];
rz(-0.1512089) q[1];
x q[2];
rz(-2.0968998) q[3];
sx q[3];
rz(-1.8013489) q[3];
sx q[3];
rz(-1.5073204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4667929) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(-1.8505081) q[2];
rz(-0.42168266) q[3];
sx q[3];
rz(-0.45160523) q[3];
sx q[3];
rz(1.5765367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52294937) q[0];
sx q[0];
rz(-0.72294253) q[0];
sx q[0];
rz(-1.500754) q[0];
rz(0.067642637) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(-2.0134451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65928946) q[0];
sx q[0];
rz(-1.3829675) q[0];
sx q[0];
rz(1.5330305) q[0];
x q[1];
rz(1.0447776) q[2];
sx q[2];
rz(-1.8785718) q[2];
sx q[2];
rz(0.36164944) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0813898) q[1];
sx q[1];
rz(-3.0160025) q[1];
sx q[1];
rz(-2.0141898) q[1];
rz(-pi) q[2];
rz(1.2913843) q[3];
sx q[3];
rz(-1.7959204) q[3];
sx q[3];
rz(-1.0042508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.03380123) q[2];
sx q[2];
rz(-2.0267603) q[2];
sx q[2];
rz(2.4647253) q[2];
rz(0.90773165) q[3];
sx q[3];
rz(-1.2264484) q[3];
sx q[3];
rz(-0.41456732) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229729) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(-0.76939097) q[0];
rz(1.2409302) q[1];
sx q[1];
rz(-2.2878094) q[1];
sx q[1];
rz(2.9262537) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071393911) q[0];
sx q[0];
rz(-2.5339273) q[0];
sx q[0];
rz(-1.2383226) q[0];
rz(-3.1245232) q[2];
sx q[2];
rz(-1.9015549) q[2];
sx q[2];
rz(2.0278553) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0907545) q[1];
sx q[1];
rz(-1.3267446) q[1];
sx q[1];
rz(-0.65907101) q[1];
rz(-pi) q[2];
rz(-0.69778473) q[3];
sx q[3];
rz(-1.8337733) q[3];
sx q[3];
rz(-1.3828204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.874959) q[2];
sx q[2];
rz(-1.5077488) q[2];
sx q[2];
rz(0.61666644) q[2];
rz(-2.941361) q[3];
sx q[3];
rz(-0.73035208) q[3];
sx q[3];
rz(-2.0088137) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1277593) q[0];
sx q[0];
rz(-0.56202373) q[0];
sx q[0];
rz(3.1264937) q[0];
rz(3.0191782) q[1];
sx q[1];
rz(-1.2950803) q[1];
sx q[1];
rz(1.0282358) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434179) q[0];
sx q[0];
rz(-0.91599303) q[0];
sx q[0];
rz(-1.5006939) q[0];
x q[1];
rz(-2.188855) q[2];
sx q[2];
rz(-1.3985123) q[2];
sx q[2];
rz(3.1069482) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15785698) q[1];
sx q[1];
rz(-2.7089556) q[1];
sx q[1];
rz(3.0745201) q[1];
rz(-pi) q[2];
rz(1.3002197) q[3];
sx q[3];
rz(-1.8755873) q[3];
sx q[3];
rz(-0.51606015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.60468173) q[2];
sx q[2];
rz(-1.9313507) q[2];
sx q[2];
rz(-2.6417285) q[2];
rz(0.31575051) q[3];
sx q[3];
rz(-0.67086589) q[3];
sx q[3];
rz(2.1226814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4278118) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(0.94386238) q[0];
rz(1.7367412) q[1];
sx q[1];
rz(-1.2662788) q[1];
sx q[1];
rz(-2.2199383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93225828) q[0];
sx q[0];
rz(-2.8186322) q[0];
sx q[0];
rz(-2.9684116) q[0];
rz(2.2711742) q[2];
sx q[2];
rz(-0.38258115) q[2];
sx q[2];
rz(-2.062881) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5068622) q[1];
sx q[1];
rz(-0.94351879) q[1];
sx q[1];
rz(2.2192713) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1373164) q[3];
sx q[3];
rz(-1.0579526) q[3];
sx q[3];
rz(0.23425929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1450119) q[2];
sx q[2];
rz(-1.8700446) q[2];
sx q[2];
rz(-1.5550295) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.8127706) q[3];
sx q[3];
rz(1.4012977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9911875) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(-2.4990668) q[0];
rz(-1.7036899) q[1];
sx q[1];
rz(-0.75690401) q[1];
sx q[1];
rz(-0.14370758) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7449581) q[0];
sx q[0];
rz(-1.8734583) q[0];
sx q[0];
rz(-0.58505262) q[0];
rz(-pi) q[1];
x q[1];
rz(0.917786) q[2];
sx q[2];
rz(-2.3812903) q[2];
sx q[2];
rz(1.1117293) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9160812) q[1];
sx q[1];
rz(-0.71174946) q[1];
sx q[1];
rz(-0.41081482) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1863352) q[3];
sx q[3];
rz(-1.9805621) q[3];
sx q[3];
rz(1.0863956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6045427) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(2.878888) q[2];
rz(-0.25775868) q[3];
sx q[3];
rz(-0.38193211) q[3];
sx q[3];
rz(0.8530544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2951374) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(-1.9211796) q[0];
rz(2.659761) q[1];
sx q[1];
rz(-0.82492963) q[1];
sx q[1];
rz(-2.2442472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7677444) q[0];
sx q[0];
rz(-1.4486827) q[0];
sx q[0];
rz(-2.6227345) q[0];
rz(-pi) q[1];
rz(1.6807031) q[2];
sx q[2];
rz(-1.8898003) q[2];
sx q[2];
rz(-1.6750592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5612302) q[1];
sx q[1];
rz(-1.8288611) q[1];
sx q[1];
rz(0.051405426) q[1];
rz(-pi) q[2];
rz(-1.4225716) q[3];
sx q[3];
rz(-0.61398849) q[3];
sx q[3];
rz(2.2999291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64409488) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(-3.0268055) q[2];
rz(-2.3980906) q[3];
sx q[3];
rz(-1.2657575) q[3];
sx q[3];
rz(-2.6928597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262309) q[0];
sx q[0];
rz(-2.3804433) q[0];
sx q[0];
rz(-0.95809715) q[0];
rz(1.505898) q[1];
sx q[1];
rz(-2.0151357) q[1];
sx q[1];
rz(-1.9912079) q[1];
rz(2.220357) q[2];
sx q[2];
rz(-2.2050646) q[2];
sx q[2];
rz(0.83940432) q[2];
rz(-2.791849) q[3];
sx q[3];
rz(-2.6833224) q[3];
sx q[3];
rz(0.53862017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
