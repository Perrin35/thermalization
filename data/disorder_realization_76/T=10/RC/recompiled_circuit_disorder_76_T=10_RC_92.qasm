OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8653712) q[0];
sx q[0];
rz(-2.2844391) q[0];
sx q[0];
rz(3.0091118) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(5.6981882) q[1];
sx q[1];
rz(11.873801) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8938457) q[0];
sx q[0];
rz(-0.5286628) q[0];
sx q[0];
rz(1.7696487) q[0];
x q[1];
rz(-2.4970826) q[2];
sx q[2];
rz(-1.1877726) q[2];
sx q[2];
rz(0.884998) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3092279) q[1];
sx q[1];
rz(-1.698306) q[1];
sx q[1];
rz(0.072673701) q[1];
x q[2];
rz(1.1107221) q[3];
sx q[3];
rz(-2.6220136) q[3];
sx q[3];
rz(0.67668623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0341558) q[2];
sx q[2];
rz(-0.52705708) q[2];
sx q[2];
rz(1.5365323) q[2];
rz(-1.6202554) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261616) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(-1.2492299) q[0];
rz(2.5800887) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(2.5610279) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97894325) q[0];
sx q[0];
rz(-1.2622841) q[0];
sx q[0];
rz(-1.4573775) q[0];
rz(1.7891907) q[2];
sx q[2];
rz(-0.94422715) q[2];
sx q[2];
rz(2.7004957) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.7564056) q[1];
sx q[1];
rz(-2.2573973) q[1];
sx q[1];
rz(-2.3910206) q[1];
rz(-pi) q[2];
rz(-2.21653) q[3];
sx q[3];
rz(-2.5930773) q[3];
sx q[3];
rz(0.95228449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7291752) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(0.87835971) q[2];
rz(2.7495524) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(2.5382606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47675258) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-2.3679249) q[0];
rz(-0.0013008612) q[1];
sx q[1];
rz(-1.5258077) q[1];
sx q[1];
rz(-0.032827854) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5531909) q[0];
sx q[0];
rz(-1.8817888) q[0];
sx q[0];
rz(-1.3283967) q[0];
rz(-1.1187395) q[2];
sx q[2];
rz(-2.5411122) q[2];
sx q[2];
rz(0.26417363) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5579538) q[1];
sx q[1];
rz(-1.7350405) q[1];
sx q[1];
rz(0.79973952) q[1];
x q[2];
rz(1.0266617) q[3];
sx q[3];
rz(-1.2715724) q[3];
sx q[3];
rz(-0.92218375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34439987) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(-0.27734217) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(0.69916454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4375777) q[0];
sx q[0];
rz(-2.8058348) q[0];
sx q[0];
rz(2.8787956) q[0];
rz(-0.2335877) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-2.3707726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.748427) q[0];
sx q[0];
rz(-1.59032) q[0];
sx q[0];
rz(0.016419134) q[0];
x q[1];
rz(0.1044354) q[2];
sx q[2];
rz(-1.0829139) q[2];
sx q[2];
rz(-1.7761531) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1781436) q[1];
sx q[1];
rz(-1.932486) q[1];
sx q[1];
rz(-0.55808918) q[1];
rz(-pi) q[2];
rz(-0.869107) q[3];
sx q[3];
rz(-2.8512555) q[3];
sx q[3];
rz(-1.3645736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(0.7129933) q[2];
rz(-1.0130079) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(2.1876984) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35904303) q[0];
sx q[0];
rz(-1.0721711) q[0];
sx q[0];
rz(-1.7011401) q[0];
rz(-0.094093181) q[1];
sx q[1];
rz(-2.4021939) q[1];
sx q[1];
rz(-2.9715911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3225587) q[0];
sx q[0];
rz(-1.337262) q[0];
sx q[0];
rz(-0.76604953) q[0];
x q[1];
rz(-1.1414358) q[2];
sx q[2];
rz(-2.6300207) q[2];
sx q[2];
rz(-2.9432952) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1483037) q[1];
sx q[1];
rz(-1.171604) q[1];
sx q[1];
rz(0.94435512) q[1];
x q[2];
rz(0.65808987) q[3];
sx q[3];
rz(-0.98539017) q[3];
sx q[3];
rz(-2.493849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(-1.4040995) q[2];
rz(-1.5480301) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(-1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(0.45853841) q[0];
rz(2.8857152) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(-2.4564254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.348649) q[0];
sx q[0];
rz(-2.2198338) q[0];
sx q[0];
rz(3.0543324) q[0];
x q[1];
rz(1.4267169) q[2];
sx q[2];
rz(-2.3428168) q[2];
sx q[2];
rz(3.0720403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31720686) q[1];
sx q[1];
rz(-0.82815352) q[1];
sx q[1];
rz(-2.2627027) q[1];
rz(-2.1567247) q[3];
sx q[3];
rz(-1.5555447) q[3];
sx q[3];
rz(1.0451942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7489862) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(-1.4292599) q[2];
rz(2.0424992) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(-2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(0.72934735) q[0];
rz(-2.8485281) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(-1.9940631) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4097737) q[0];
sx q[0];
rz(-1.6454576) q[0];
sx q[0];
rz(2.7274107) q[0];
x q[1];
rz(-2.7717934) q[2];
sx q[2];
rz(-2.1714006) q[2];
sx q[2];
rz(-2.415654) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8116248) q[1];
sx q[1];
rz(-2.6553272) q[1];
sx q[1];
rz(-0.087841308) q[1];
rz(-pi) q[2];
rz(-0.64671867) q[3];
sx q[3];
rz(-2.5463856) q[3];
sx q[3];
rz(1.7531542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78836936) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(-0.74907556) q[2];
rz(-0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1483243) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(-1.3759026) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-0.39852279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24005323) q[0];
sx q[0];
rz(-1.1050637) q[0];
sx q[0];
rz(-0.17850152) q[0];
rz(-pi) q[1];
rz(-2.987791) q[2];
sx q[2];
rz(-2.5349344) q[2];
sx q[2];
rz(2.3369044) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9498082) q[1];
sx q[1];
rz(-0.59760082) q[1];
sx q[1];
rz(-1.6114086) q[1];
rz(-pi) q[2];
rz(-1.0141482) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(2.2198912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9514256) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(-1.8438967) q[2];
rz(-2.0166345) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(-2.4979533) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012638906) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(-1.7522316) q[0];
rz(-1.6268436) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(1.0983889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5441355) q[0];
sx q[0];
rz(-1.3381334) q[0];
sx q[0];
rz(1.0730037) q[0];
rz(-pi) q[1];
rz(2.3214925) q[2];
sx q[2];
rz(-1.2250049) q[2];
sx q[2];
rz(1.0320013) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2201982) q[1];
sx q[1];
rz(-2.1500593) q[1];
sx q[1];
rz(-1.5701576) q[1];
rz(1.7421726) q[3];
sx q[3];
rz(-1.3203353) q[3];
sx q[3];
rz(-0.40303883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2417458) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(-0.51952726) q[2];
rz(-1.8064921) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7463503) q[0];
sx q[0];
rz(-2.0210176) q[0];
sx q[0];
rz(-0.46646068) q[0];
rz(-0.17164224) q[1];
sx q[1];
rz(-1.9263093) q[1];
sx q[1];
rz(2.5126273) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0590203) q[0];
sx q[0];
rz(-1.9452381) q[0];
sx q[0];
rz(-3.0890205) q[0];
rz(-pi) q[1];
rz(-0.92832698) q[2];
sx q[2];
rz(-1.7648342) q[2];
sx q[2];
rz(2.5460668) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.92470691) q[1];
sx q[1];
rz(-2.6994884) q[1];
sx q[1];
rz(3.1367723) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4217989) q[3];
sx q[3];
rz(-1.3070953) q[3];
sx q[3];
rz(1.3849474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(1.0591327) q[2];
rz(-2.4641666) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(-2.2476851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28329904) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(2.8876866) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(-1.5268486) q[2];
sx q[2];
rz(-2.2938001) q[2];
sx q[2];
rz(-3.0707035) q[2];
rz(-1.4678636) q[3];
sx q[3];
rz(-2.5522305) q[3];
sx q[3];
rz(1.1141368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];