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
rz(3.717489) q[0];
sx q[0];
rz(5.1562638) q[0];
rz(-0.95070401) q[1];
sx q[1];
rz(-0.60125142) q[1];
sx q[1];
rz(-0.33369219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55328773) q[0];
sx q[0];
rz(-1.1516363) q[0];
sx q[0];
rz(-0.16417154) q[0];
rz(-pi) q[1];
rz(2.8744281) q[2];
sx q[2];
rz(-1.1295302) q[2];
sx q[2];
rz(-1.6187606) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0036949) q[1];
sx q[1];
rz(-2.7337791) q[1];
sx q[1];
rz(2.5812056) q[1];
rz(-pi) q[2];
rz(-0.80218704) q[3];
sx q[3];
rz(-1.7628551) q[3];
sx q[3];
rz(-0.036969846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2600962) q[2];
sx q[2];
rz(-2.0659955) q[2];
sx q[2];
rz(-1.0502226) q[2];
rz(2.8460734) q[3];
sx q[3];
rz(-2.3727356) q[3];
sx q[3];
rz(-1.3692921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1693901) q[0];
sx q[0];
rz(-1.0193595) q[0];
sx q[0];
rz(-2.4743359) q[0];
rz(0.15070209) q[1];
sx q[1];
rz(-1.0548016) q[1];
sx q[1];
rz(-0.31164718) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0586747) q[0];
sx q[0];
rz(-0.5265412) q[0];
sx q[0];
rz(-1.4301926) q[0];
rz(2.0097334) q[2];
sx q[2];
rz(-1.2905741) q[2];
sx q[2];
rz(0.2187905) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5087392) q[1];
sx q[1];
rz(-1.3876186) q[1];
sx q[1];
rz(-1.7126717) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4979579) q[3];
sx q[3];
rz(-1.1741116) q[3];
sx q[3];
rz(1.9997627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47505891) q[2];
sx q[2];
rz(-2.2481613) q[2];
sx q[2];
rz(2.3739572) q[2];
rz(2.660699) q[3];
sx q[3];
rz(-0.77156639) q[3];
sx q[3];
rz(-1.5868928) q[3];
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
rz(0.84592205) q[0];
sx q[0];
rz(-1.6039811) q[0];
sx q[0];
rz(0.77958244) q[0];
rz(-2.7698611) q[1];
sx q[1];
rz(-1.0075684) q[1];
sx q[1];
rz(-0.76098162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0944989) q[0];
sx q[0];
rz(-1.2505184) q[0];
sx q[0];
rz(0.30517942) q[0];
rz(-pi) q[1];
x q[1];
rz(2.089431) q[2];
sx q[2];
rz(-2.385879) q[2];
sx q[2];
rz(-1.1417127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0798988) q[1];
sx q[1];
rz(-1.736228) q[1];
sx q[1];
rz(2.295619) q[1];
rz(-pi) q[2];
rz(3.1366421) q[3];
sx q[3];
rz(-0.85270665) q[3];
sx q[3];
rz(-0.45482054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.80428213) q[2];
sx q[2];
rz(-0.90039841) q[2];
sx q[2];
rz(2.6864181) q[2];
rz(0.69194397) q[3];
sx q[3];
rz(-0.71440905) q[3];
sx q[3];
rz(2.3327995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96875018) q[0];
sx q[0];
rz(-1.7946365) q[0];
sx q[0];
rz(-2.8853048) q[0];
rz(-1.434727) q[1];
sx q[1];
rz(-1.4204357) q[1];
sx q[1];
rz(3.0228379) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2560121) q[0];
sx q[0];
rz(-0.47614723) q[0];
sx q[0];
rz(2.7976967) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2776472) q[2];
sx q[2];
rz(-1.4481232) q[2];
sx q[2];
rz(2.6495767) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.12417158) q[1];
sx q[1];
rz(-1.802449) q[1];
sx q[1];
rz(-0.79504658) q[1];
rz(1.1944653) q[3];
sx q[3];
rz(-0.85187618) q[3];
sx q[3];
rz(3.1270936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.28748301) q[2];
sx q[2];
rz(-2.867925) q[2];
sx q[2];
rz(2.0859094) q[2];
rz(-1.4500729) q[3];
sx q[3];
rz(-1.6943211) q[3];
sx q[3];
rz(-2.292574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56181041) q[0];
sx q[0];
rz(-2.9636443) q[0];
sx q[0];
rz(1.6135038) q[0];
rz(-1.1771419) q[1];
sx q[1];
rz(-1.2700932) q[1];
sx q[1];
rz(-1.5632163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3731454) q[0];
sx q[0];
rz(-2.5119436) q[0];
sx q[0];
rz(-3.0238053) q[0];
rz(-0.41069846) q[2];
sx q[2];
rz(-1.8735527) q[2];
sx q[2];
rz(-1.8524285) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10341748) q[1];
sx q[1];
rz(-0.6486054) q[1];
sx q[1];
rz(0.23596779) q[1];
rz(1.568119) q[3];
sx q[3];
rz(-1.1613881) q[3];
sx q[3];
rz(2.1684627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0057808) q[2];
sx q[2];
rz(-1.1567189) q[2];
sx q[2];
rz(-1.3021674) q[2];
rz(-2.8016413) q[3];
sx q[3];
rz(-2.3348742) q[3];
sx q[3];
rz(2.4792041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0078916773) q[0];
sx q[0];
rz(-1.5971203) q[0];
sx q[0];
rz(1.9418465) q[0];
rz(1.3613191) q[1];
sx q[1];
rz(-0.97313762) q[1];
sx q[1];
rz(0.57788411) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91980714) q[0];
sx q[0];
rz(-0.01520201) q[0];
sx q[0];
rz(2.0073246) q[0];
rz(-2.7087542) q[2];
sx q[2];
rz(-1.7146401) q[2];
sx q[2];
rz(-1.8707616) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84915811) q[1];
sx q[1];
rz(-2.6657988) q[1];
sx q[1];
rz(0.30026309) q[1];
x q[2];
rz(-1.8697) q[3];
sx q[3];
rz(-1.5953957) q[3];
sx q[3];
rz(0.92211039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63152385) q[2];
sx q[2];
rz(-2.5670467) q[2];
sx q[2];
rz(-0.4734545) q[2];
rz(-1.2724642) q[3];
sx q[3];
rz(-1.3926287) q[3];
sx q[3];
rz(0.22245358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94721395) q[0];
sx q[0];
rz(-0.65776238) q[0];
sx q[0];
rz(0.33313242) q[0];
rz(-0.22854742) q[1];
sx q[1];
rz(-1.3401778) q[1];
sx q[1];
rz(-2.5659335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4116019) q[0];
sx q[0];
rz(-1.6766929) q[0];
sx q[0];
rz(-0.81231711) q[0];
rz(-pi) q[1];
rz(-1.2387462) q[2];
sx q[2];
rz(-1.1162469) q[2];
sx q[2];
rz(-3.0177096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0488494) q[1];
sx q[1];
rz(-1.7705838) q[1];
sx q[1];
rz(1.2380935) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8262528) q[3];
sx q[3];
rz(-2.0386348) q[3];
sx q[3];
rz(-1.6088736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.26611844) q[2];
sx q[2];
rz(-0.55142752) q[2];
sx q[2];
rz(-1.4761338) q[2];
rz(1.1009781) q[3];
sx q[3];
rz(-2.0783547) q[3];
sx q[3];
rz(2.6235918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5440893) q[0];
sx q[0];
rz(-0.88633716) q[0];
sx q[0];
rz(2.8344717) q[0];
rz(2.8111474) q[1];
sx q[1];
rz(-1.5978866) q[1];
sx q[1];
rz(-1.7764567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4348286) q[0];
sx q[0];
rz(-1.4894823) q[0];
sx q[0];
rz(1.5887194) q[0];
rz(-pi) q[1];
rz(-2.6583985) q[2];
sx q[2];
rz(-0.38887923) q[2];
sx q[2];
rz(2.4318397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.04330895) q[1];
sx q[1];
rz(-1.829147) q[1];
sx q[1];
rz(-1.1275379) q[1];
rz(-1.4151814) q[3];
sx q[3];
rz(-0.74374108) q[3];
sx q[3];
rz(-1.2228325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.27552989) q[2];
sx q[2];
rz(-0.58774647) q[2];
sx q[2];
rz(-0.24869871) q[2];
rz(-1.9281049) q[3];
sx q[3];
rz(-1.4444193) q[3];
sx q[3];
rz(0.061633751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1520749) q[0];
sx q[0];
rz(-0.9062506) q[0];
sx q[0];
rz(-0.686598) q[0];
rz(1.1129414) q[1];
sx q[1];
rz(-2.3390892) q[1];
sx q[1];
rz(-0.31563219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9292727) q[0];
sx q[0];
rz(-1.2480191) q[0];
sx q[0];
rz(0.048521532) q[0];
rz(-pi) q[1];
rz(2.8807345) q[2];
sx q[2];
rz(-0.77702287) q[2];
sx q[2];
rz(2.379247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.823608) q[1];
sx q[1];
rz(-0.66901842) q[1];
sx q[1];
rz(-1.7229592) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0635438) q[3];
sx q[3];
rz(-0.73906088) q[3];
sx q[3];
rz(2.259425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.45372066) q[2];
sx q[2];
rz(-2.5747955) q[2];
sx q[2];
rz(-2.4570214) q[2];
rz(-1.2001996) q[3];
sx q[3];
rz(-2.0122416) q[3];
sx q[3];
rz(0.17957345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0612563) q[0];
sx q[0];
rz(-2.6604524) q[0];
sx q[0];
rz(-0.0048333724) q[0];
rz(1.5176516) q[1];
sx q[1];
rz(-2.3976517) q[1];
sx q[1];
rz(-2.8878816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4986202) q[0];
sx q[0];
rz(-1.7812087) q[0];
sx q[0];
rz(0.087281373) q[0];
rz(-2.0433304) q[2];
sx q[2];
rz(-0.8923549) q[2];
sx q[2];
rz(2.0022165) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1349134) q[1];
sx q[1];
rz(-0.90613922) q[1];
sx q[1];
rz(-0.24941872) q[1];
x q[2];
rz(-1.3718448) q[3];
sx q[3];
rz(-1.0291489) q[3];
sx q[3];
rz(-1.0197848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1634875) q[2];
sx q[2];
rz(-1.2349962) q[2];
sx q[2];
rz(1.9592436) q[2];
rz(-0.67665082) q[3];
sx q[3];
rz(-2.7550321) q[3];
sx q[3];
rz(2.1277908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61426281) q[0];
sx q[0];
rz(-0.78582055) q[0];
sx q[0];
rz(0.2926122) q[0];
rz(-1.0489427) q[1];
sx q[1];
rz(-1.7221778) q[1];
sx q[1];
rz(-2.4377951) q[1];
rz(2.537773) q[2];
sx q[2];
rz(-2.5868844) q[2];
sx q[2];
rz(-1.2617688) q[2];
rz(-1.3181237) q[3];
sx q[3];
rz(-2.7636486) q[3];
sx q[3];
rz(-1.6259125) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
