OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7991601) q[0];
sx q[0];
rz(-2.7437796) q[0];
sx q[0];
rz(1.2678658) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(-0.78650147) q[1];
sx q[1];
rz(1.9222577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.08377) q[0];
sx q[0];
rz(-1.738213) q[0];
sx q[0];
rz(2.1517702) q[0];
rz(-0.54852416) q[2];
sx q[2];
rz(-1.7146352) q[2];
sx q[2];
rz(-1.583969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3883512) q[1];
sx q[1];
rz(-1.8495535) q[1];
sx q[1];
rz(-2.9294347) q[1];
rz(1.8699617) q[3];
sx q[3];
rz(-2.4623021) q[3];
sx q[3];
rz(-2.3878179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5775602) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(1.3684028) q[2];
rz(2.2300301) q[3];
sx q[3];
rz(-1.3699968) q[3];
sx q[3];
rz(-0.024356775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0381222) q[0];
sx q[0];
rz(-2.7439674) q[0];
sx q[0];
rz(0.11904112) q[0];
rz(-2.0114404) q[1];
sx q[1];
rz(-2.1739013) q[1];
sx q[1];
rz(1.4281323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9790065) q[0];
sx q[0];
rz(-1.1631199) q[0];
sx q[0];
rz(-0.56637052) q[0];
rz(2.7409254) q[2];
sx q[2];
rz(-1.180449) q[2];
sx q[2];
rz(1.1125178) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2547639) q[1];
sx q[1];
rz(-0.48476754) q[1];
sx q[1];
rz(1.280275) q[1];
x q[2];
rz(-1.9088332) q[3];
sx q[3];
rz(-1.2814953) q[3];
sx q[3];
rz(-2.3551331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6156562) q[2];
sx q[2];
rz(-1.6776513) q[2];
sx q[2];
rz(2.3770135) q[2];
rz(-2.6584451) q[3];
sx q[3];
rz(-1.8254447) q[3];
sx q[3];
rz(0.84883261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0371542) q[0];
sx q[0];
rz(-0.26532441) q[0];
sx q[0];
rz(-1.4720488) q[0];
rz(2.5813591) q[1];
sx q[1];
rz(-1.3872223) q[1];
sx q[1];
rz(-1.4405506) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.773303) q[0];
sx q[0];
rz(-1.76923) q[0];
sx q[0];
rz(-1.8097359) q[0];
rz(-pi) q[1];
rz(2.4666489) q[2];
sx q[2];
rz(-1.5663354) q[2];
sx q[2];
rz(0.5677815) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8941108) q[1];
sx q[1];
rz(-1.2318582) q[1];
sx q[1];
rz(-2.1268658) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23594956) q[3];
sx q[3];
rz(-1.7914158) q[3];
sx q[3];
rz(1.5853796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6560087) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(-0.34240016) q[2];
rz(-1.7229236) q[3];
sx q[3];
rz(-2.4125621) q[3];
sx q[3];
rz(0.5595783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90958428) q[0];
sx q[0];
rz(-2.7689458) q[0];
sx q[0];
rz(1.6147856) q[0];
rz(2.3341663) q[1];
sx q[1];
rz(-1.5680983) q[1];
sx q[1];
rz(-1.2015013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75185173) q[0];
sx q[0];
rz(-2.3065615) q[0];
sx q[0];
rz(0.40919183) q[0];
rz(-0.92113481) q[2];
sx q[2];
rz(-1.0277205) q[2];
sx q[2];
rz(-1.4411826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0240942) q[1];
sx q[1];
rz(-1.9650241) q[1];
sx q[1];
rz(2.3181163) q[1];
rz(-pi) q[2];
rz(1.5189999) q[3];
sx q[3];
rz(-1.8267617) q[3];
sx q[3];
rz(-1.761375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3379007) q[2];
sx q[2];
rz(-2.1521229) q[2];
sx q[2];
rz(-1.9039512) q[2];
rz(-1.8303653) q[3];
sx q[3];
rz(-1.5179736) q[3];
sx q[3];
rz(-1.1782014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0230947) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(2.232724) q[0];
rz(2.2591649) q[1];
sx q[1];
rz(-0.71417037) q[1];
sx q[1];
rz(2.4095101) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1650989) q[0];
sx q[0];
rz(-2.0400751) q[0];
sx q[0];
rz(0.83013541) q[0];
rz(-pi) q[1];
rz(-2.3734748) q[2];
sx q[2];
rz(-1.7752194) q[2];
sx q[2];
rz(2.8724175) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.656658) q[1];
sx q[1];
rz(-1.1371433) q[1];
sx q[1];
rz(2.0051433) q[1];
rz(-0.74589269) q[3];
sx q[3];
rz(-2.159366) q[3];
sx q[3];
rz(-0.70927995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.070179209) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(1.8522235) q[2];
rz(1.6728801) q[3];
sx q[3];
rz(-1.9332935) q[3];
sx q[3];
rz(2.7220791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5759204) q[0];
sx q[0];
rz(-1.9603632) q[0];
sx q[0];
rz(2.0400203) q[0];
rz(-0.22483243) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(-0.98519957) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3168047) q[0];
sx q[0];
rz(-1.7508306) q[0];
sx q[0];
rz(0.92702528) q[0];
rz(-pi) q[1];
rz(2.5271068) q[2];
sx q[2];
rz(-1.859815) q[2];
sx q[2];
rz(-2.7151263) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26347736) q[1];
sx q[1];
rz(-2.2190071) q[1];
sx q[1];
rz(-0.99261673) q[1];
x q[2];
rz(0.52805488) q[3];
sx q[3];
rz(-1.5925358) q[3];
sx q[3];
rz(-2.4491212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3605986) q[2];
sx q[2];
rz(-0.66086078) q[2];
sx q[2];
rz(-1.9471656) q[2];
rz(-0.099695168) q[3];
sx q[3];
rz(-1.3705148) q[3];
sx q[3];
rz(0.97894198) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7128971) q[0];
sx q[0];
rz(-0.45658699) q[0];
sx q[0];
rz(2.0275443) q[0];
rz(1.3975337) q[1];
sx q[1];
rz(-1.7618529) q[1];
sx q[1];
rz(-0.31375113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.926696) q[0];
sx q[0];
rz(-0.29698661) q[0];
sx q[0];
rz(-1.7420705) q[0];
rz(2.9767562) q[2];
sx q[2];
rz(-1.4408752) q[2];
sx q[2];
rz(2.8292953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.113768) q[1];
sx q[1];
rz(-1.8701487) q[1];
sx q[1];
rz(-0.76211263) q[1];
rz(-pi) q[2];
rz(-1.3082771) q[3];
sx q[3];
rz(-2.2394925) q[3];
sx q[3];
rz(1.4895542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0745915) q[2];
sx q[2];
rz(-0.40412298) q[2];
sx q[2];
rz(2.5355549) q[2];
rz(-2.0634985) q[3];
sx q[3];
rz(-2.2793016) q[3];
sx q[3];
rz(-1.3585453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17230497) q[0];
sx q[0];
rz(-2.2242039) q[0];
sx q[0];
rz(-2.0148) q[0];
rz(2.2276095) q[1];
sx q[1];
rz(-0.4387478) q[1];
sx q[1];
rz(-0.21805683) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0261966) q[0];
sx q[0];
rz(-0.65803981) q[0];
sx q[0];
rz(2.3283919) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70773196) q[2];
sx q[2];
rz(-0.92365217) q[2];
sx q[2];
rz(0.92491481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40849028) q[1];
sx q[1];
rz(-2.4839851) q[1];
sx q[1];
rz(-1.3642743) q[1];
x q[2];
rz(0.117649) q[3];
sx q[3];
rz(-1.434064) q[3];
sx q[3];
rz(0.6499539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5128936) q[2];
sx q[2];
rz(-2.4662374) q[2];
sx q[2];
rz(-0.2891573) q[2];
rz(-2.1469927) q[3];
sx q[3];
rz(-1.1887487) q[3];
sx q[3];
rz(0.4579671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7116123) q[0];
sx q[0];
rz(-1.0884322) q[0];
sx q[0];
rz(0.76643884) q[0];
rz(-1.6038766) q[1];
sx q[1];
rz(-1.8285373) q[1];
sx q[1];
rz(2.2235353) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34167463) q[0];
sx q[0];
rz(-0.8091439) q[0];
sx q[0];
rz(2.2412712) q[0];
rz(0.70785825) q[2];
sx q[2];
rz(-2.9212657) q[2];
sx q[2];
rz(0.79263055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5919216) q[1];
sx q[1];
rz(-1.3240314) q[1];
sx q[1];
rz(0.73149577) q[1];
rz(1.9647801) q[3];
sx q[3];
rz(-2.8716762) q[3];
sx q[3];
rz(0.22701015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90091577) q[2];
sx q[2];
rz(-2.675481) q[2];
sx q[2];
rz(3.0898928) q[2];
rz(2.2895571) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(0.26028546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8155415) q[0];
sx q[0];
rz(-1.6764078) q[0];
sx q[0];
rz(3.0503804) q[0];
rz(-0.7971898) q[1];
sx q[1];
rz(-1.7557095) q[1];
sx q[1];
rz(2.7104654) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53724683) q[0];
sx q[0];
rz(-1.6291999) q[0];
sx q[0];
rz(1.6184774) q[0];
x q[1];
rz(0.92208879) q[2];
sx q[2];
rz(-1.6493634) q[2];
sx q[2];
rz(2.9804413) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9111002) q[1];
sx q[1];
rz(-0.5020895) q[1];
sx q[1];
rz(2.6761495) q[1];
rz(-pi) q[2];
rz(-0.96980181) q[3];
sx q[3];
rz(-1.2088115) q[3];
sx q[3];
rz(-0.34011823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5690696) q[2];
sx q[2];
rz(-1.5670245) q[2];
sx q[2];
rz(-1.8756867) q[2];
rz(-1.2931394) q[3];
sx q[3];
rz(-1.061941) q[3];
sx q[3];
rz(2.7354447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010392808) q[0];
sx q[0];
rz(-0.69080234) q[0];
sx q[0];
rz(0.77558415) q[0];
rz(-2.5776183) q[1];
sx q[1];
rz(-2.0717944) q[1];
sx q[1];
rz(2.0800128) q[1];
rz(1.3971267) q[2];
sx q[2];
rz(-0.43890719) q[2];
sx q[2];
rz(-0.70499805) q[2];
rz(2.8108786) q[3];
sx q[3];
rz(-2.0277818) q[3];
sx q[3];
rz(-0.51206931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
