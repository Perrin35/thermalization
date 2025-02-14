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
rz(0.97857082) q[0];
sx q[0];
rz(-0.11712722) q[0];
sx q[0];
rz(0.22476619) q[0];
rz(-0.4966785) q[1];
sx q[1];
rz(-1.574312) q[1];
sx q[1];
rz(-1.9615251) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4311572) q[0];
sx q[0];
rz(-1.707516) q[0];
sx q[0];
rz(-1.7140634) q[0];
rz(-1.5772555) q[2];
sx q[2];
rz(-1.1273555) q[2];
sx q[2];
rz(-0.214516) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22332033) q[1];
sx q[1];
rz(-1.8448571) q[1];
sx q[1];
rz(2.5943701) q[1];
rz(0.080022607) q[3];
sx q[3];
rz(-2.0014187) q[3];
sx q[3];
rz(2.9324556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67285377) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(1.2057745) q[2];
rz(-2.2852211) q[3];
sx q[3];
rz(-0.93256336) q[3];
sx q[3];
rz(2.2642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4641651) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(2.933266) q[0];
rz(1.9147929) q[1];
sx q[1];
rz(-1.0700763) q[1];
sx q[1];
rz(2.5770381) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78075942) q[0];
sx q[0];
rz(-1.6134904) q[0];
sx q[0];
rz(-1.4111931) q[0];
x q[1];
rz(0.54148285) q[2];
sx q[2];
rz(-0.80831203) q[2];
sx q[2];
rz(1.186651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98726749) q[1];
sx q[1];
rz(-0.72817737) q[1];
sx q[1];
rz(2.8716692) q[1];
rz(-pi) q[2];
rz(0.019314841) q[3];
sx q[3];
rz(-2.1105511) q[3];
sx q[3];
rz(-3.0125953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5047001) q[2];
sx q[2];
rz(-1.8919287) q[2];
sx q[2];
rz(1.4862109) q[2];
rz(-0.49556035) q[3];
sx q[3];
rz(-1.7335745) q[3];
sx q[3];
rz(-1.0068309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9946852) q[0];
sx q[0];
rz(-0.72145307) q[0];
sx q[0];
rz(-1.7578693) q[0];
rz(1.9231298) q[1];
sx q[1];
rz(-2.0530901) q[1];
sx q[1];
rz(2.7168435) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0828511) q[0];
sx q[0];
rz(-0.070151873) q[0];
sx q[0];
rz(0.46746107) q[0];
rz(-pi) q[1];
rz(1.8239543) q[2];
sx q[2];
rz(-1.3382398) q[2];
sx q[2];
rz(-0.81040224) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.220881) q[1];
sx q[1];
rz(-1.494981) q[1];
sx q[1];
rz(-2.7644509) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5307755) q[3];
sx q[3];
rz(-2.1870625) q[3];
sx q[3];
rz(-1.54824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4212627) q[2];
sx q[2];
rz(-2.0739372) q[2];
sx q[2];
rz(3.1171411) q[2];
rz(1.1484185) q[3];
sx q[3];
rz(-2.0423856) q[3];
sx q[3];
rz(-2.059977) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352683) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(2.9277053) q[0];
rz(2.3230486) q[1];
sx q[1];
rz(-1.3571502) q[1];
sx q[1];
rz(0.68665409) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.548107) q[0];
sx q[0];
rz(-1.0229551) q[0];
sx q[0];
rz(2.0440897) q[0];
rz(-pi) q[1];
rz(1.114421) q[2];
sx q[2];
rz(-0.59774146) q[2];
sx q[2];
rz(-2.1591612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5221646) q[1];
sx q[1];
rz(-2.3261737) q[1];
sx q[1];
rz(-1.700089) q[1];
x q[2];
rz(-2.7044933) q[3];
sx q[3];
rz(-1.0954787) q[3];
sx q[3];
rz(-0.10224414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0246747) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(0.68946687) q[2];
rz(1.3245448) q[3];
sx q[3];
rz(-1.9246512) q[3];
sx q[3];
rz(-2.569258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784742) q[0];
sx q[0];
rz(-3.0480223) q[0];
sx q[0];
rz(-2.1043188) q[0];
rz(-2.4689238) q[1];
sx q[1];
rz(-1.9692407) q[1];
sx q[1];
rz(0.79576463) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8032629) q[0];
sx q[0];
rz(-1.065155) q[0];
sx q[0];
rz(-0.36257669) q[0];
x q[1];
rz(0.3032844) q[2];
sx q[2];
rz(-1.0192843) q[2];
sx q[2];
rz(1.1597275) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2356253) q[1];
sx q[1];
rz(-2.2722031) q[1];
sx q[1];
rz(0.4404621) q[1];
x q[2];
rz(-0.66821258) q[3];
sx q[3];
rz(-1.8132121) q[3];
sx q[3];
rz(2.6072864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1611288) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(-2.3058092) q[2];
rz(-0.25660822) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-0.49032828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68818727) q[0];
sx q[0];
rz(-1.1510993) q[0];
sx q[0];
rz(1.9292462) q[0];
rz(-1.0351828) q[1];
sx q[1];
rz(-1.4915219) q[1];
sx q[1];
rz(-1.0119247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31549227) q[0];
sx q[0];
rz(-1.4619215) q[0];
sx q[0];
rz(-1.7427765) q[0];
x q[1];
rz(3.1306562) q[2];
sx q[2];
rz(-1.8746398) q[2];
sx q[2];
rz(-2.2774709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3287836) q[1];
sx q[1];
rz(-1.145383) q[1];
sx q[1];
rz(2.0219457) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3009146) q[3];
sx q[3];
rz(-1.9414475) q[3];
sx q[3];
rz(2.5605367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47784352) q[2];
sx q[2];
rz(-2.9426844) q[2];
sx q[2];
rz(2.193023) q[2];
rz(0.47032022) q[3];
sx q[3];
rz(-1.0985273) q[3];
sx q[3];
rz(-1.8868014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54742852) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(-1.8773361) q[0];
rz(0.7803548) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(-0.19187127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.980968) q[0];
sx q[0];
rz(-1.1066529) q[0];
sx q[0];
rz(-1.9689409) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7816747) q[2];
sx q[2];
rz(-1.9267757) q[2];
sx q[2];
rz(2.6448665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.807232) q[1];
sx q[1];
rz(-1.8340115) q[1];
sx q[1];
rz(2.6802147) q[1];
x q[2];
rz(0.34911134) q[3];
sx q[3];
rz(-2.1421332) q[3];
sx q[3];
rz(3.0219363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4912305) q[2];
sx q[2];
rz(-2.0872842) q[2];
sx q[2];
rz(-2.3336156) q[2];
rz(-0.84336495) q[3];
sx q[3];
rz(-1.1542412) q[3];
sx q[3];
rz(1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.657044) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(-0.76883823) q[0];
rz(0.27278849) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(-1.3974894) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5834806) q[0];
sx q[0];
rz(-1.8761823) q[0];
sx q[0];
rz(-2.7598513) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0447803) q[2];
sx q[2];
rz(-2.5616841) q[2];
sx q[2];
rz(-1.2666262) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8500415) q[1];
sx q[1];
rz(-1.5102036) q[1];
sx q[1];
rz(1.4539976) q[1];
rz(2.0040183) q[3];
sx q[3];
rz(-2.1931291) q[3];
sx q[3];
rz(0.64309263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.93691319) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(-1.5020471) q[2];
rz(-1.8223193) q[3];
sx q[3];
rz(-0.84984142) q[3];
sx q[3];
rz(1.5892861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69916344) q[0];
sx q[0];
rz(-2.2971239) q[0];
sx q[0];
rz(-0.97802877) q[0];
rz(2.6507822) q[1];
sx q[1];
rz(-1.699479) q[1];
sx q[1];
rz(1.6455654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.967088) q[0];
sx q[0];
rz(-1.3347111) q[0];
sx q[0];
rz(1.4437136) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7285887) q[2];
sx q[2];
rz(-1.4133412) q[2];
sx q[2];
rz(0.26034376) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5900951) q[1];
sx q[1];
rz(-1.8935154) q[1];
sx q[1];
rz(1.8264052) q[1];
rz(-2.6096212) q[3];
sx q[3];
rz(-1.6364179) q[3];
sx q[3];
rz(-0.73871368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28828037) q[2];
sx q[2];
rz(-2.2744982) q[2];
sx q[2];
rz(1.3221928) q[2];
rz(0.36618048) q[3];
sx q[3];
rz(-1.7145232) q[3];
sx q[3];
rz(-2.0100994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91728297) q[0];
sx q[0];
rz(-1.6095105) q[0];
sx q[0];
rz(0.073534615) q[0];
rz(1.8247983) q[1];
sx q[1];
rz(-1.8837181) q[1];
sx q[1];
rz(1.8427461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.038485) q[0];
sx q[0];
rz(-1.847205) q[0];
sx q[0];
rz(-1.7544708) q[0];
x q[1];
rz(0.13330524) q[2];
sx q[2];
rz(-2.0553737) q[2];
sx q[2];
rz(-2.7078748) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3042481) q[1];
sx q[1];
rz(-0.44884822) q[1];
sx q[1];
rz(-1.0034837) q[1];
rz(-pi) q[2];
rz(3.0231944) q[3];
sx q[3];
rz(-2.3984635) q[3];
sx q[3];
rz(1.2606674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1765882) q[2];
sx q[2];
rz(-2.1829288) q[2];
sx q[2];
rz(-2.5625572) q[2];
rz(1.6252919) q[3];
sx q[3];
rz(-2.5987891) q[3];
sx q[3];
rz(1.496284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2753006) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(-0.91118377) q[1];
sx q[1];
rz(-1.3180399) q[1];
sx q[1];
rz(-1.3175189) q[1];
rz(-1.5742792) q[2];
sx q[2];
rz(-1.9949253) q[2];
sx q[2];
rz(-0.16535769) q[2];
rz(1.1796307) q[3];
sx q[3];
rz(-0.65476553) q[3];
sx q[3];
rz(-2.7162566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
