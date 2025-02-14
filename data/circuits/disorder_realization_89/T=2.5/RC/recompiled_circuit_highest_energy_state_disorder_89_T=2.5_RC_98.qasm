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
rz(-2.9168265) q[0];
rz(2.6449142) q[1];
sx q[1];
rz(4.7159046) q[1];
sx q[1];
rz(11.386303) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5242033) q[0];
sx q[0];
rz(-2.943882) q[0];
sx q[0];
rz(2.3376505) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5772555) q[2];
sx q[2];
rz(-1.1273555) q[2];
sx q[2];
rz(2.9270767) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6306873) q[1];
sx q[1];
rz(-2.0954335) q[1];
sx q[1];
rz(1.8888297) q[1];
rz(-pi) q[2];
rz(0.080022607) q[3];
sx q[3];
rz(-1.140174) q[3];
sx q[3];
rz(0.20913707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67285377) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(1.9358181) q[2];
rz(-0.85637158) q[3];
sx q[3];
rz(-0.93256336) q[3];
sx q[3];
rz(-2.2642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641651) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(-2.933266) q[0];
rz(-1.2267998) q[1];
sx q[1];
rz(-2.0715163) q[1];
sx q[1];
rz(0.56455451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78316654) q[0];
sx q[0];
rz(-1.7302528) q[0];
sx q[0];
rz(-0.043242976) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0759841) q[2];
sx q[2];
rz(-0.90246614) q[2];
sx q[2];
rz(-0.47023103) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7617088) q[1];
sx q[1];
rz(-1.7492047) q[1];
sx q[1];
rz(0.70990035) q[1];
x q[2];
rz(-0.019314841) q[3];
sx q[3];
rz(-1.0310415) q[3];
sx q[3];
rz(-3.0125953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5047001) q[2];
sx q[2];
rz(-1.8919287) q[2];
sx q[2];
rz(1.4862109) q[2];
rz(-0.49556035) q[3];
sx q[3];
rz(-1.4080181) q[3];
sx q[3];
rz(1.0068309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14690742) q[0];
sx q[0];
rz(-0.72145307) q[0];
sx q[0];
rz(1.7578693) q[0];
rz(-1.9231298) q[1];
sx q[1];
rz(-2.0530901) q[1];
sx q[1];
rz(0.4247492) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1871763) q[0];
sx q[0];
rz(-1.5392051) q[0];
sx q[0];
rz(0.0626465) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3282503) q[2];
sx q[2];
rz(-2.7995297) q[2];
sx q[2];
rz(-1.4881211) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.220881) q[1];
sx q[1];
rz(-1.494981) q[1];
sx q[1];
rz(2.7644509) q[1];
rz(-0.056428595) q[3];
sx q[3];
rz(-0.61739576) q[3];
sx q[3];
rz(1.4790725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4212627) q[2];
sx q[2];
rz(-1.0676554) q[2];
sx q[2];
rz(-0.024451582) q[2];
rz(1.9931741) q[3];
sx q[3];
rz(-1.0992071) q[3];
sx q[3];
rz(1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20632437) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(-2.9277053) q[0];
rz(2.3230486) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(-0.68665409) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4249466) q[0];
sx q[0];
rz(-1.1711297) q[0];
sx q[0];
rz(-0.60092385) q[0];
x q[1];
rz(0.29149194) q[2];
sx q[2];
rz(-1.0412058) q[2];
sx q[2];
rz(-1.6232217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.80685341) q[1];
sx q[1];
rz(-2.377393) q[1];
sx q[1];
rz(0.13607009) q[1];
rz(-pi) q[2];
rz(-2.0874168) q[3];
sx q[3];
rz(-1.1849019) q[3];
sx q[3];
rz(-1.8836881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0246747) q[2];
sx q[2];
rz(-1.4093829) q[2];
sx q[2];
rz(-2.4521258) q[2];
rz(-1.8170478) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(2.569258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.4631184) q[0];
sx q[0];
rz(-0.093570396) q[0];
sx q[0];
rz(-1.0372739) q[0];
rz(-2.4689238) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(-0.79576463) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.138863) q[0];
sx q[0];
rz(-0.61289591) q[0];
sx q[0];
rz(-1.0010368) q[0];
x q[1];
rz(2.8383083) q[2];
sx q[2];
rz(-1.0192843) q[2];
sx q[2];
rz(-1.1597275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5368294) q[1];
sx q[1];
rz(-0.80789551) q[1];
sx q[1];
rz(2.038234) q[1];
rz(-pi) q[2];
rz(1.2656059) q[3];
sx q[3];
rz(-0.92545907) q[3];
sx q[3];
rz(0.84922817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98046389) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(-2.3058092) q[2];
rz(2.8849844) q[3];
sx q[3];
rz(-1.1107239) q[3];
sx q[3];
rz(0.49032828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68818727) q[0];
sx q[0];
rz(-1.1510993) q[0];
sx q[0];
rz(-1.2123464) q[0];
rz(-1.0351828) q[1];
sx q[1];
rz(-1.4915219) q[1];
sx q[1];
rz(-1.0119247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4453011) q[0];
sx q[0];
rz(-2.9383349) q[0];
sx q[0];
rz(2.1392034) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6056608) q[2];
sx q[2];
rz(-2.8375585) q[2];
sx q[2];
rz(-0.82758289) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1862667) q[1];
sx q[1];
rz(-1.9791934) q[1];
sx q[1];
rz(-0.46640654) q[1];
rz(-2.5402734) q[3];
sx q[3];
rz(-0.45479247) q[3];
sx q[3];
rz(-3.0704344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5941641) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(-1.2642566) q[0];
rz(2.3612379) q[1];
sx q[1];
rz(-2.5006313) q[1];
sx q[1];
rz(-0.19187127) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3656552) q[0];
sx q[0];
rz(-1.2167551) q[0];
sx q[0];
rz(2.6440622) q[0];
rz(-pi) q[1];
rz(0.35991798) q[2];
sx q[2];
rz(-1.214817) q[2];
sx q[2];
rz(0.49672613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75429186) q[1];
sx q[1];
rz(-0.52644345) q[1];
sx q[1];
rz(-2.597288) q[1];
x q[2];
rz(-1.0818008) q[3];
sx q[3];
rz(-2.4823175) q[3];
sx q[3];
rz(-0.47283543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65036217) q[2];
sx q[2];
rz(-1.0543084) q[2];
sx q[2];
rz(-0.80797705) q[2];
rz(-2.2982277) q[3];
sx q[3];
rz(-1.9873514) q[3];
sx q[3];
rz(1.5850867) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.657044) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(2.3727544) q[0];
rz(2.8688042) q[1];
sx q[1];
rz(-1.5130679) q[1];
sx q[1];
rz(1.7441033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63015712) q[0];
sx q[0];
rz(-0.4841329) q[0];
sx q[0];
rz(-0.70229437) q[0];
rz(-2.0862085) q[2];
sx q[2];
rz(-1.8495108) q[2];
sx q[2];
rz(2.9936522) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2721401) q[1];
sx q[1];
rz(-1.6873797) q[1];
sx q[1];
rz(-3.0805853) q[1];
rz(-pi) q[2];
rz(0.66889735) q[3];
sx q[3];
rz(-1.9188768) q[3];
sx q[3];
rz(0.66434464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2046795) q[2];
sx q[2];
rz(-1.0656554) q[2];
sx q[2];
rz(-1.5020471) q[2];
rz(1.8223193) q[3];
sx q[3];
rz(-0.84984142) q[3];
sx q[3];
rz(1.5523065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.4424292) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(-2.1635639) q[0];
rz(2.6507822) q[1];
sx q[1];
rz(-1.4421137) q[1];
sx q[1];
rz(1.4960272) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.467062) q[0];
sx q[0];
rz(-0.26755324) q[0];
sx q[0];
rz(-2.6567099) q[0];
x q[1];
rz(-2.7649058) q[2];
sx q[2];
rz(-0.44038195) q[2];
sx q[2];
rz(2.1747766) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5514976) q[1];
sx q[1];
rz(-1.2480772) q[1];
sx q[1];
rz(-1.3151874) q[1];
rz(0.12883993) q[3];
sx q[3];
rz(-2.6059756) q[3];
sx q[3];
rz(-2.420466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8533123) q[2];
sx q[2];
rz(-2.2744982) q[2];
sx q[2];
rz(-1.8193998) q[2];
rz(0.36618048) q[3];
sx q[3];
rz(-1.4270695) q[3];
sx q[3];
rz(2.0100994) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91728297) q[0];
sx q[0];
rz(-1.5320822) q[0];
sx q[0];
rz(0.073534615) q[0];
rz(1.8247983) q[1];
sx q[1];
rz(-1.8837181) q[1];
sx q[1];
rz(-1.2988466) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.038485) q[0];
sx q[0];
rz(-1.2943876) q[0];
sx q[0];
rz(-1.3871219) q[0];
rz(-pi) q[1];
rz(2.0590563) q[2];
sx q[2];
rz(-1.6886782) q[2];
sx q[2];
rz(-2.0669017) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6886814) q[1];
sx q[1];
rz(-1.1961403) q[1];
sx q[1];
rz(-0.25325798) q[1];
rz(-pi) q[2];
rz(1.462684) q[3];
sx q[3];
rz(-0.834081) q[3];
sx q[3];
rz(-2.0410868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9650044) q[2];
sx q[2];
rz(-0.95866385) q[2];
sx q[2];
rz(-2.5625572) q[2];
rz(-1.5163007) q[3];
sx q[3];
rz(-0.54280353) q[3];
sx q[3];
rz(1.6453086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2753006) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(2.2304089) q[1];
sx q[1];
rz(-1.3180399) q[1];
sx q[1];
rz(-1.3175189) q[1];
rz(-1.5673135) q[2];
sx q[2];
rz(-1.1466673) q[2];
sx q[2];
rz(2.976235) q[2];
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
