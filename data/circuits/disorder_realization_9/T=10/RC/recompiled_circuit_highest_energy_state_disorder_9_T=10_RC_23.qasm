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
rz(0.23586805) q[0];
sx q[0];
rz(1.9960825) q[0];
sx q[0];
rz(9.3762015) q[0];
rz(1.5161169) q[1];
sx q[1];
rz(-2.0487509) q[1];
sx q[1];
rz(-0.71576524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9484723) q[0];
sx q[0];
rz(-0.17331757) q[0];
sx q[0];
rz(1.3993708) q[0];
x q[1];
rz(0.94004905) q[2];
sx q[2];
rz(-2.7995281) q[2];
sx q[2];
rz(-0.91724455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44385829) q[1];
sx q[1];
rz(-2.3045973) q[1];
sx q[1];
rz(-0.55056449) q[1];
rz(-pi) q[2];
rz(1.9747026) q[3];
sx q[3];
rz(-2.0527186) q[3];
sx q[3];
rz(-1.9306077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6526661) q[2];
sx q[2];
rz(-0.55997866) q[2];
sx q[2];
rz(1.0555335) q[2];
rz(-0.40600285) q[3];
sx q[3];
rz(-1.3136274) q[3];
sx q[3];
rz(2.3199911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838298) q[0];
sx q[0];
rz(-1.0394179) q[0];
sx q[0];
rz(-2.5854172) q[0];
rz(1.3293386) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(-0.083273085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8502959) q[0];
sx q[0];
rz(-1.2477739) q[0];
sx q[0];
rz(2.4730541) q[0];
rz(-pi) q[1];
rz(-2.3131392) q[2];
sx q[2];
rz(-2.3745685) q[2];
sx q[2];
rz(-1.3924862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.00035380414) q[1];
sx q[1];
rz(-1.7906931) q[1];
sx q[1];
rz(-2.6882437) q[1];
rz(-pi) q[2];
rz(2.2733685) q[3];
sx q[3];
rz(-0.69528841) q[3];
sx q[3];
rz(-0.27588683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1978153) q[2];
sx q[2];
rz(-1.752172) q[2];
sx q[2];
rz(2.6488292) q[2];
rz(-1.0265776) q[3];
sx q[3];
rz(-0.30288282) q[3];
sx q[3];
rz(1.4125642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3781085) q[0];
sx q[0];
rz(-2.5788827) q[0];
sx q[0];
rz(-0.43877959) q[0];
rz(1.6890866) q[1];
sx q[1];
rz(-2.8197598) q[1];
sx q[1];
rz(-0.39295331) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2344491) q[0];
sx q[0];
rz(-2.6496072) q[0];
sx q[0];
rz(1.6299295) q[0];
x q[1];
rz(1.5147171) q[2];
sx q[2];
rz(-2.6125557) q[2];
sx q[2];
rz(-2.7521486) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4247228) q[1];
sx q[1];
rz(-0.42415127) q[1];
sx q[1];
rz(-1.7440026) q[1];
rz(-pi) q[2];
rz(2.4923101) q[3];
sx q[3];
rz(-1.6341795) q[3];
sx q[3];
rz(1.5017623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68985933) q[2];
sx q[2];
rz(-2.2971575) q[2];
sx q[2];
rz(-2.8680657) q[2];
rz(2.5290329) q[3];
sx q[3];
rz(-1.7275683) q[3];
sx q[3];
rz(2.2851022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81371152) q[0];
sx q[0];
rz(-2.5484945) q[0];
sx q[0];
rz(0.80075049) q[0];
rz(2.4294991) q[1];
sx q[1];
rz(-2.021603) q[1];
sx q[1];
rz(-0.48316479) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8022435) q[0];
sx q[0];
rz(-0.74269811) q[0];
sx q[0];
rz(-0.00089641103) q[0];
rz(2.0676548) q[2];
sx q[2];
rz(-2.0570847) q[2];
sx q[2];
rz(0.72975791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9312716) q[1];
sx q[1];
rz(-0.64669791) q[1];
sx q[1];
rz(-3.0937322) q[1];
rz(-pi) q[2];
rz(-0.68928257) q[3];
sx q[3];
rz(-1.3069671) q[3];
sx q[3];
rz(1.0466913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98372769) q[2];
sx q[2];
rz(-0.71061504) q[2];
sx q[2];
rz(-1.857081) q[2];
rz(0.48786783) q[3];
sx q[3];
rz(-1.7043461) q[3];
sx q[3];
rz(1.6871066) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81216413) q[0];
sx q[0];
rz(-0.85947961) q[0];
sx q[0];
rz(0.46464768) q[0];
rz(-2.1931785) q[1];
sx q[1];
rz(-2.3034838) q[1];
sx q[1];
rz(1.4533739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8226979) q[0];
sx q[0];
rz(-1.5359211) q[0];
sx q[0];
rz(-0.025392763) q[0];
rz(-2.4011432) q[2];
sx q[2];
rz(-1.3794823) q[2];
sx q[2];
rz(-0.96093169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.075045295) q[1];
sx q[1];
rz(-2.5270224) q[1];
sx q[1];
rz(-0.11654186) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0876168) q[3];
sx q[3];
rz(-0.91435104) q[3];
sx q[3];
rz(2.2051728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1571618) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(-2.7471527) q[2];
rz(-1.1557584) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(-1.8331147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406469) q[0];
sx q[0];
rz(-0.41579682) q[0];
sx q[0];
rz(-1.0619324) q[0];
rz(2.1901219) q[1];
sx q[1];
rz(-2.2679592) q[1];
sx q[1];
rz(-1.6208614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3979231) q[0];
sx q[0];
rz(-2.8594404) q[0];
sx q[0];
rz(-1.762853) q[0];
rz(-2.4302519) q[2];
sx q[2];
rz(-0.78186505) q[2];
sx q[2];
rz(3.1202673) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9545867) q[1];
sx q[1];
rz(-2.329064) q[1];
sx q[1];
rz(-1.345326) q[1];
x q[2];
rz(2.7476375) q[3];
sx q[3];
rz(-1.3666913) q[3];
sx q[3];
rz(0.98180727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1157397) q[2];
sx q[2];
rz(-1.9151177) q[2];
sx q[2];
rz(2.8958877) q[2];
rz(1.687441) q[3];
sx q[3];
rz(-1.629963) q[3];
sx q[3];
rz(2.3112442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65569735) q[0];
sx q[0];
rz(-0.73807722) q[0];
sx q[0];
rz(-0.6947211) q[0];
rz(-0.10063902) q[1];
sx q[1];
rz(-1.0357608) q[1];
sx q[1];
rz(-2.2763841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4072587) q[0];
sx q[0];
rz(-1.6122804) q[0];
sx q[0];
rz(-2.7232552) q[0];
x q[1];
rz(2.2312282) q[2];
sx q[2];
rz(-1.5857134) q[2];
sx q[2];
rz(2.2098324) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45564902) q[1];
sx q[1];
rz(-2.5091617) q[1];
sx q[1];
rz(1.8258105) q[1];
rz(0.67460028) q[3];
sx q[3];
rz(-2.018542) q[3];
sx q[3];
rz(-2.646776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9575017) q[2];
sx q[2];
rz(-0.75347334) q[2];
sx q[2];
rz(0.91020477) q[2];
rz(1.6893859) q[3];
sx q[3];
rz(-1.9096749) q[3];
sx q[3];
rz(0.46428251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.8978187) q[0];
sx q[0];
rz(-2.0938566) q[0];
sx q[0];
rz(0.62328231) q[0];
rz(2.9858164) q[1];
sx q[1];
rz(-2.364295) q[1];
sx q[1];
rz(-0.26330858) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1415188) q[0];
sx q[0];
rz(-0.21528582) q[0];
sx q[0];
rz(1.5207284) q[0];
rz(-pi) q[1];
rz(-1.6313577) q[2];
sx q[2];
rz(-1.4944296) q[2];
sx q[2];
rz(1.5431984) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3015107) q[1];
sx q[1];
rz(-2.5599944) q[1];
sx q[1];
rz(-0.14029293) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5978388) q[3];
sx q[3];
rz(-2.2787146) q[3];
sx q[3];
rz(1.6342362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2498563) q[2];
sx q[2];
rz(-0.24954924) q[2];
sx q[2];
rz(0.88805324) q[2];
rz(-0.92787162) q[3];
sx q[3];
rz(-1.1279305) q[3];
sx q[3];
rz(-2.1685062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71858281) q[0];
sx q[0];
rz(-0.50964481) q[0];
sx q[0];
rz(-0.50073671) q[0];
rz(-1.1205193) q[1];
sx q[1];
rz(-2.3605533) q[1];
sx q[1];
rz(2.3401071) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716176) q[0];
sx q[0];
rz(-1.5688251) q[0];
sx q[0];
rz(3.0706186) q[0];
rz(-pi) q[1];
rz(-1.0455442) q[2];
sx q[2];
rz(-2.0727951) q[2];
sx q[2];
rz(-1.0848622) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5902596) q[1];
sx q[1];
rz(-2.1864744) q[1];
sx q[1];
rz(-2.1034296) q[1];
rz(0.19480199) q[3];
sx q[3];
rz(-1.8240989) q[3];
sx q[3];
rz(-1.4089438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72006172) q[2];
sx q[2];
rz(-1.3214107) q[2];
sx q[2];
rz(-2.9113972) q[2];
rz(-3.1318393) q[3];
sx q[3];
rz(-1.6246656) q[3];
sx q[3];
rz(0.15596381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941876) q[0];
sx q[0];
rz(-1.367584) q[0];
sx q[0];
rz(0.76960027) q[0];
rz(-2.1790478) q[1];
sx q[1];
rz(-1.3178408) q[1];
sx q[1];
rz(-2.1544971) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5819823) q[0];
sx q[0];
rz(-1.6405237) q[0];
sx q[0];
rz(1.0435253) q[0];
rz(-pi) q[1];
rz(-1.8136421) q[2];
sx q[2];
rz(-1.4220024) q[2];
sx q[2];
rz(0.75243581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7915114) q[1];
sx q[1];
rz(-1.529585) q[1];
sx q[1];
rz(-3.1084395) q[1];
rz(-pi) q[2];
rz(-0.77642595) q[3];
sx q[3];
rz(-1.2033495) q[3];
sx q[3];
rz(-2.020379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4276245) q[2];
sx q[2];
rz(-2.2043113) q[2];
sx q[2];
rz(-3.0734708) q[2];
rz(-2.6918329) q[3];
sx q[3];
rz(-0.41523784) q[3];
sx q[3];
rz(1.4382039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9925256) q[0];
sx q[0];
rz(-1.5328007) q[0];
sx q[0];
rz(-1.3971064) q[0];
rz(2.8515011) q[1];
sx q[1];
rz(-2.7128704) q[1];
sx q[1];
rz(-2.6147978) q[1];
rz(-2.5505603) q[2];
sx q[2];
rz(-2.5606511) q[2];
sx q[2];
rz(-2.64369) q[2];
rz(1.5896593) q[3];
sx q[3];
rz(-1.1723779) q[3];
sx q[3];
rz(1.0642024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
