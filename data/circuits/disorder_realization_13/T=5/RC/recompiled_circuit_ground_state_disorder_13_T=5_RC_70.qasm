OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.00081113022) q[0];
sx q[0];
rz(-0.82511628) q[0];
sx q[0];
rz(-2.7197279) q[0];
rz(-0.74692625) q[1];
sx q[1];
rz(-2.5513625) q[1];
sx q[1];
rz(2.3205369) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29971545) q[0];
sx q[0];
rz(-1.6546735) q[0];
sx q[0];
rz(-1.6158577) q[0];
x q[1];
rz(-1.7763605) q[2];
sx q[2];
rz(-1.7920307) q[2];
sx q[2];
rz(-1.6541755) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6118879) q[1];
sx q[1];
rz(-0.48445736) q[1];
sx q[1];
rz(2.2147708) q[1];
x q[2];
rz(-0.13387891) q[3];
sx q[3];
rz(-1.9124096) q[3];
sx q[3];
rz(0.70646369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7464298) q[2];
sx q[2];
rz(-0.91150993) q[2];
sx q[2];
rz(0.54473031) q[2];
rz(1.644545) q[3];
sx q[3];
rz(-1.112273) q[3];
sx q[3];
rz(-1.3397608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4561975) q[0];
sx q[0];
rz(-2.39769) q[0];
sx q[0];
rz(-2.5681382) q[0];
rz(0.042512976) q[1];
sx q[1];
rz(-1.7492234) q[1];
sx q[1];
rz(-1.858985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3488878) q[0];
sx q[0];
rz(-0.12785873) q[0];
sx q[0];
rz(-0.16221817) q[0];
x q[1];
rz(-0.71252026) q[2];
sx q[2];
rz(-2.7606568) q[2];
sx q[2];
rz(2.088758) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6635884) q[1];
sx q[1];
rz(-1.4299031) q[1];
sx q[1];
rz(-2.6114467) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3815879) q[3];
sx q[3];
rz(-0.78746057) q[3];
sx q[3];
rz(1.2889047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87806988) q[2];
sx q[2];
rz(-1.9922549) q[2];
sx q[2];
rz(0.44352201) q[2];
rz(1.589132) q[3];
sx q[3];
rz(-0.3905206) q[3];
sx q[3];
rz(-0.21903567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1294915) q[0];
sx q[0];
rz(-2.2624367) q[0];
sx q[0];
rz(0.40170676) q[0];
rz(-1.7237639) q[1];
sx q[1];
rz(-2.6938853) q[1];
sx q[1];
rz(-2.0254501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48211574) q[0];
sx q[0];
rz(-1.8748613) q[0];
sx q[0];
rz(0.67993865) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7557451) q[2];
sx q[2];
rz(-2.4413337) q[2];
sx q[2];
rz(-1.2724233) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9520719) q[1];
sx q[1];
rz(-2.4671989) q[1];
sx q[1];
rz(-2.2728069) q[1];
rz(-pi) q[2];
rz(-0.54858923) q[3];
sx q[3];
rz(-1.635713) q[3];
sx q[3];
rz(1.3260076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22905722) q[2];
sx q[2];
rz(-2.8450862) q[2];
sx q[2];
rz(2.1336446) q[2];
rz(1.5063162) q[3];
sx q[3];
rz(-1.1791891) q[3];
sx q[3];
rz(-0.87872163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85663831) q[0];
sx q[0];
rz(-3.1143739) q[0];
sx q[0];
rz(-0.83754367) q[0];
rz(-1.2614999) q[1];
sx q[1];
rz(-1.7439758) q[1];
sx q[1];
rz(1.2417485) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1364131) q[0];
sx q[0];
rz(-1.2997928) q[0];
sx q[0];
rz(2.6684929) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9292232) q[2];
sx q[2];
rz(-2.1687733) q[2];
sx q[2];
rz(2.6670307) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0251522) q[1];
sx q[1];
rz(-0.63092996) q[1];
sx q[1];
rz(-2.4063944) q[1];
x q[2];
rz(2.1513274) q[3];
sx q[3];
rz(-1.255569) q[3];
sx q[3];
rz(-2.6516556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1602829) q[2];
sx q[2];
rz(-0.95429388) q[2];
sx q[2];
rz(-1.0120288) q[2];
rz(-2.1435598) q[3];
sx q[3];
rz(-2.4055552) q[3];
sx q[3];
rz(-1.6691104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0618458) q[0];
sx q[0];
rz(-2.340402) q[0];
sx q[0];
rz(-0.27859846) q[0];
rz(-0.31696907) q[1];
sx q[1];
rz(-1.7520889) q[1];
sx q[1];
rz(2.680079) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46733785) q[0];
sx q[0];
rz(-1.5227571) q[0];
sx q[0];
rz(-1.7397299) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95038173) q[2];
sx q[2];
rz(-1.912552) q[2];
sx q[2];
rz(-0.9859964) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.52109226) q[1];
sx q[1];
rz(-0.25687718) q[1];
sx q[1];
rz(1.1939373) q[1];
x q[2];
rz(1.4123437) q[3];
sx q[3];
rz(-0.55084832) q[3];
sx q[3];
rz(2.2864263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7757814) q[2];
sx q[2];
rz(-1.8041939) q[2];
sx q[2];
rz(-0.27654761) q[2];
rz(-0.60025275) q[3];
sx q[3];
rz(-0.7163896) q[3];
sx q[3];
rz(-2.6071809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.9944331) q[0];
sx q[0];
rz(-2.5614547) q[0];
sx q[0];
rz(-2.4603727) q[0];
rz(0.45411202) q[1];
sx q[1];
rz(-1.157016) q[1];
sx q[1];
rz(-0.62201321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91151448) q[0];
sx q[0];
rz(-2.921836) q[0];
sx q[0];
rz(0.4498555) q[0];
rz(-pi) q[1];
rz(-2.5037982) q[2];
sx q[2];
rz(-2.0908815) q[2];
sx q[2];
rz(-0.74842036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6890735) q[1];
sx q[1];
rz(-1.9770685) q[1];
sx q[1];
rz(2.3157332) q[1];
rz(-1.9196904) q[3];
sx q[3];
rz(-1.318228) q[3];
sx q[3];
rz(3.0608197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6260208) q[2];
sx q[2];
rz(-0.61352789) q[2];
sx q[2];
rz(0.75508368) q[2];
rz(-0.10609047) q[3];
sx q[3];
rz(-1.875149) q[3];
sx q[3];
rz(2.6809926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945011) q[0];
sx q[0];
rz(-2.5300808) q[0];
sx q[0];
rz(0.066135429) q[0];
rz(3.0905837) q[1];
sx q[1];
rz(-0.74791932) q[1];
sx q[1];
rz(1.8775108) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1544113) q[0];
sx q[0];
rz(-2.4307152) q[0];
sx q[0];
rz(2.1412503) q[0];
x q[1];
rz(2.7876652) q[2];
sx q[2];
rz(-2.3945237) q[2];
sx q[2];
rz(2.9198271) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6173827) q[1];
sx q[1];
rz(-2.2508989) q[1];
sx q[1];
rz(-1.0587803) q[1];
rz(1.1751647) q[3];
sx q[3];
rz(-1.1671178) q[3];
sx q[3];
rz(2.5772571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2833726) q[2];
sx q[2];
rz(-2.0141352) q[2];
sx q[2];
rz(-0.14632012) q[2];
rz(-0.60394168) q[3];
sx q[3];
rz(-0.54655176) q[3];
sx q[3];
rz(1.4124136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42501763) q[0];
sx q[0];
rz(-1.1973493) q[0];
sx q[0];
rz(-3.0498411) q[0];
rz(-0.07235202) q[1];
sx q[1];
rz(-1.8232583) q[1];
sx q[1];
rz(2.4718557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7726752) q[0];
sx q[0];
rz(-1.8631422) q[0];
sx q[0];
rz(1.4293074) q[0];
x q[1];
rz(-0.40939949) q[2];
sx q[2];
rz(-1.690662) q[2];
sx q[2];
rz(1.9425962) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.32423985) q[1];
sx q[1];
rz(-1.8756556) q[1];
sx q[1];
rz(1.7225811) q[1];
x q[2];
rz(2.9961139) q[3];
sx q[3];
rz(-2.0091741) q[3];
sx q[3];
rz(-2.9117025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3503795) q[2];
sx q[2];
rz(-0.71724856) q[2];
sx q[2];
rz(-2.4199602) q[2];
rz(2.3676938) q[3];
sx q[3];
rz(-2.1589203) q[3];
sx q[3];
rz(0.53988808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72188193) q[0];
sx q[0];
rz(-1.1996491) q[0];
sx q[0];
rz(-2.6019959) q[0];
rz(-2.3609912) q[1];
sx q[1];
rz(-1.8949948) q[1];
sx q[1];
rz(-2.7194729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1273859) q[0];
sx q[0];
rz(-0.84833586) q[0];
sx q[0];
rz(-2.3425472) q[0];
rz(-pi) q[1];
rz(-1.6306179) q[2];
sx q[2];
rz(-2.0622396) q[2];
sx q[2];
rz(0.37183842) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6680702) q[1];
sx q[1];
rz(-2.0501185) q[1];
sx q[1];
rz(3.098686) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5255181) q[3];
sx q[3];
rz(-2.1589734) q[3];
sx q[3];
rz(-0.28376337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0427986) q[2];
sx q[2];
rz(-1.021794) q[2];
sx q[2];
rz(-0.54879028) q[2];
rz(2.9889066) q[3];
sx q[3];
rz(-1.0278253) q[3];
sx q[3];
rz(-2.4299183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76950908) q[0];
sx q[0];
rz(-2.7476269) q[0];
sx q[0];
rz(-2.5373996) q[0];
rz(1.6744042) q[1];
sx q[1];
rz(-1.3507651) q[1];
sx q[1];
rz(-1.1473354) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98367386) q[0];
sx q[0];
rz(-2.3566453) q[0];
sx q[0];
rz(-0.60529937) q[0];
x q[1];
rz(-1.2560955) q[2];
sx q[2];
rz(-2.2813548) q[2];
sx q[2];
rz(1.1473555) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5903476) q[1];
sx q[1];
rz(-2.6585124) q[1];
sx q[1];
rz(0.03309588) q[1];
rz(1.9705521) q[3];
sx q[3];
rz(-1.2467524) q[3];
sx q[3];
rz(-0.43982616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6818105) q[2];
sx q[2];
rz(-1.9777538) q[2];
sx q[2];
rz(0.60010827) q[2];
rz(2.4188304) q[3];
sx q[3];
rz(-0.99812752) q[3];
sx q[3];
rz(-0.37989894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.4031274) q[0];
sx q[0];
rz(-1.5140139) q[0];
sx q[0];
rz(-1.4366666) q[0];
rz(1.5557095) q[1];
sx q[1];
rz(-1.7441505) q[1];
sx q[1];
rz(2.2081262) q[1];
rz(2.6513097) q[2];
sx q[2];
rz(-0.96302196) q[2];
sx q[2];
rz(-1.8204126) q[2];
rz(-0.24642701) q[3];
sx q[3];
rz(-1.3613614) q[3];
sx q[3];
rz(-2.341933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
