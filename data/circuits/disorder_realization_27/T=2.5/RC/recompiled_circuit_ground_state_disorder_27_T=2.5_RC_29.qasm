OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5840924) q[0];
sx q[0];
rz(3.1182365) q[0];
sx q[0];
rz(11.630848) q[0];
rz(1.525653) q[1];
sx q[1];
rz(-1.6211809) q[1];
sx q[1];
rz(0.27443767) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6870428) q[0];
sx q[0];
rz(-1.6419715) q[0];
sx q[0];
rz(-3.073259) q[0];
x q[1];
rz(2.4440632) q[2];
sx q[2];
rz(-2.4781057) q[2];
sx q[2];
rz(1.0349719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5554065) q[1];
sx q[1];
rz(-1.5342536) q[1];
sx q[1];
rz(3.1312185) q[1];
rz(-1.1660444) q[3];
sx q[3];
rz(-0.18530986) q[3];
sx q[3];
rz(2.6626793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2157669) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(2.0634148) q[2];
rz(-0.82873851) q[3];
sx q[3];
rz(-1.5107892) q[3];
sx q[3];
rz(-2.3327648) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0593798) q[0];
sx q[0];
rz(-1.8635211) q[0];
sx q[0];
rz(1.7510121) q[0];
rz(1.7104205) q[1];
sx q[1];
rz(-3.1371959) q[1];
sx q[1];
rz(1.4336525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2947049) q[0];
sx q[0];
rz(-1.4661015) q[0];
sx q[0];
rz(-0.93331915) q[0];
rz(0.016489224) q[2];
sx q[2];
rz(-1.1270583) q[2];
sx q[2];
rz(-1.6042142) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2880651) q[1];
sx q[1];
rz(-1.57987) q[1];
sx q[1];
rz(-0.00041043343) q[1];
rz(-pi) q[2];
rz(0.065304718) q[3];
sx q[3];
rz(-2.3529144) q[3];
sx q[3];
rz(-2.1826133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7626875) q[2];
sx q[2];
rz(-1.5451558) q[2];
sx q[2];
rz(-1.5691266) q[2];
rz(-0.76111859) q[3];
sx q[3];
rz(-0.053839024) q[3];
sx q[3];
rz(-1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222027) q[0];
sx q[0];
rz(-2.5313105) q[0];
sx q[0];
rz(0.56184226) q[0];
rz(1.5678844) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(-0.017339658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2527232) q[0];
sx q[0];
rz(-0.66735744) q[0];
sx q[0];
rz(0.86020893) q[0];
rz(2.1465535) q[2];
sx q[2];
rz(-2.3684566) q[2];
sx q[2];
rz(2.7553158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9063532) q[1];
sx q[1];
rz(-1.5548551) q[1];
sx q[1];
rz(1.2083183) q[1];
rz(-0.51029499) q[3];
sx q[3];
rz(-1.5507574) q[3];
sx q[3];
rz(0.35773548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6277546) q[2];
sx q[2];
rz(-1.4189812) q[2];
sx q[2];
rz(-0.55297744) q[2];
rz(-1.9364457) q[3];
sx q[3];
rz(-1.5670992) q[3];
sx q[3];
rz(-1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7503081) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(-1.3543825) q[0];
rz(1.3722108) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(-1.7599958) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171661) q[0];
sx q[0];
rz(-1.9726511) q[0];
sx q[0];
rz(-1.0665994) q[0];
rz(1.5731029) q[2];
sx q[2];
rz(-1.5735642) q[2];
sx q[2];
rz(1.6244013) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.957037) q[1];
sx q[1];
rz(-0.77691764) q[1];
sx q[1];
rz(-2.0691815) q[1];
rz(0.28691157) q[3];
sx q[3];
rz(-0.82334703) q[3];
sx q[3];
rz(-0.2940184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0682721) q[2];
sx q[2];
rz(-3.1219411) q[2];
sx q[2];
rz(1.2400631) q[2];
rz(-2.9114919) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(0.41964644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0536026) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(1.2552274) q[0];
rz(-0.0093731006) q[1];
sx q[1];
rz(-1.3690989) q[1];
sx q[1];
rz(3.110041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2173715) q[0];
sx q[0];
rz(-1.517202) q[0];
sx q[0];
rz(1.0384667) q[0];
x q[1];
rz(3.0537082) q[2];
sx q[2];
rz(-1.8586718) q[2];
sx q[2];
rz(-0.56729546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7158791) q[1];
sx q[1];
rz(-1.3495933) q[1];
sx q[1];
rz(-2.1137266) q[1];
rz(-pi) q[2];
rz(0.95057861) q[3];
sx q[3];
rz(-1.0744541) q[3];
sx q[3];
rz(-0.84211189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3343398) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(-0.80017153) q[2];
rz(-2.3470894) q[3];
sx q[3];
rz(-3.1092293) q[3];
sx q[3];
rz(-0.86875027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.054258) q[0];
sx q[0];
rz(-0.15905173) q[0];
sx q[0];
rz(1.4999088) q[0];
rz(0.17290393) q[1];
sx q[1];
rz(-3.0992295) q[1];
sx q[1];
rz(-0.049887966) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5373851) q[0];
sx q[0];
rz(-2.4695463) q[0];
sx q[0];
rz(-1.9792569) q[0];
rz(1.6880269) q[2];
sx q[2];
rz(-0.77146009) q[2];
sx q[2];
rz(-1.7457123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3059002) q[1];
sx q[1];
rz(-2.3819807) q[1];
sx q[1];
rz(1.2051969) q[1];
x q[2];
rz(0.98007794) q[3];
sx q[3];
rz(-0.81390611) q[3];
sx q[3];
rz(0.68099672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40338966) q[2];
sx q[2];
rz(-0.047813606) q[2];
sx q[2];
rz(-1.2460463) q[2];
rz(-1.7854779) q[3];
sx q[3];
rz(-3.1062283) q[3];
sx q[3];
rz(-1.7252007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7127011) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(1.7190546) q[0];
rz(1.2369583) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(-2.9388156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3966438) q[0];
sx q[0];
rz(-2.2988964) q[0];
sx q[0];
rz(-2.6999695) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70901633) q[2];
sx q[2];
rz(-1.0142027) q[2];
sx q[2];
rz(0.32278827) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8494871) q[1];
sx q[1];
rz(-1.6260481) q[1];
sx q[1];
rz(-0.34088742) q[1];
x q[2];
rz(-0.15722991) q[3];
sx q[3];
rz(-2.0766602) q[3];
sx q[3];
rz(-2.0802405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5286336) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(-2.4670777) q[2];
rz(1.7965192) q[3];
sx q[3];
rz(-0.14480545) q[3];
sx q[3];
rz(-1.8176746) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26789185) q[0];
sx q[0];
rz(-2.38509) q[0];
sx q[0];
rz(-0.85195136) q[0];
rz(2.9497228) q[1];
sx q[1];
rz(-3.1286897) q[1];
sx q[1];
rz(-0.26564863) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0582377) q[0];
sx q[0];
rz(-1.7449656) q[0];
sx q[0];
rz(0.41634286) q[0];
rz(-2.526409) q[2];
sx q[2];
rz(-1.8447723) q[2];
sx q[2];
rz(-0.09466234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77924624) q[1];
sx q[1];
rz(-1.4873449) q[1];
sx q[1];
rz(-1.6298184) q[1];
rz(-2.5784078) q[3];
sx q[3];
rz(-0.63473375) q[3];
sx q[3];
rz(-0.70574443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77136451) q[2];
sx q[2];
rz(-3.0493272) q[2];
sx q[2];
rz(-2.6293758) q[2];
rz(0.15906119) q[3];
sx q[3];
rz(-3.1064807) q[3];
sx q[3];
rz(-1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8856186) q[0];
sx q[0];
rz(-1.4775448) q[0];
sx q[0];
rz(1.0783827) q[0];
rz(-1.501561) q[1];
sx q[1];
rz(-0.20137943) q[1];
sx q[1];
rz(-1.5578425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5974332) q[0];
sx q[0];
rz(-2.2578116) q[0];
sx q[0];
rz(1.5907955) q[0];
rz(-pi) q[1];
rz(-0.89084004) q[2];
sx q[2];
rz(-1.8906697) q[2];
sx q[2];
rz(0.25242463) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.16405995) q[1];
sx q[1];
rz(-1.5674453) q[1];
sx q[1];
rz(-3.1195669) q[1];
x q[2];
rz(-2.0191865) q[3];
sx q[3];
rz(-2.4351154) q[3];
sx q[3];
rz(-2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.885159) q[2];
sx q[2];
rz(-0.014545518) q[2];
sx q[2];
rz(3.0719768) q[2];
rz(0.066667892) q[3];
sx q[3];
rz(-2.127141) q[3];
sx q[3];
rz(2.4623509) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081864) q[0];
sx q[0];
rz(-1.2766159) q[0];
sx q[0];
rz(-2.8896914) q[0];
rz(1.6690286) q[1];
sx q[1];
rz(-2.9258969) q[1];
sx q[1];
rz(3.0676945) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65739261) q[0];
sx q[0];
rz(-1.9354685) q[0];
sx q[0];
rz(-0.075693746) q[0];
rz(-pi) q[1];
rz(1.6021219) q[2];
sx q[2];
rz(-1.5534473) q[2];
sx q[2];
rz(0.67048873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.042797814) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(1.2513729) q[1];
x q[2];
rz(0.68040397) q[3];
sx q[3];
rz(-0.26016737) q[3];
sx q[3];
rz(2.5870067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68022388) q[2];
sx q[2];
rz(-0.0071439925) q[2];
sx q[2];
rz(2.3738677) q[2];
rz(-1.7420306) q[3];
sx q[3];
rz(-3.1407686) q[3];
sx q[3];
rz(-2.6373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298228) q[0];
sx q[0];
rz(-0.98291021) q[0];
sx q[0];
rz(1.7194189) q[0];
rz(3.1172251) q[1];
sx q[1];
rz(-2.9822646) q[1];
sx q[1];
rz(0.23039625) q[1];
rz(-2.2461535) q[2];
sx q[2];
rz(-1.2661305) q[2];
sx q[2];
rz(0.2998395) q[2];
rz(2.92504) q[3];
sx q[3];
rz(-0.42600507) q[3];
sx q[3];
rz(1.101839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
