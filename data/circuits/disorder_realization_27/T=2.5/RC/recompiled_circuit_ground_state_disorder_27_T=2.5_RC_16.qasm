OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55750027) q[0];
sx q[0];
rz(-3.1182365) q[0];
sx q[0];
rz(-0.93552247) q[0];
rz(-1.6159396) q[1];
sx q[1];
rz(-1.5204117) q[1];
sx q[1];
rz(-0.27443767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0204791) q[0];
sx q[0];
rz(-1.6389567) q[0];
sx q[0];
rz(1.6421374) q[0];
rz(-pi) q[1];
rz(-2.0361314) q[2];
sx q[2];
rz(-1.0792152) q[2];
sx q[2];
rz(0.21869379) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9842311) q[1];
sx q[1];
rz(-1.5604291) q[1];
sx q[1];
rz(-1.5342516) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7414344) q[3];
sx q[3];
rz(-1.4981761) q[3];
sx q[3];
rz(-0.69334465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9258257) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(2.0634148) q[2];
rz(-2.3128541) q[3];
sx q[3];
rz(-1.5107892) q[3];
sx q[3];
rz(-0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0593798) q[0];
sx q[0];
rz(-1.8635211) q[0];
sx q[0];
rz(1.7510121) q[0];
rz(1.4311721) q[1];
sx q[1];
rz(-3.1371959) q[1];
sx q[1];
rz(1.7079401) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9427458) q[0];
sx q[0];
rz(-2.2042243) q[0];
sx q[0];
rz(-0.1300227) q[0];
rz(-pi) q[1];
rz(1.1270056) q[2];
sx q[2];
rz(-1.5559042) q[2];
sx q[2];
rz(-0.040497517) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8535276) q[1];
sx q[1];
rz(-1.5617227) q[1];
sx q[1];
rz(-0.00041043343) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0762879) q[3];
sx q[3];
rz(-0.78867824) q[3];
sx q[3];
rz(-2.1826133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7626875) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(1.572466) q[2];
rz(-0.76111859) q[3];
sx q[3];
rz(-3.0877536) q[3];
sx q[3];
rz(-1.2083453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222027) q[0];
sx q[0];
rz(-0.61028218) q[0];
sx q[0];
rz(-0.56184226) q[0];
rz(1.5678844) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(3.124253) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8888694) q[0];
sx q[0];
rz(-0.66735744) q[0];
sx q[0];
rz(-0.86020893) q[0];
rz(-pi) q[1];
rz(0.88490136) q[2];
sx q[2];
rz(-1.1807311) q[2];
sx q[2];
rz(-1.5221514) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9063532) q[1];
sx q[1];
rz(-1.5548551) q[1];
sx q[1];
rz(-1.9332744) q[1];
rz(-pi) q[2];
rz(-1.5937599) q[3];
sx q[3];
rz(-2.0809789) q[3];
sx q[3];
rz(1.9397473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6277546) q[2];
sx q[2];
rz(-1.4189812) q[2];
sx q[2];
rz(2.5886152) q[2];
rz(-1.9364457) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(-1.5470101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7503081) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(1.7872101) q[0];
rz(-1.7693819) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(-1.7599958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8244265) q[0];
sx q[0];
rz(-1.9726511) q[0];
sx q[0];
rz(2.0749932) q[0];
x q[1];
rz(2.4468378) q[2];
sx q[2];
rz(-3.1379897) q[2];
sx q[2];
rz(2.3191593) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18455566) q[1];
sx q[1];
rz(-2.364675) q[1];
sx q[1];
rz(-2.0691815) q[1];
rz(-pi) q[2];
rz(-0.28691157) q[3];
sx q[3];
rz(-2.3182456) q[3];
sx q[3];
rz(-0.2940184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0733205) q[2];
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
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0536026) q[0];
sx q[0];
rz(-2.3542861) q[0];
sx q[0];
rz(-1.2552274) q[0];
rz(-0.0093731006) q[1];
sx q[1];
rz(-1.3690989) q[1];
sx q[1];
rz(-0.031551687) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9242212) q[0];
sx q[0];
rz(-1.517202) q[0];
sx q[0];
rz(-1.0384667) q[0];
x q[1];
rz(-1.8597262) q[2];
sx q[2];
rz(-1.4865371) q[2];
sx q[2];
rz(-2.1631028) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1281368) q[1];
sx q[1];
rz(-1.0425048) q[1];
sx q[1];
rz(-2.8847413) q[1];
rz(-pi) q[2];
rz(0.58720354) q[3];
sx q[3];
rz(-2.1072344) q[3];
sx q[3];
rz(2.0850542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8072529) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(0.80017153) q[2];
rz(2.3470894) q[3];
sx q[3];
rz(-3.1092293) q[3];
sx q[3];
rz(-2.2728424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0873347) q[0];
sx q[0];
rz(-0.15905173) q[0];
sx q[0];
rz(1.6416838) q[0];
rz(-2.9686887) q[1];
sx q[1];
rz(-0.042363107) q[1];
sx q[1];
rz(-3.0917047) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.848441) q[0];
sx q[0];
rz(-1.8206791) q[0];
sx q[0];
rz(0.94012733) q[0];
x q[1];
rz(1.6880269) q[2];
sx q[2];
rz(-0.77146009) q[2];
sx q[2];
rz(-1.7457123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8356925) q[1];
sx q[1];
rz(-2.3819807) q[1];
sx q[1];
rz(-1.9363957) q[1];
rz(-pi) q[2];
rz(-0.53277758) q[3];
sx q[3];
rz(-0.92255892) q[3];
sx q[3];
rz(1.6870354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40338966) q[2];
sx q[2];
rz(-0.047813606) q[2];
sx q[2];
rz(-1.2460463) q[2];
rz(-1.3561148) q[3];
sx q[3];
rz(-0.035364371) q[3];
sx q[3];
rz(-1.7252007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7127011) q[0];
sx q[0];
rz(-0.8242979) q[0];
sx q[0];
rz(1.4225381) q[0];
rz(-1.9046344) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(-2.9388156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0143429) q[0];
sx q[0];
rz(-0.83006751) q[0];
sx q[0];
rz(1.1237445) q[0];
x q[1];
rz(-2.3788664) q[2];
sx q[2];
rz(-2.2710851) q[2];
sx q[2];
rz(1.3415847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8494871) q[1];
sx q[1];
rz(-1.6260481) q[1];
sx q[1];
rz(-0.34088742) q[1];
rz(-pi) q[2];
rz(2.0819386) q[3];
sx q[3];
rz(-1.7082001) q[3];
sx q[3];
rz(-0.58611548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5286336) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(-0.67451492) q[2];
rz(-1.3450735) q[3];
sx q[3];
rz(-2.9967872) q[3];
sx q[3];
rz(1.8176746) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26789185) q[0];
sx q[0];
rz(-2.38509) q[0];
sx q[0];
rz(-2.2896413) q[0];
rz(-0.1918699) q[1];
sx q[1];
rz(-0.012902915) q[1];
sx q[1];
rz(-2.875944) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56391956) q[0];
sx q[0];
rz(-1.1611337) q[0];
sx q[0];
rz(1.760861) q[0];
rz(-pi) q[1];
rz(-2.688411) q[2];
sx q[2];
rz(-0.66614775) q[2];
sx q[2];
rz(-1.11048) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3959256) q[1];
sx q[1];
rz(-3.0394181) q[1];
sx q[1];
rz(2.5273782) q[1];
rz(-2.5784078) q[3];
sx q[3];
rz(-2.5068589) q[3];
sx q[3];
rz(0.70574443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3702281) q[2];
sx q[2];
rz(-3.0493272) q[2];
sx q[2];
rz(-2.6293758) q[2];
rz(2.9825315) q[3];
sx q[3];
rz(-0.03511196) q[3];
sx q[3];
rz(-1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2559741) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(1.0783827) q[0];
rz(1.501561) q[1];
sx q[1];
rz(-2.9402132) q[1];
sx q[1];
rz(1.5837502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1276409) q[0];
sx q[0];
rz(-1.5862582) q[0];
sx q[0];
rz(-0.68711335) q[0];
rz(-pi) q[1];
rz(2.2507526) q[2];
sx q[2];
rz(-1.250923) q[2];
sx q[2];
rz(2.889168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7347825) q[1];
sx q[1];
rz(-1.592822) q[1];
sx q[1];
rz(1.5674445) q[1];
x q[2];
rz(-2.0191865) q[3];
sx q[3];
rz(-0.70647722) q[3];
sx q[3];
rz(-1.1170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.885159) q[2];
sx q[2];
rz(-0.014545518) q[2];
sx q[2];
rz(-0.069615901) q[2];
rz(-0.066667892) q[3];
sx q[3];
rz(-2.127141) q[3];
sx q[3];
rz(-2.4623509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2081864) q[0];
sx q[0];
rz(-1.2766159) q[0];
sx q[0];
rz(-0.25190121) q[0];
rz(-1.6690286) q[1];
sx q[1];
rz(-0.21569574) q[1];
sx q[1];
rz(3.0676945) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4842) q[0];
sx q[0];
rz(-1.9354685) q[0];
sx q[0];
rz(-3.0658989) q[0];
rz(-pi) q[1];
rz(-1.0649071) q[2];
sx q[2];
rz(-0.035807583) q[2];
sx q[2];
rz(-0.39469013) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0987948) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(1.2513729) q[1];
rz(1.7367253) q[3];
sx q[3];
rz(-1.7721126) q[3];
sx q[3];
rz(-1.2517868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68022388) q[2];
sx q[2];
rz(-3.1344487) q[2];
sx q[2];
rz(2.3738677) q[2];
rz(1.399562) q[3];
sx q[3];
rz(-3.1407686) q[3];
sx q[3];
rz(-2.6373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298228) q[0];
sx q[0];
rz(-0.98291021) q[0];
sx q[0];
rz(1.7194189) q[0];
rz(0.024367532) q[1];
sx q[1];
rz(-0.15932803) q[1];
sx q[1];
rz(-2.9111964) q[1];
rz(-2.2461535) q[2];
sx q[2];
rz(-1.2661305) q[2];
sx q[2];
rz(0.2998395) q[2];
rz(0.41718116) q[3];
sx q[3];
rz(-1.659703) q[3];
sx q[3];
rz(-0.27123527) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
