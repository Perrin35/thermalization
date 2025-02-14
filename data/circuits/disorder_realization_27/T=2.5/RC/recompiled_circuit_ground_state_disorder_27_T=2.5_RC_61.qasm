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
rz(1.525653) q[1];
sx q[1];
rz(-1.6211809) q[1];
sx q[1];
rz(0.27443767) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12111353) q[0];
sx q[0];
rz(-1.6389567) q[0];
sx q[0];
rz(-1.6421374) q[0];
rz(-pi) q[1];
rz(0.53977592) q[2];
sx q[2];
rz(-1.9775043) q[2];
sx q[2];
rz(1.1193898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2787267) q[1];
sx q[1];
rz(-3.1036065) q[1];
sx q[1];
rz(-1.8472865) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9755483) q[3];
sx q[3];
rz(-0.18530986) q[3];
sx q[3];
rz(0.47891339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9258257) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(-2.0634148) q[2];
rz(-2.3128541) q[3];
sx q[3];
rz(-1.5107892) q[3];
sx q[3];
rz(2.3327648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0593798) q[0];
sx q[0];
rz(-1.2780715) q[0];
sx q[0];
rz(1.3905806) q[0];
rz(-1.4311721) q[1];
sx q[1];
rz(-3.1371959) q[1];
sx q[1];
rz(1.4336525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2947049) q[0];
sx q[0];
rz(-1.4661015) q[0];
sx q[0];
rz(-2.2082735) q[0];
rz(-pi) q[1];
rz(2.0145871) q[2];
sx q[2];
rz(-1.5559042) q[2];
sx q[2];
rz(0.040497517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8588577) q[1];
sx q[1];
rz(-1.5712067) q[1];
sx q[1];
rz(1.5617227) q[1];
rz(-pi) q[2];
rz(3.0762879) q[3];
sx q[3];
rz(-0.78867824) q[3];
sx q[3];
rz(0.95897934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7626875) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(-1.5691266) q[2];
rz(2.3804741) q[3];
sx q[3];
rz(-3.0877536) q[3];
sx q[3];
rz(1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91956562) q[0];
sx q[0];
rz(-0.61028218) q[0];
sx q[0];
rz(-2.5797504) q[0];
rz(-1.5678844) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(-3.124253) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8888694) q[0];
sx q[0];
rz(-0.66735744) q[0];
sx q[0];
rz(2.2813837) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2566913) q[2];
sx q[2];
rz(-1.1807311) q[2];
sx q[2];
rz(1.6194413) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3295118) q[1];
sx q[1];
rz(-1.9332262) q[1];
sx q[1];
rz(3.1245438) q[1];
x q[2];
rz(1.5937599) q[3];
sx q[3];
rz(-2.0809789) q[3];
sx q[3];
rz(1.2018454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6277546) q[2];
sx q[2];
rz(-1.7226115) q[2];
sx q[2];
rz(0.55297744) q[2];
rz(-1.2051469) q[3];
sx q[3];
rz(-1.5670992) q[3];
sx q[3];
rz(-1.5470101) q[3];
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
rz(0.39128458) q[0];
sx q[0];
rz(-1.0306083) q[0];
sx q[0];
rz(1.7872101) q[0];
rz(1.7693819) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(1.7599958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87028044) q[0];
sx q[0];
rz(-0.63374937) q[0];
sx q[0];
rz(2.2922754) q[0];
x q[1];
rz(0.69475485) q[2];
sx q[2];
rz(-3.1379897) q[2];
sx q[2];
rz(0.82243332) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83643944) q[1];
sx q[1];
rz(-0.90739319) q[1];
sx q[1];
rz(-0.43933489) q[1];
x q[2];
rz(-2.8546811) q[3];
sx q[3];
rz(-0.82334703) q[3];
sx q[3];
rz(-0.2940184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0682721) q[2];
sx q[2];
rz(-3.1219411) q[2];
sx q[2];
rz(-1.9015296) q[2];
rz(0.23010075) q[3];
sx q[3];
rz(-3.1374044) q[3];
sx q[3];
rz(-0.41964644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087990046) q[0];
sx q[0];
rz(-2.3542861) q[0];
sx q[0];
rz(-1.8863652) q[0];
rz(-0.0093731006) q[1];
sx q[1];
rz(-1.3690989) q[1];
sx q[1];
rz(-0.031551687) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7566236) q[0];
sx q[0];
rz(-2.1022804) q[0];
sx q[0];
rz(3.079412) q[0];
rz(-1.2826233) q[2];
sx q[2];
rz(-0.30063486) q[2];
sx q[2];
rz(-2.2733781) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
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
rz(-2.8072529) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(2.3414211) q[2];
rz(2.3470894) q[3];
sx q[3];
rz(-0.032363351) q[3];
sx q[3];
rz(-0.86875027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0873347) q[0];
sx q[0];
rz(-2.9825409) q[0];
sx q[0];
rz(1.4999088) q[0];
rz(-2.9686887) q[1];
sx q[1];
rz(-3.0992295) q[1];
sx q[1];
rz(-0.049887966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2931516) q[0];
sx q[0];
rz(-1.8206791) q[0];
sx q[0];
rz(0.94012733) q[0];
x q[1];
rz(-0.80277819) q[2];
sx q[2];
rz(-1.4891616) q[2];
sx q[2];
rz(2.8824474) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3499455) q[1];
sx q[1];
rz(-2.269372) q[1];
sx q[1];
rz(-0.32731815) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1615147) q[3];
sx q[3];
rz(-0.81390611) q[3];
sx q[3];
rz(0.68099672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4288915) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(1.7190546) q[0];
rz(1.9046344) q[1];
sx q[1];
rz(-3.1045541) q[1];
sx q[1];
rz(-2.9388156) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0109459) q[0];
sx q[0];
rz(-1.8955064) q[0];
sx q[0];
rz(-2.3492011) q[0];
x q[1];
rz(-0.76272623) q[2];
sx q[2];
rz(-2.2710851) q[2];
sx q[2];
rz(-1.3415847) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2591035) q[1];
sx q[1];
rz(-1.2304502) q[1];
sx q[1];
rz(-1.5121786) q[1];
rz(-pi) q[2];
rz(1.0596541) q[3];
sx q[3];
rz(-1.7082001) q[3];
sx q[3];
rz(0.58611548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.612959) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(-0.67451492) q[2];
rz(1.3450735) q[3];
sx q[3];
rz(-0.14480545) q[3];
sx q[3];
rz(1.8176746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26789185) q[0];
sx q[0];
rz(-2.38509) q[0];
sx q[0];
rz(-2.2896413) q[0];
rz(0.1918699) q[1];
sx q[1];
rz(-3.1286897) q[1];
sx q[1];
rz(0.26564863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083355) q[0];
sx q[0];
rz(-1.3966271) q[0];
sx q[0];
rz(2.7252498) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9022361) q[2];
sx q[2];
rz(-2.1598615) q[2];
sx q[2];
rz(1.6650496) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3959256) q[1];
sx q[1];
rz(-3.0394181) q[1];
sx q[1];
rz(2.5273782) q[1];
rz(2.5846768) q[3];
sx q[3];
rz(-1.2486826) q[3];
sx q[3];
rz(-0.39469257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77136451) q[2];
sx q[2];
rz(-0.092265487) q[2];
sx q[2];
rz(-2.6293758) q[2];
rz(2.9825315) q[3];
sx q[3];
rz(-3.1064807) q[3];
sx q[3];
rz(1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8856186) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(-2.0632099) q[0];
rz(-1.6400317) q[1];
sx q[1];
rz(-2.9402132) q[1];
sx q[1];
rz(-1.5578425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1276409) q[0];
sx q[0];
rz(-1.5553345) q[0];
sx q[0];
rz(0.68711335) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0859162) q[2];
sx q[2];
rz(-0.74046248) q[2];
sx q[2];
rz(-1.452335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7347825) q[1];
sx q[1];
rz(-1.5487707) q[1];
sx q[1];
rz(1.5674445) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35435481) q[3];
sx q[3];
rz(-2.1956596) q[3];
sx q[3];
rz(-0.55312485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25643361) q[2];
sx q[2];
rz(-0.014545518) q[2];
sx q[2];
rz(3.0719768) q[2];
rz(3.0749248) q[3];
sx q[3];
rz(-2.127141) q[3];
sx q[3];
rz(0.67924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081864) q[0];
sx q[0];
rz(-1.8649768) q[0];
sx q[0];
rz(2.8896914) q[0];
rz(-1.6690286) q[1];
sx q[1];
rz(-0.21569574) q[1];
sx q[1];
rz(-0.073898166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65739261) q[0];
sx q[0];
rz(-1.2061242) q[0];
sx q[0];
rz(3.0658989) q[0];
x q[1];
rz(0.017357512) q[2];
sx q[2];
rz(-1.6021172) q[2];
sx q[2];
rz(2.2407414) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.042797814) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(1.8902197) q[1];
rz(-1.7367253) q[3];
sx q[3];
rz(-1.36948) q[3];
sx q[3];
rz(1.8898058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68022388) q[2];
sx q[2];
rz(-0.0071439925) q[2];
sx q[2];
rz(0.76772493) q[2];
rz(-1.7420306) q[3];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
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
rz(0.89543912) q[2];
sx q[2];
rz(-1.2661305) q[2];
sx q[2];
rz(0.2998395) q[2];
rz(0.21655269) q[3];
sx q[3];
rz(-2.7155876) q[3];
sx q[3];
rz(-2.0397536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
