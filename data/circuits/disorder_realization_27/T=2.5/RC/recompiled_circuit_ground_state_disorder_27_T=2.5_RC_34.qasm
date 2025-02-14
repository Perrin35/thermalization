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
rz(-0.023356181) q[0];
sx q[0];
rz(-2.2060702) q[0];
rz(-1.6159396) q[1];
sx q[1];
rz(-1.5204117) q[1];
sx q[1];
rz(2.867155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4545499) q[0];
sx q[0];
rz(-1.4996212) q[0];
sx q[0];
rz(-0.068333621) q[0];
rz(2.0361314) q[2];
sx q[2];
rz(-1.0792152) q[2];
sx q[2];
rz(-0.21869379) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1573616) q[1];
sx q[1];
rz(-1.5811635) q[1];
sx q[1];
rz(-1.5342516) q[1];
rz(1.4001582) q[3];
sx q[3];
rz(-1.4981761) q[3];
sx q[3];
rz(2.448248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2157669) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(2.0634148) q[2];
rz(0.82873851) q[3];
sx q[3];
rz(-1.6308035) q[3];
sx q[3];
rz(0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082212903) q[0];
sx q[0];
rz(-1.2780715) q[0];
sx q[0];
rz(-1.3905806) q[0];
rz(1.4311721) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(-1.7079401) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4162671) q[0];
sx q[0];
rz(-0.64483374) q[0];
sx q[0];
rz(-1.3960442) q[0];
rz(1.6054691) q[2];
sx q[2];
rz(-0.44402396) q[2];
sx q[2];
rz(-1.498986) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8588577) q[1];
sx q[1];
rz(-1.5712067) q[1];
sx q[1];
rz(1.5617227) q[1];
rz(-pi) q[2];
rz(-3.0762879) q[3];
sx q[3];
rz(-2.3529144) q[3];
sx q[3];
rz(0.95897934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37890515) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(1.572466) q[2];
rz(-0.76111859) q[3];
sx q[3];
rz(-0.053839024) q[3];
sx q[3];
rz(-1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222027) q[0];
sx q[0];
rz(-2.5313105) q[0];
sx q[0];
rz(-2.5797504) q[0];
rz(1.5678844) q[1];
sx q[1];
rz(-1.5627197) q[1];
sx q[1];
rz(-3.124253) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0837096) q[0];
sx q[0];
rz(-2.0590933) q[0];
sx q[0];
rz(-2.6668307) q[0];
x q[1];
rz(0.99503912) q[2];
sx q[2];
rz(-0.77313609) q[2];
sx q[2];
rz(2.7553158) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3775656) q[1];
sx q[1];
rz(-2.7787797) q[1];
sx q[1];
rz(1.6157263) q[1];
rz(-pi) q[2];
rz(1.5937599) q[3];
sx q[3];
rz(-1.0606137) q[3];
sx q[3];
rz(1.9397473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6277546) q[2];
sx q[2];
rz(-1.4189812) q[2];
sx q[2];
rz(-0.55297744) q[2];
rz(1.2051469) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7503081) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(-1.7872101) q[0];
rz(1.7693819) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(-1.3815968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171661) q[0];
sx q[0];
rz(-1.9726511) q[0];
sx q[0];
rz(-1.0665994) q[0];
rz(-pi) q[1];
rz(1.5684897) q[2];
sx q[2];
rz(-1.5680285) q[2];
sx q[2];
rz(1.6244013) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18455566) q[1];
sx q[1];
rz(-2.364675) q[1];
sx q[1];
rz(-2.0691815) q[1];
rz(-1.2744585) q[3];
sx q[3];
rz(-0.79056406) q[3];
sx q[3];
rz(-2.4380655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0682721) q[2];
sx q[2];
rz(-0.019651532) q[2];
sx q[2];
rz(1.2400631) q[2];
rz(-0.23010075) q[3];
sx q[3];
rz(-3.1374044) q[3];
sx q[3];
rz(-2.7219462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087990046) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(1.8863652) q[0];
rz(3.1322196) q[1];
sx q[1];
rz(-1.7724937) q[1];
sx q[1];
rz(0.031551687) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9242212) q[0];
sx q[0];
rz(-1.6243906) q[0];
sx q[0];
rz(-1.0384667) q[0];
rz(1.8597262) q[2];
sx q[2];
rz(-1.6550555) q[2];
sx q[2];
rz(0.97848985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42571354) q[1];
sx q[1];
rz(-1.3495933) q[1];
sx q[1];
rz(1.0278661) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95057861) q[3];
sx q[3];
rz(-1.0744541) q[3];
sx q[3];
rz(-0.84211189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3343398) q[2];
sx q[2];
rz(-0.0060609239) q[2];
sx q[2];
rz(2.3414211) q[2];
rz(-2.3470894) q[3];
sx q[3];
rz(-3.1092293) q[3];
sx q[3];
rz(2.2728424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0873347) q[0];
sx q[0];
rz(-2.9825409) q[0];
sx q[0];
rz(-1.6416838) q[0];
rz(-2.9686887) q[1];
sx q[1];
rz(-3.0992295) q[1];
sx q[1];
rz(-0.049887966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5373851) q[0];
sx q[0];
rz(-2.4695463) q[0];
sx q[0];
rz(1.1623357) q[0];
x q[1];
rz(2.3388145) q[2];
sx q[2];
rz(-1.652431) q[2];
sx q[2];
rz(-2.8824474) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3059002) q[1];
sx q[1];
rz(-2.3819807) q[1];
sx q[1];
rz(1.2051969) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2920554) q[3];
sx q[3];
rz(-1.153933) q[3];
sx q[3];
rz(-0.45826926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40338966) q[2];
sx q[2];
rz(-0.047813606) q[2];
sx q[2];
rz(-1.8955463) q[2];
rz(1.3561148) q[3];
sx q[3];
rz(-0.035364371) q[3];
sx q[3];
rz(1.7252007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.4288915) q[0];
sx q[0];
rz(-0.8242979) q[0];
sx q[0];
rz(1.7190546) q[0];
rz(1.9046344) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(2.9388156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3966438) q[0];
sx q[0];
rz(-2.2988964) q[0];
sx q[0];
rz(-0.44162314) q[0];
rz(-pi) q[1];
rz(-2.3788664) q[2];
sx q[2];
rz(-0.87050754) q[2];
sx q[2];
rz(-1.3415847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8824892) q[1];
sx q[1];
rz(-1.9111425) q[1];
sx q[1];
rz(-1.629414) q[1];
x q[2];
rz(-0.15722991) q[3];
sx q[3];
rz(-2.0766602) q[3];
sx q[3];
rz(-2.0802405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.612959) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(2.4670777) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8737008) q[0];
sx q[0];
rz(-0.75650263) q[0];
sx q[0];
rz(-0.85195136) q[0];
rz(2.9497228) q[1];
sx q[1];
rz(-0.012902915) q[1];
sx q[1];
rz(-2.875944) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0582377) q[0];
sx q[0];
rz(-1.7449656) q[0];
sx q[0];
rz(2.7252498) q[0];
rz(1.9022361) q[2];
sx q[2];
rz(-0.98173117) q[2];
sx q[2];
rz(1.6650496) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3959256) q[1];
sx q[1];
rz(-0.10217459) q[1];
sx q[1];
rz(2.5273782) q[1];
x q[2];
rz(1.9453796) q[3];
sx q[3];
rz(-2.0959955) q[3];
sx q[3];
rz(-1.370726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77136451) q[2];
sx q[2];
rz(-3.0493272) q[2];
sx q[2];
rz(2.6293758) q[2];
rz(0.15906119) q[3];
sx q[3];
rz(-0.03511196) q[3];
sx q[3];
rz(1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8856186) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(2.0632099) q[0];
rz(-1.501561) q[1];
sx q[1];
rz(-0.20137943) q[1];
sx q[1];
rz(1.5837502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5441594) q[0];
sx q[0];
rz(-0.88378105) q[0];
sx q[0];
rz(-1.5507971) q[0];
rz(2.0556765) q[2];
sx q[2];
rz(-2.4011302) q[2];
sx q[2];
rz(-1.6892576) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9775327) q[1];
sx q[1];
rz(-1.5674453) q[1];
sx q[1];
rz(0.022025755) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35435481) q[3];
sx q[3];
rz(-2.1956596) q[3];
sx q[3];
rz(-2.5884678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.885159) q[2];
sx q[2];
rz(-0.014545518) q[2];
sx q[2];
rz(-0.069615901) q[2];
rz(-3.0749248) q[3];
sx q[3];
rz(-2.127141) q[3];
sx q[3];
rz(-0.67924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93340623) q[0];
sx q[0];
rz(-1.8649768) q[0];
sx q[0];
rz(-2.8896914) q[0];
rz(1.6690286) q[1];
sx q[1];
rz(-2.9258969) q[1];
sx q[1];
rz(3.0676945) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4842) q[0];
sx q[0];
rz(-1.9354685) q[0];
sx q[0];
rz(-0.075693746) q[0];
rz(1.0649071) q[2];
sx q[2];
rz(-0.035807583) q[2];
sx q[2];
rz(0.39469013) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3162176) q[1];
sx q[1];
rz(-1.3298522) q[1];
sx q[1];
rz(-0.73338738) q[1];
rz(-pi) q[2];
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
rz(-2.4613688) q[2];
sx q[2];
rz(-0.0071439925) q[2];
sx q[2];
rz(-0.76772493) q[2];
rz(1.399562) q[3];
sx q[3];
rz(-0.00082409516) q[3];
sx q[3];
rz(2.6373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117699) q[0];
sx q[0];
rz(-0.98291021) q[0];
sx q[0];
rz(1.7194189) q[0];
rz(3.1172251) q[1];
sx q[1];
rz(-2.9822646) q[1];
sx q[1];
rz(0.23039625) q[1];
rz(0.89543912) q[2];
sx q[2];
rz(-1.2661305) q[2];
sx q[2];
rz(0.2998395) q[2];
rz(-1.4735994) q[3];
sx q[3];
rz(-1.1553649) q[3];
sx q[3];
rz(-1.8027007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
