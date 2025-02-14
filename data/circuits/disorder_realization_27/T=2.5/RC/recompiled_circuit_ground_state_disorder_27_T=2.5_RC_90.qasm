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
rz(0.93552247) q[0];
rz(1.525653) q[1];
sx q[1];
rz(4.6620044) q[1];
sx q[1];
rz(9.6992156) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0204791) q[0];
sx q[0];
rz(-1.6389567) q[0];
sx q[0];
rz(1.6421374) q[0];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-0.86286592) q[1];
sx q[1];
rz(-0.037986156) q[1];
sx q[1];
rz(-1.8472865) q[1];
rz(-1.1660444) q[3];
sx q[3];
rz(-0.18530986) q[3];
sx q[3];
rz(-0.47891339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2157669) q[2];
sx q[2];
rz(-0.01130686) q[2];
sx q[2];
rz(-2.0634148) q[2];
rz(2.3128541) q[3];
sx q[3];
rz(-1.6308035) q[3];
sx q[3];
rz(2.3327648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082212903) q[0];
sx q[0];
rz(-1.2780715) q[0];
sx q[0];
rz(1.7510121) q[0];
rz(-1.7104205) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(1.4336525) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8468877) q[0];
sx q[0];
rz(-1.4661015) q[0];
sx q[0];
rz(-0.93331915) q[0];
rz(-pi) q[1];
rz(1.1270056) q[2];
sx q[2];
rz(-1.5559042) q[2];
sx q[2];
rz(3.1010951) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8588577) q[1];
sx q[1];
rz(-1.5712067) q[1];
sx q[1];
rz(1.57987) q[1];
rz(-pi) q[2];
rz(0.78761132) q[3];
sx q[3];
rz(-1.5244841) q[3];
sx q[3];
rz(-2.5758343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7626875) q[2];
sx q[2];
rz(-1.5451558) q[2];
sx q[2];
rz(1.5691266) q[2];
rz(-2.3804741) q[3];
sx q[3];
rz(-0.053839024) q[3];
sx q[3];
rz(-1.2083453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91956562) q[0];
sx q[0];
rz(-2.5313105) q[0];
sx q[0];
rz(-0.56184226) q[0];
rz(-1.5737083) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(3.124253) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2527232) q[0];
sx q[0];
rz(-0.66735744) q[0];
sx q[0];
rz(-0.86020893) q[0];
rz(0.48835619) q[2];
sx q[2];
rz(-2.1966509) q[2];
sx q[2];
rz(-0.35036119) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.764027) q[1];
sx q[1];
rz(-2.7787797) q[1];
sx q[1];
rz(-1.6157263) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1005834) q[3];
sx q[3];
rz(-2.6309391) q[3];
sx q[3];
rz(1.8927495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6277546) q[2];
sx q[2];
rz(-1.4189812) q[2];
sx q[2];
rz(-0.55297744) q[2];
rz(1.2051469) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(-1.5470101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.39128458) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(-1.3543825) q[0];
rz(1.3722108) q[1];
sx q[1];
rz(-3.1387699) q[1];
sx q[1];
rz(-1.3815968) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1005122) q[0];
sx q[0];
rz(-1.1100475) q[0];
sx q[0];
rz(2.6897088) q[0];
rz(-3.1388248) q[2];
sx q[2];
rz(-1.5684897) q[2];
sx q[2];
rz(0.053611343) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1255252) q[1];
sx q[1];
rz(-1.9125332) q[1];
sx q[1];
rz(0.85846947) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8546811) q[3];
sx q[3];
rz(-2.3182456) q[3];
sx q[3];
rz(0.2940184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0733205) q[2];
sx q[2];
rz(-0.019651532) q[2];
sx q[2];
rz(-1.2400631) q[2];
rz(-0.23010075) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(-0.41964644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(3.0536026) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(1.8863652) q[0];
rz(3.1322196) q[1];
sx q[1];
rz(-1.7724937) q[1];
sx q[1];
rz(-3.110041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.878859) q[0];
sx q[0];
rz(-0.53476214) q[0];
sx q[0];
rz(-1.6761024) q[0];
x q[1];
rz(-1.8597262) q[2];
sx q[2];
rz(-1.4865371) q[2];
sx q[2];
rz(0.97848985) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6477709) q[1];
sx q[1];
rz(-2.5595409) q[1];
sx q[1];
rz(1.1602559) q[1];
rz(-pi) q[2];
x q[2];
rz(2.191014) q[3];
sx q[3];
rz(-1.0744541) q[3];
sx q[3];
rz(0.84211189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3343398) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(-0.80017153) q[2];
rz(2.3470894) q[3];
sx q[3];
rz(-3.1092293) q[3];
sx q[3];
rz(-2.2728424) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.054258) q[0];
sx q[0];
rz(-0.15905173) q[0];
sx q[0];
rz(-1.4999088) q[0];
rz(2.9686887) q[1];
sx q[1];
rz(-0.042363107) q[1];
sx q[1];
rz(-0.049887966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5373851) q[0];
sx q[0];
rz(-2.4695463) q[0];
sx q[0];
rz(-1.1623357) q[0];
rz(-pi) q[1];
rz(-0.11325963) q[2];
sx q[2];
rz(-2.335603) q[2];
sx q[2];
rz(-1.2330556) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0058635423) q[1];
sx q[1];
rz(-1.3220437) q[1];
sx q[1];
rz(-2.2963524) q[1];
x q[2];
rz(0.84953725) q[3];
sx q[3];
rz(-1.153933) q[3];
sx q[3];
rz(2.6833234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40338966) q[2];
sx q[2];
rz(-0.047813606) q[2];
sx q[2];
rz(-1.8955463) q[2];
rz(1.3561148) q[3];
sx q[3];
rz(-3.1062283) q[3];
sx q[3];
rz(1.416392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-1.7127011) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(-1.4225381) q[0];
rz(1.2369583) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(-2.9388156) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13064676) q[0];
sx q[0];
rz(-1.8955064) q[0];
sx q[0];
rz(0.79239158) q[0];
rz(0.70901633) q[2];
sx q[2];
rz(-2.1273899) q[2];
sx q[2];
rz(-2.8188044) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8824892) q[1];
sx q[1];
rz(-1.9111425) q[1];
sx q[1];
rz(1.5121786) q[1];
rz(-pi) q[2];
rz(0.15722991) q[3];
sx q[3];
rz(-2.0766602) q[3];
sx q[3];
rz(-1.0613522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5286336) q[2];
sx q[2];
rz(-0.10047675) q[2];
sx q[2];
rz(-0.67451492) q[2];
rz(1.3450735) q[3];
sx q[3];
rz(-2.9967872) q[3];
sx q[3];
rz(1.323918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26789185) q[0];
sx q[0];
rz(-0.75650263) q[0];
sx q[0];
rz(2.2896413) q[0];
rz(0.1918699) q[1];
sx q[1];
rz(-0.012902915) q[1];
sx q[1];
rz(2.875944) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56391956) q[0];
sx q[0];
rz(-1.1611337) q[0];
sx q[0];
rz(1.760861) q[0];
x q[1];
rz(2.688411) q[2];
sx q[2];
rz(-0.66614775) q[2];
sx q[2];
rz(-2.0311126) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78662465) q[1];
sx q[1];
rz(-1.5119799) q[1];
sx q[1];
rz(-0.083596283) q[1];
rz(1.1962131) q[3];
sx q[3];
rz(-1.0455971) q[3];
sx q[3];
rz(1.7708667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77136451) q[2];
sx q[2];
rz(-0.092265487) q[2];
sx q[2];
rz(-0.51221687) q[2];
rz(-2.9825315) q[3];
sx q[3];
rz(-3.1064807) q[3];
sx q[3];
rz(1.7062645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2559741) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(1.0783827) q[0];
rz(-1.501561) q[1];
sx q[1];
rz(-0.20137943) q[1];
sx q[1];
rz(1.5837502) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5974332) q[0];
sx q[0];
rz(-0.88378105) q[0];
sx q[0];
rz(-1.5907955) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7388849) q[2];
sx q[2];
rz(-0.93120775) q[2];
sx q[2];
rz(2.072203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4068102) q[1];
sx q[1];
rz(-1.592822) q[1];
sx q[1];
rz(-1.5741482) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0191865) q[3];
sx q[3];
rz(-0.70647722) q[3];
sx q[3];
rz(-2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25643361) q[2];
sx q[2];
rz(-3.1270471) q[2];
sx q[2];
rz(0.069615901) q[2];
rz(-3.0749248) q[3];
sx q[3];
rz(-2.127141) q[3];
sx q[3];
rz(-0.67924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93340623) q[0];
sx q[0];
rz(-1.8649768) q[0];
sx q[0];
rz(-2.8896914) q[0];
rz(-1.6690286) q[1];
sx q[1];
rz(-0.21569574) q[1];
sx q[1];
rz(3.0676945) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6937249) q[0];
sx q[0];
rz(-0.37210074) q[0];
sx q[0];
rz(1.3752346) q[0];
x q[1];
rz(-1.5394707) q[2];
sx q[2];
rz(-1.5534473) q[2];
sx q[2];
rz(0.67048873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.51346362) q[1];
sx q[1];
rz(-2.3766915) q[1];
sx q[1];
rz(-2.7897863) q[1];
x q[2];
rz(2.9375495) q[3];
sx q[3];
rz(-1.7333441) q[3];
sx q[3];
rz(-2.7891087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4613688) q[2];
sx q[2];
rz(-0.0071439925) q[2];
sx q[2];
rz(2.3738677) q[2];
rz(1.399562) q[3];
sx q[3];
rz(-3.1407686) q[3];
sx q[3];
rz(0.50423938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(-3.1172251) q[1];
sx q[1];
rz(-0.15932803) q[1];
sx q[1];
rz(-2.9111964) q[1];
rz(-0.89543912) q[2];
sx q[2];
rz(-1.8754621) q[2];
sx q[2];
rz(-2.8417532) q[2];
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
