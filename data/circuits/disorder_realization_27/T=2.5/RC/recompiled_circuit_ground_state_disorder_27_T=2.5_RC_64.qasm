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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12111353) q[0];
sx q[0];
rz(-1.6389567) q[0];
sx q[0];
rz(-1.6421374) q[0];
rz(-0.53977592) q[2];
sx q[2];
rz(-1.9775043) q[2];
sx q[2];
rz(-1.1193898) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1573616) q[1];
sx q[1];
rz(-1.5811635) q[1];
sx q[1];
rz(1.6073411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9755483) q[3];
sx q[3];
rz(-0.18530986) q[3];
sx q[3];
rz(2.6626793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2157669) q[2];
sx q[2];
rz(-0.01130686) q[2];
sx q[2];
rz(-1.0781778) q[2];
rz(-0.82873851) q[3];
sx q[3];
rz(-1.5107892) q[3];
sx q[3];
rz(0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0593798) q[0];
sx q[0];
rz(-1.8635211) q[0];
sx q[0];
rz(1.7510121) q[0];
rz(-1.7104205) q[1];
sx q[1];
rz(-3.1371959) q[1];
sx q[1];
rz(-1.4336525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4162671) q[0];
sx q[0];
rz(-2.4967589) q[0];
sx q[0];
rz(-1.3960442) q[0];
x q[1];
rz(-2.0145871) q[2];
sx q[2];
rz(-1.5559042) q[2];
sx q[2];
rz(-0.040497517) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.28273496) q[1];
sx q[1];
rz(-1.5703859) q[1];
sx q[1];
rz(1.5617227) q[1];
x q[2];
rz(-1.5052027) q[3];
sx q[3];
rz(-0.7842614) q[3];
sx q[3];
rz(-2.0900871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37890515) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(1.572466) q[2];
rz(2.3804741) q[3];
sx q[3];
rz(-0.053839024) q[3];
sx q[3];
rz(1.2083453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.222027) q[0];
sx q[0];
rz(-0.61028218) q[0];
sx q[0];
rz(0.56184226) q[0];
rz(1.5678844) q[1];
sx q[1];
rz(-1.5627197) q[1];
sx q[1];
rz(-3.124253) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0837096) q[0];
sx q[0];
rz(-2.0590933) q[0];
sx q[0];
rz(0.47476193) q[0];
x q[1];
rz(0.88490136) q[2];
sx q[2];
rz(-1.1807311) q[2];
sx q[2];
rz(1.6194413) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8120809) q[1];
sx q[1];
rz(-1.9332262) q[1];
sx q[1];
rz(-0.017048841) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6312977) q[3];
sx q[3];
rz(-1.5908352) q[3];
sx q[3];
rz(-0.35773548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6277546) q[2];
sx q[2];
rz(-1.7226115) q[2];
sx q[2];
rz(-2.5886152) q[2];
rz(-1.2051469) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(-1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39128458) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(1.3543825) q[0];
rz(1.3722108) q[1];
sx q[1];
rz(-3.1387699) q[1];
sx q[1];
rz(1.7599958) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8244265) q[0];
sx q[0];
rz(-1.9726511) q[0];
sx q[0];
rz(2.0749932) q[0];
rz(3.1388248) q[2];
sx q[2];
rz(-1.5684897) q[2];
sx q[2];
rz(3.0879813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.83643944) q[1];
sx q[1];
rz(-2.2341995) q[1];
sx q[1];
rz(2.7022578) q[1];
rz(-pi) q[2];
rz(-2.3390807) q[3];
sx q[3];
rz(-1.3617235) q[3];
sx q[3];
rz(1.0788364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0682721) q[2];
sx q[2];
rz(-0.019651532) q[2];
sx q[2];
rz(-1.9015296) q[2];
rz(-2.9114919) q[3];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087990046) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(1.8863652) q[0];
rz(-0.0093731006) q[1];
sx q[1];
rz(-1.3690989) q[1];
sx q[1];
rz(-0.031551687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2627336) q[0];
sx q[0];
rz(-0.53476214) q[0];
sx q[0];
rz(1.6761024) q[0];
rz(-1.8589694) q[2];
sx q[2];
rz(-0.30063486) q[2];
sx q[2];
rz(-0.86821454) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7158791) q[1];
sx q[1];
rz(-1.3495933) q[1];
sx q[1];
rz(-1.0278661) q[1];
x q[2];
rz(-0.82069759) q[3];
sx q[3];
rz(-0.77337556) q[3];
sx q[3];
rz(-0.14062961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8072529) q[2];
sx q[2];
rz(-0.0060609239) q[2];
sx q[2];
rz(0.80017153) q[2];
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
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0873347) q[0];
sx q[0];
rz(-0.15905173) q[0];
sx q[0];
rz(-1.6416838) q[0];
rz(-0.17290393) q[1];
sx q[1];
rz(-3.0992295) q[1];
sx q[1];
rz(0.049887966) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5373851) q[0];
sx q[0];
rz(-2.4695463) q[0];
sx q[0];
rz(1.9792569) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.028333) q[2];
sx q[2];
rz(-2.335603) q[2];
sx q[2];
rz(-1.9085371) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3499455) q[1];
sx q[1];
rz(-0.87222067) q[1];
sx q[1];
rz(-0.32731815) q[1];
rz(-pi) q[2];
rz(-0.84953725) q[3];
sx q[3];
rz(-1.153933) q[3];
sx q[3];
rz(0.45826926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.738203) q[2];
sx q[2];
rz(-3.093779) q[2];
sx q[2];
rz(-1.8955463) q[2];
rz(1.7854779) q[3];
sx q[3];
rz(-0.035364371) q[3];
sx q[3];
rz(1.416392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-1.4225381) q[0];
rz(1.2369583) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(0.2027771) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1272498) q[0];
sx q[0];
rz(-2.3115251) q[0];
sx q[0];
rz(1.1237445) q[0];
x q[1];
rz(-0.88411096) q[2];
sx q[2];
rz(-2.1565314) q[2];
sx q[2];
rz(-2.3190829) q[2];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.629414) q[1];
x q[2];
rz(-1.2953128) q[3];
sx q[3];
rz(-2.6138895) q[3];
sx q[3];
rz(-2.3964411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5286336) q[2];
sx q[2];
rz(-0.10047675) q[2];
sx q[2];
rz(-2.4670777) q[2];
rz(1.7965192) q[3];
sx q[3];
rz(-2.9967872) q[3];
sx q[3];
rz(1.8176746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26789185) q[0];
sx q[0];
rz(-2.38509) q[0];
sx q[0];
rz(-2.2896413) q[0];
rz(-2.9497228) q[1];
sx q[1];
rz(-0.012902915) q[1];
sx q[1];
rz(2.875944) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56391956) q[0];
sx q[0];
rz(-1.980459) q[0];
sx q[0];
rz(1.760861) q[0];
rz(-2.526409) q[2];
sx q[2];
rz(-1.8447723) q[2];
sx q[2];
rz(3.0469303) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.354968) q[1];
sx q[1];
rz(-1.5119799) q[1];
sx q[1];
rz(-3.0579964) q[1];
rz(0.55691584) q[3];
sx q[3];
rz(-1.89291) q[3];
sx q[3];
rz(2.7469001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3702281) q[2];
sx q[2];
rz(-0.092265487) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8856186) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(1.0783827) q[0];
rz(1.6400317) q[1];
sx q[1];
rz(-0.20137943) q[1];
sx q[1];
rz(-1.5578425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013951741) q[0];
sx q[0];
rz(-1.5553345) q[0];
sx q[0];
rz(-2.4544793) q[0];
rz(-pi) q[1];
rz(-1.0859162) q[2];
sx q[2];
rz(-0.74046248) q[2];
sx q[2];
rz(-1.452335) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8858151) q[1];
sx q[1];
rz(-0.022279169) q[1];
sx q[1];
rz(-0.15099578) q[1];
rz(-pi) q[2];
rz(-1.1224062) q[3];
sx q[3];
rz(-2.4351154) q[3];
sx q[3];
rz(2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.885159) q[2];
sx q[2];
rz(-0.014545518) q[2];
sx q[2];
rz(0.069615901) q[2];
rz(3.0749248) q[3];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2081864) q[0];
sx q[0];
rz(-1.8649768) q[0];
sx q[0];
rz(0.25190121) q[0];
rz(1.4725641) q[1];
sx q[1];
rz(-0.21569574) q[1];
sx q[1];
rz(3.0676945) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6937249) q[0];
sx q[0];
rz(-2.7694919) q[0];
sx q[0];
rz(1.7663581) q[0];
rz(-pi) q[1];
rz(3.1242351) q[2];
sx q[2];
rz(-1.6021172) q[2];
sx q[2];
rz(0.90085122) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.042797814) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(-1.2513729) q[1];
rz(2.4611887) q[3];
sx q[3];
rz(-2.8814253) q[3];
sx q[3];
rz(2.5870067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4613688) q[2];
sx q[2];
rz(-0.0071439925) q[2];
sx q[2];
rz(2.3738677) q[2];
rz(-1.7420306) q[3];
sx q[3];
rz(-3.1407686) q[3];
sx q[3];
rz(0.50423938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-2.7585898) q[2];
sx q[2];
rz(-0.93180626) q[2];
sx q[2];
rz(2.1064482) q[2];
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
