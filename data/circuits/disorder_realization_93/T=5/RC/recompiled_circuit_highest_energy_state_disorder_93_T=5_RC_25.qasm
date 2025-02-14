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
rz(2.6858202) q[0];
sx q[0];
rz(-2.0210285) q[0];
sx q[0];
rz(1.7280581) q[0];
rz(-2.7780374) q[1];
sx q[1];
rz(-1.8283365) q[1];
sx q[1];
rz(0.29120905) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.04714) q[0];
sx q[0];
rz(-1.5507076) q[0];
sx q[0];
rz(1.8909341) q[0];
x q[1];
rz(1.2313263) q[2];
sx q[2];
rz(-1.6874773) q[2];
sx q[2];
rz(-1.9872023) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5346966) q[1];
sx q[1];
rz(-1.7217858) q[1];
sx q[1];
rz(0.57735898) q[1];
rz(-pi) q[2];
x q[2];
rz(2.324027) q[3];
sx q[3];
rz(-1.0888087) q[3];
sx q[3];
rz(0.40460247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1353961) q[2];
sx q[2];
rz(-1.5781032) q[2];
sx q[2];
rz(3.1175933) q[2];
rz(2.7855347) q[3];
sx q[3];
rz(-2.4672716) q[3];
sx q[3];
rz(-3.1168028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559219) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(0.63496494) q[0];
rz(-2.3786646) q[1];
sx q[1];
rz(-1.4632016) q[1];
sx q[1];
rz(2.2244804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35058366) q[0];
sx q[0];
rz(-1.6054594) q[0];
sx q[0];
rz(-3.123218) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6659662) q[2];
sx q[2];
rz(-0.94671072) q[2];
sx q[2];
rz(-1.2211354) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8134384) q[1];
sx q[1];
rz(-1.1309237) q[1];
sx q[1];
rz(-3.0564223) q[1];
rz(-pi) q[2];
rz(-1.8648326) q[3];
sx q[3];
rz(-0.64084478) q[3];
sx q[3];
rz(-1.8185735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1941173) q[2];
sx q[2];
rz(-0.52299356) q[2];
sx q[2];
rz(2.4083162) q[2];
rz(2.0745847) q[3];
sx q[3];
rz(-1.4209483) q[3];
sx q[3];
rz(2.9616621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.5697524) q[0];
sx q[0];
rz(-1.1886007) q[0];
sx q[0];
rz(0.52823129) q[0];
rz(-2.275548) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(2.6285062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24241867) q[0];
sx q[0];
rz(-2.0089825) q[0];
sx q[0];
rz(-3.0946381) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31103525) q[2];
sx q[2];
rz(-0.94488178) q[2];
sx q[2];
rz(-3.1300822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21973038) q[1];
sx q[1];
rz(-0.37223909) q[1];
sx q[1];
rz(-1.9277186) q[1];
rz(-pi) q[2];
rz(0.49969466) q[3];
sx q[3];
rz(-1.5958188) q[3];
sx q[3];
rz(-0.5145038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1055866) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(3.1028683) q[2];
rz(0.45768467) q[3];
sx q[3];
rz(-0.80515146) q[3];
sx q[3];
rz(1.3409486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7078581) q[0];
sx q[0];
rz(-2.7271294) q[0];
sx q[0];
rz(-1.1210972) q[0];
rz(2.3838249) q[1];
sx q[1];
rz(-1.517375) q[1];
sx q[1];
rz(-0.68414348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30267495) q[0];
sx q[0];
rz(-0.7644628) q[0];
sx q[0];
rz(1.8643127) q[0];
rz(-pi) q[1];
rz(1.563638) q[2];
sx q[2];
rz(-1.560692) q[2];
sx q[2];
rz(-1.9310538) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0348548) q[1];
sx q[1];
rz(-1.8632006) q[1];
sx q[1];
rz(-0.3180983) q[1];
x q[2];
rz(-2.1675893) q[3];
sx q[3];
rz(-1.1349548) q[3];
sx q[3];
rz(1.3485933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6571558) q[2];
sx q[2];
rz(-1.5269205) q[2];
sx q[2];
rz(-0.26051513) q[2];
rz(0.83812964) q[3];
sx q[3];
rz(-1.9939491) q[3];
sx q[3];
rz(-0.014613541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92358661) q[0];
sx q[0];
rz(-1.3777233) q[0];
sx q[0];
rz(0.51061428) q[0];
rz(-1.8922197) q[1];
sx q[1];
rz(-1.9793972) q[1];
sx q[1];
rz(1.4543264) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.175486) q[0];
sx q[0];
rz(-1.5425872) q[0];
sx q[0];
rz(2.9278446) q[0];
rz(-0.13975039) q[2];
sx q[2];
rz(-1.6183766) q[2];
sx q[2];
rz(2.6222677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0328503) q[1];
sx q[1];
rz(-1.3128442) q[1];
sx q[1];
rz(-0.56405073) q[1];
rz(-pi) q[2];
rz(-0.84564836) q[3];
sx q[3];
rz(-0.21580869) q[3];
sx q[3];
rz(-0.25419054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.55383596) q[2];
sx q[2];
rz(-0.38148701) q[2];
sx q[2];
rz(-0.39056632) q[2];
rz(0.25911123) q[3];
sx q[3];
rz(-1.6389537) q[3];
sx q[3];
rz(0.34671569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7873586) q[0];
sx q[0];
rz(-0.43208313) q[0];
sx q[0];
rz(0.64467347) q[0];
rz(-1.3020172) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(-2.7929746) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88953274) q[0];
sx q[0];
rz(-0.67399287) q[0];
sx q[0];
rz(-0.78788449) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9863527) q[2];
sx q[2];
rz(-2.0698498) q[2];
sx q[2];
rz(2.3504194) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.28305123) q[1];
sx q[1];
rz(-1.9894454) q[1];
sx q[1];
rz(2.8178697) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8556267) q[3];
sx q[3];
rz(-3.0197201) q[3];
sx q[3];
rz(-1.414145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1031441) q[2];
sx q[2];
rz(-0.95547533) q[2];
sx q[2];
rz(3.1032584) q[2];
rz(-2.1740055) q[3];
sx q[3];
rz(-1.7318232) q[3];
sx q[3];
rz(-0.35965317) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2751145) q[0];
sx q[0];
rz(-0.51946467) q[0];
sx q[0];
rz(3.0679829) q[0];
rz(1.7346802) q[1];
sx q[1];
rz(-0.55714566) q[1];
sx q[1];
rz(-2.4085192) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.276537) q[0];
sx q[0];
rz(-2.3650794) q[0];
sx q[0];
rz(0.76188253) q[0];
rz(0.91821155) q[2];
sx q[2];
rz(-0.84647734) q[2];
sx q[2];
rz(1.0670964) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8658856) q[1];
sx q[1];
rz(-1.2162677) q[1];
sx q[1];
rz(1.9237296) q[1];
rz(-pi) q[2];
rz(-1.8215976) q[3];
sx q[3];
rz(-1.7989252) q[3];
sx q[3];
rz(-3.1056014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1072032) q[2];
sx q[2];
rz(-1.4185536) q[2];
sx q[2];
rz(-0.011073152) q[2];
rz(-2.8360046) q[3];
sx q[3];
rz(-2.0162851) q[3];
sx q[3];
rz(1.8886458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3333862) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(2.407684) q[0];
rz(2.5607196) q[1];
sx q[1];
rz(-2.1323233) q[1];
sx q[1];
rz(0.098800585) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042511333) q[0];
sx q[0];
rz(-0.83957273) q[0];
sx q[0];
rz(0.73584105) q[0];
rz(-2.4634741) q[2];
sx q[2];
rz(-2.0653915) q[2];
sx q[2];
rz(-1.8910318) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72291127) q[1];
sx q[1];
rz(-1.1288112) q[1];
sx q[1];
rz(1.2936472) q[1];
rz(-1.4905761) q[3];
sx q[3];
rz(-1.7003683) q[3];
sx q[3];
rz(0.097565325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29844478) q[2];
sx q[2];
rz(-1.2848022) q[2];
sx q[2];
rz(-2.0419545) q[2];
rz(-3.1067276) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(0.59665027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027503969) q[0];
sx q[0];
rz(-1.1570258) q[0];
sx q[0];
rz(1.9644568) q[0];
rz(-1.6674532) q[1];
sx q[1];
rz(-2.2087704) q[1];
sx q[1];
rz(-1.4449545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53430492) q[0];
sx q[0];
rz(-1.8518847) q[0];
sx q[0];
rz(0.40977733) q[0];
x q[1];
rz(-0.99924008) q[2];
sx q[2];
rz(-2.2543467) q[2];
sx q[2];
rz(-2.7627166) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34945378) q[1];
sx q[1];
rz(-1.0568403) q[1];
sx q[1];
rz(-0.07899125) q[1];
rz(1.2304495) q[3];
sx q[3];
rz(-0.7045463) q[3];
sx q[3];
rz(-1.6753538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16704796) q[2];
sx q[2];
rz(-1.5571152) q[2];
sx q[2];
rz(-0.17189279) q[2];
rz(0.78957549) q[3];
sx q[3];
rz(-1.1372477) q[3];
sx q[3];
rz(-1.3353698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0563141) q[0];
sx q[0];
rz(-0.21509898) q[0];
sx q[0];
rz(2.6000182) q[0];
rz(-1.9821292) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(-2.3084739) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43521066) q[0];
sx q[0];
rz(-2.4738389) q[0];
sx q[0];
rz(-2.6866962) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6037349) q[2];
sx q[2];
rz(-0.61213697) q[2];
sx q[2];
rz(-0.15370856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7685582) q[1];
sx q[1];
rz(-1.7972065) q[1];
sx q[1];
rz(0.86707244) q[1];
rz(-pi) q[2];
rz(-2.2984357) q[3];
sx q[3];
rz(-1.5693446) q[3];
sx q[3];
rz(-0.42556083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1200166) q[2];
sx q[2];
rz(-0.77777255) q[2];
sx q[2];
rz(0.94338256) q[2];
rz(2.0639482) q[3];
sx q[3];
rz(-1.5543944) q[3];
sx q[3];
rz(0.32850346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42644603) q[0];
sx q[0];
rz(-2.1052512) q[0];
sx q[0];
rz(-0.25682009) q[0];
rz(2.9019451) q[1];
sx q[1];
rz(-0.9050723) q[1];
sx q[1];
rz(-1.1368652) q[1];
rz(-2.4478365) q[2];
sx q[2];
rz(-2.0907331) q[2];
sx q[2];
rz(1.5246684) q[2];
rz(1.8244459) q[3];
sx q[3];
rz(-1.2696878) q[3];
sx q[3];
rz(2.9326143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
