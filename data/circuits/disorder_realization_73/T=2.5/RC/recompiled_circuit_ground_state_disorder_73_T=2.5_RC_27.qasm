OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(-2.8889416) q[0];
sx q[0];
rz(1.2055612) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(4.4681273) q[1];
sx q[1];
rz(11.820643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4843895) q[0];
sx q[0];
rz(-1.6307388) q[0];
sx q[0];
rz(1.7933229) q[0];
rz(-pi) q[1];
rz(1.5775852) q[2];
sx q[2];
rz(-1.3464084) q[2];
sx q[2];
rz(-0.16358384) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2066951) q[1];
sx q[1];
rz(-1.7803915) q[1];
sx q[1];
rz(-1.2278201) q[1];
rz(-pi) q[2];
rz(1.0966572) q[3];
sx q[3];
rz(-1.385194) q[3];
sx q[3];
rz(1.4193648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2241609) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(1.6708299) q[2];
rz(-1.6338232) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(-0.9233709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(1.5731328) q[0];
rz(0.16954999) q[1];
sx q[1];
rz(-3.0270271) q[1];
sx q[1];
rz(-0.13470185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8226263) q[0];
sx q[0];
rz(-1.9958226) q[0];
sx q[0];
rz(-1.951773) q[0];
rz(-1.5021654) q[2];
sx q[2];
rz(-1.5830212) q[2];
sx q[2];
rz(-3.1032012) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.811482) q[1];
sx q[1];
rz(-2.5977784) q[1];
sx q[1];
rz(0.61119975) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.991113) q[3];
sx q[3];
rz(-2.947071) q[3];
sx q[3];
rz(-1.9357301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0160825) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(0.15277319) q[2];
rz(1.3656535) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(-2.9735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50853866) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(-2.659944) q[0];
rz(2.956849) q[1];
sx q[1];
rz(-1.3667204) q[1];
sx q[1];
rz(-0.99536037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48993123) q[0];
sx q[0];
rz(-1.7999987) q[0];
sx q[0];
rz(1.9263173) q[0];
rz(-pi) q[1];
rz(-0.048599343) q[2];
sx q[2];
rz(-1.6311797) q[2];
sx q[2];
rz(-0.27602613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0479483) q[1];
sx q[1];
rz(-1.8694832) q[1];
sx q[1];
rz(2.4574952) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5257201) q[3];
sx q[3];
rz(-0.59642506) q[3];
sx q[3];
rz(-3.0395122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92622009) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(-0.064662956) q[2];
rz(2.0926545) q[3];
sx q[3];
rz(-3.1147396) q[3];
sx q[3];
rz(1.8612727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8365086) q[0];
sx q[0];
rz(-0.11798141) q[0];
sx q[0];
rz(-2.3205561) q[0];
rz(3.0631284) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(-0.99748126) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4560735) q[0];
sx q[0];
rz(-0.711383) q[0];
sx q[0];
rz(2.6472241) q[0];
rz(-pi) q[1];
rz(-1.6299963) q[2];
sx q[2];
rz(-1.534605) q[2];
sx q[2];
rz(2.9179887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.32179896) q[1];
sx q[1];
rz(-0.60056486) q[1];
sx q[1];
rz(-1.0933881) q[1];
rz(-pi) q[2];
rz(-1.5134345) q[3];
sx q[3];
rz(-1.093075) q[3];
sx q[3];
rz(0.78335947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0492101) q[2];
sx q[2];
rz(-3.1123078) q[2];
sx q[2];
rz(1.9878261) q[2];
rz(-2.9554101) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(0.33128273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368197) q[0];
sx q[0];
rz(-2.3240219) q[0];
sx q[0];
rz(1.1432884) q[0];
rz(-1.1413057) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(-2.6236261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9310936) q[0];
sx q[0];
rz(-2.2552975) q[0];
sx q[0];
rz(-0.71523358) q[0];
rz(-pi) q[1];
rz(-1.5557655) q[2];
sx q[2];
rz(-1.5698729) q[2];
sx q[2];
rz(2.8944588) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3822381) q[1];
sx q[1];
rz(-1.8439981) q[1];
sx q[1];
rz(1.9135936) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7180743) q[3];
sx q[3];
rz(-1.275835) q[3];
sx q[3];
rz(1.1671821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3190069) q[2];
sx q[2];
rz(-0.95390445) q[2];
sx q[2];
rz(-2.8138568) q[2];
rz(-1.94708) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(-0.85668844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1389393) q[0];
sx q[0];
rz(-0.26873538) q[0];
sx q[0];
rz(0.56185454) q[0];
rz(-1.6506763) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(0.098310016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033103) q[0];
sx q[0];
rz(-1.9925963) q[0];
sx q[0];
rz(-0.038859239) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1409522) q[2];
sx q[2];
rz(-1.5709236) q[2];
sx q[2];
rz(-1.887111) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5206388) q[1];
sx q[1];
rz(-2.6700182) q[1];
sx q[1];
rz(0.10271272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6410511) q[3];
sx q[3];
rz(-1.1574452) q[3];
sx q[3];
rz(-0.40412892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.085122846) q[2];
sx q[2];
rz(-2.9462908) q[2];
sx q[2];
rz(-2.0758212) q[2];
rz(0.39984518) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(-1.2679509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9999076) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(1.1204002) q[1];
sx q[1];
rz(-0.1381865) q[1];
sx q[1];
rz(-0.33946005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65050478) q[0];
sx q[0];
rz(-1.5575214) q[0];
sx q[0];
rz(3.0646532) q[0];
rz(-pi) q[1];
rz(2.5869114) q[2];
sx q[2];
rz(-1.5825869) q[2];
sx q[2];
rz(1.5761216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.016638674) q[1];
sx q[1];
rz(-1.5691688) q[1];
sx q[1];
rz(-1.6563936) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27642997) q[3];
sx q[3];
rz(-0.35110858) q[3];
sx q[3];
rz(1.029226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7772943) q[2];
sx q[2];
rz(-3.1396781) q[2];
sx q[2];
rz(0.36336362) q[2];
rz(1.0904788) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(-2.0232078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(0.61984396) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(-0.92292619) q[0];
rz(1.482831) q[1];
sx q[1];
rz(-0.62983477) q[1];
sx q[1];
rz(-0.0042075687) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18989604) q[0];
sx q[0];
rz(-1.0605264) q[0];
sx q[0];
rz(0.6999335) q[0];
x q[1];
rz(1.9816859) q[2];
sx q[2];
rz(-1.5762323) q[2];
sx q[2];
rz(-1.5856575) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5397268) q[1];
sx q[1];
rz(-1.5063573) q[1];
sx q[1];
rz(1.0947202) q[1];
rz(-1.5245304) q[3];
sx q[3];
rz(-0.95887254) q[3];
sx q[3];
rz(2.9718252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5620455) q[2];
sx q[2];
rz(-1.6041218) q[2];
sx q[2];
rz(-1.9407678) q[2];
rz(-1.3935401) q[3];
sx q[3];
rz(-0.0037007185) q[3];
sx q[3];
rz(0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414108) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(-1.2196983) q[0];
rz(-1.3297184) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(-2.9776998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.208928) q[0];
sx q[0];
rz(-2.744356) q[0];
sx q[0];
rz(-2.0276638) q[0];
x q[1];
rz(2.5570325) q[2];
sx q[2];
rz(-1.8001846) q[2];
sx q[2];
rz(-0.10525119) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5245318) q[1];
sx q[1];
rz(-1.6304029) q[1];
sx q[1];
rz(-1.6888794) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4285942) q[3];
sx q[3];
rz(-0.87942356) q[3];
sx q[3];
rz(-0.52934968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4280052) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(0.62291992) q[2];
rz(3.0981787) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(-2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49122214) q[0];
sx q[0];
rz(-0.0064042052) q[0];
sx q[0];
rz(2.6553335) q[0];
rz(-0.69475118) q[1];
sx q[1];
rz(-0.34725747) q[1];
sx q[1];
rz(0.4756701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2025288) q[0];
sx q[0];
rz(-0.049204218) q[0];
sx q[0];
rz(-0.81728191) q[0];
rz(2.2769502) q[2];
sx q[2];
rz(-1.5448031) q[2];
sx q[2];
rz(-3.1033422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5359395) q[1];
sx q[1];
rz(-1.5772784) q[1];
sx q[1];
rz(-2.5340213) q[1];
rz(2.6127897) q[3];
sx q[3];
rz(-0.35548726) q[3];
sx q[3];
rz(2.5606009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67705578) q[2];
sx q[2];
rz(-0.060364351) q[2];
sx q[2];
rz(-1.2976868) q[2];
rz(-1.248598) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(-0.36778522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1260592) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(-0.88232782) q[1];
sx q[1];
rz(-0.066451646) q[1];
sx q[1];
rz(0.8393504) q[1];
rz(1.5346943) q[2];
sx q[2];
rz(-1.6629728) q[2];
sx q[2];
rz(-2.865553) q[2];
rz(-1.5757379) q[3];
sx q[3];
rz(-1.8099804) q[3];
sx q[3];
rz(-3.119131) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
