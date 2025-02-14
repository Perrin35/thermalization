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
rz(2.7316982) q[0];
sx q[0];
rz(-1.7562261) q[0];
sx q[0];
rz(-1.0710427) q[0];
rz(1.8154124) q[1];
sx q[1];
rz(-0.92547995) q[1];
sx q[1];
rz(-0.34223908) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5919094) q[0];
sx q[0];
rz(-0.91381493) q[0];
sx q[0];
rz(2.6651938) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1360004) q[2];
sx q[2];
rz(-2.2965429) q[2];
sx q[2];
rz(1.0623887) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46821478) q[1];
sx q[1];
rz(-1.6771375) q[1];
sx q[1];
rz(2.0061226) q[1];
rz(-pi) q[2];
rz(1.6116404) q[3];
sx q[3];
rz(-2.7451635) q[3];
sx q[3];
rz(1.886823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.60275117) q[2];
sx q[2];
rz(-2.4138236) q[2];
sx q[2];
rz(1.7068498) q[2];
rz(0.96628609) q[3];
sx q[3];
rz(-1.1937701) q[3];
sx q[3];
rz(2.5626591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12427881) q[0];
sx q[0];
rz(-0.74324981) q[0];
sx q[0];
rz(0.055835128) q[0];
rz(2.3827379) q[1];
sx q[1];
rz(-2.0003624) q[1];
sx q[1];
rz(-0.32078823) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43335981) q[0];
sx q[0];
rz(-0.77455168) q[0];
sx q[0];
rz(-0.21277986) q[0];
rz(-0.35050138) q[2];
sx q[2];
rz(-2.3308518) q[2];
sx q[2];
rz(1.8271918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1244679) q[1];
sx q[1];
rz(-1.6001646) q[1];
sx q[1];
rz(-1.8046124) q[1];
rz(-pi) q[2];
rz(2.6296162) q[3];
sx q[3];
rz(-0.80080355) q[3];
sx q[3];
rz(-1.2146571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.32030973) q[2];
sx q[2];
rz(-1.2592969) q[2];
sx q[2];
rz(-2.8584976) q[2];
rz(-0.29575944) q[3];
sx q[3];
rz(-1.7167973) q[3];
sx q[3];
rz(-1.504771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-0.54842424) q[0];
sx q[0];
rz(-1.230509) q[0];
sx q[0];
rz(3.044627) q[0];
rz(-1.4292258) q[1];
sx q[1];
rz(-1.2173419) q[1];
sx q[1];
rz(0.355535) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.36847) q[0];
sx q[0];
rz(-2.5036044) q[0];
sx q[0];
rz(0.49686749) q[0];
x q[1];
rz(-2.7266302) q[2];
sx q[2];
rz(-1.4769161) q[2];
sx q[2];
rz(2.0866657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1946845) q[1];
sx q[1];
rz(-0.42153639) q[1];
sx q[1];
rz(1.3409019) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16014647) q[3];
sx q[3];
rz(-2.886842) q[3];
sx q[3];
rz(-1.3046978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7045483) q[2];
sx q[2];
rz(-1.9883678) q[2];
sx q[2];
rz(2.8538749) q[2];
rz(-2.8904166) q[3];
sx q[3];
rz(-1.0947248) q[3];
sx q[3];
rz(1.2742961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9525725) q[0];
sx q[0];
rz(-2.548521) q[0];
sx q[0];
rz(0.75743842) q[0];
rz(-0.82713741) q[1];
sx q[1];
rz(-0.65991455) q[1];
sx q[1];
rz(1.8338902) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73440512) q[0];
sx q[0];
rz(-1.4467753) q[0];
sx q[0];
rz(-3.0888686) q[0];
rz(-pi) q[1];
x q[1];
rz(2.04129) q[2];
sx q[2];
rz(-1.2963106) q[2];
sx q[2];
rz(-0.16484552) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8825619) q[1];
sx q[1];
rz(-0.81444959) q[1];
sx q[1];
rz(-0.37268727) q[1];
rz(-pi) q[2];
rz(1.3321628) q[3];
sx q[3];
rz(-0.25128579) q[3];
sx q[3];
rz(-0.3005614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7149675) q[2];
sx q[2];
rz(-1.7968618) q[2];
sx q[2];
rz(-1.3955383) q[2];
rz(1.7914145) q[3];
sx q[3];
rz(-1.6387286) q[3];
sx q[3];
rz(-0.028039886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79513079) q[0];
sx q[0];
rz(-2.6030354) q[0];
sx q[0];
rz(-1.7002456) q[0];
rz(-0.92789188) q[1];
sx q[1];
rz(-2.3065232) q[1];
sx q[1];
rz(-0.79383129) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2071811) q[0];
sx q[0];
rz(-2.7156805) q[0];
sx q[0];
rz(1.8239198) q[0];
rz(-1.687811) q[2];
sx q[2];
rz(-1.4037366) q[2];
sx q[2];
rz(0.88167215) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87463996) q[1];
sx q[1];
rz(-0.62332223) q[1];
sx q[1];
rz(-2.1165352) q[1];
rz(-pi) q[2];
rz(-0.50516523) q[3];
sx q[3];
rz(-1.9117711) q[3];
sx q[3];
rz(-2.7974432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64357197) q[2];
sx q[2];
rz(-1.4683495) q[2];
sx q[2];
rz(-2.7872861) q[2];
rz(1.0962567) q[3];
sx q[3];
rz(-0.26429629) q[3];
sx q[3];
rz(1.8335584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7588014) q[0];
sx q[0];
rz(-2.379874) q[0];
sx q[0];
rz(-2.288901) q[0];
rz(-0.26607749) q[1];
sx q[1];
rz(-2.1177025) q[1];
sx q[1];
rz(-1.3583677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32222267) q[0];
sx q[0];
rz(-2.529105) q[0];
sx q[0];
rz(1.9507381) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2439315) q[2];
sx q[2];
rz(-1.7836469) q[2];
sx q[2];
rz(2.0124679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0315288) q[1];
sx q[1];
rz(-2.4569521) q[1];
sx q[1];
rz(-2.5972874) q[1];
rz(-pi) q[2];
rz(1.2262906) q[3];
sx q[3];
rz(-2.4427569) q[3];
sx q[3];
rz(1.4188053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69409662) q[2];
sx q[2];
rz(-1.2022377) q[2];
sx q[2];
rz(0.48074943) q[2];
rz(-0.9350183) q[3];
sx q[3];
rz(-0.98116773) q[3];
sx q[3];
rz(2.7267314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9347436) q[0];
sx q[0];
rz(-0.54731363) q[0];
sx q[0];
rz(0.68688399) q[0];
rz(-2.1737449) q[1];
sx q[1];
rz(-1.5279488) q[1];
sx q[1];
rz(-2.9639249) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0570014) q[0];
sx q[0];
rz(-0.54750204) q[0];
sx q[0];
rz(-2.4390319) q[0];
rz(-2.4685079) q[2];
sx q[2];
rz(-2.1161352) q[2];
sx q[2];
rz(2.2506491) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8502151) q[1];
sx q[1];
rz(-1.134877) q[1];
sx q[1];
rz(-2.8368783) q[1];
rz(0.92063825) q[3];
sx q[3];
rz(-3.0282109) q[3];
sx q[3];
rz(-2.6760898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0873506) q[2];
sx q[2];
rz(-0.44173104) q[2];
sx q[2];
rz(2.1596215) q[2];
rz(-0.25518498) q[3];
sx q[3];
rz(-1.1533777) q[3];
sx q[3];
rz(2.065778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0690689) q[0];
sx q[0];
rz(-0.60966063) q[0];
sx q[0];
rz(-0.12741086) q[0];
rz(1.6690856) q[1];
sx q[1];
rz(-1.9576879) q[1];
sx q[1];
rz(-1.4471819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0873647) q[0];
sx q[0];
rz(-1.3854829) q[0];
sx q[0];
rz(-2.8186356) q[0];
rz(-2.5650756) q[2];
sx q[2];
rz(-0.62453068) q[2];
sx q[2];
rz(1.3006398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6622611) q[1];
sx q[1];
rz(-1.4331967) q[1];
sx q[1];
rz(1.496973) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98663892) q[3];
sx q[3];
rz(-1.7800864) q[3];
sx q[3];
rz(-1.926146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.24826) q[2];
sx q[2];
rz(-0.97951639) q[2];
sx q[2];
rz(-0.40929201) q[2];
rz(-2.4457757) q[3];
sx q[3];
rz(-2.4455363) q[3];
sx q[3];
rz(2.5223562) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5291587) q[0];
sx q[0];
rz(-2.9127064) q[0];
sx q[0];
rz(2.9869475) q[0];
rz(2.0507428) q[1];
sx q[1];
rz(-0.88328528) q[1];
sx q[1];
rz(-1.8152016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94448254) q[0];
sx q[0];
rz(-1.8354699) q[0];
sx q[0];
rz(-0.84828429) q[0];
x q[1];
rz(-1.0716076) q[2];
sx q[2];
rz(-1.6130924) q[2];
sx q[2];
rz(2.7214839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1005724) q[1];
sx q[1];
rz(-2.3902767) q[1];
sx q[1];
rz(-1.643341) q[1];
x q[2];
rz(-2.7560968) q[3];
sx q[3];
rz(-0.49043819) q[3];
sx q[3];
rz(0.78133327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.078330127) q[2];
sx q[2];
rz(-0.41389725) q[2];
sx q[2];
rz(2.4388893) q[2];
rz(3.0472158) q[3];
sx q[3];
rz(-0.97797147) q[3];
sx q[3];
rz(0.92488658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6150045) q[0];
sx q[0];
rz(-2.1115392) q[0];
sx q[0];
rz(-2.7870542) q[0];
rz(-1.7189369) q[1];
sx q[1];
rz(-1.4965897) q[1];
sx q[1];
rz(2.6197701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3799135) q[0];
sx q[0];
rz(-1.9229793) q[0];
sx q[0];
rz(-2.9814979) q[0];
rz(-0.13846878) q[2];
sx q[2];
rz(-1.2684141) q[2];
sx q[2];
rz(2.2669528) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21510151) q[1];
sx q[1];
rz(-2.1128078) q[1];
sx q[1];
rz(0.78185977) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4168968) q[3];
sx q[3];
rz(-1.0130289) q[3];
sx q[3];
rz(-2.9036221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78174138) q[2];
sx q[2];
rz(-1.968911) q[2];
sx q[2];
rz(-2.4594128) q[2];
rz(1.3683246) q[3];
sx q[3];
rz(-1.5409639) q[3];
sx q[3];
rz(2.1413596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2405887) q[0];
sx q[0];
rz(-0.95745845) q[0];
sx q[0];
rz(-0.29314713) q[0];
rz(2.3121569) q[1];
sx q[1];
rz(-1.71143) q[1];
sx q[1];
rz(-1.5477187) q[1];
rz(-1.8315567) q[2];
sx q[2];
rz(-1.793486) q[2];
sx q[2];
rz(-1.7785265) q[2];
rz(0.10503332) q[3];
sx q[3];
rz(-1.5337957) q[3];
sx q[3];
rz(0.93883088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
