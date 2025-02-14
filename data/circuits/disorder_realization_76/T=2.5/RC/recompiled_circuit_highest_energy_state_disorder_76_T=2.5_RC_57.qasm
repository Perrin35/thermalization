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
rz(-2.2657356) q[0];
sx q[0];
rz(-2.1729204) q[0];
sx q[0];
rz(0.92585603) q[0];
rz(0.61697382) q[1];
sx q[1];
rz(8.7595688) q[1];
sx q[1];
rz(8.1002357) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2331794) q[0];
sx q[0];
rz(-0.97199134) q[0];
sx q[0];
rz(1.6480867) q[0];
rz(0.1466897) q[2];
sx q[2];
rz(-2.5518199) q[2];
sx q[2];
rz(1.9751994) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.032925978) q[1];
sx q[1];
rz(-0.98331988) q[1];
sx q[1];
rz(-1.0125748) q[1];
rz(-2.376287) q[3];
sx q[3];
rz(-1.0648784) q[3];
sx q[3];
rz(-1.4368865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8806184) q[2];
sx q[2];
rz(-1.3857434) q[2];
sx q[2];
rz(-0.79208881) q[2];
rz(-0.875862) q[3];
sx q[3];
rz(-0.10871092) q[3];
sx q[3];
rz(0.66332269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51139128) q[0];
sx q[0];
rz(-1.317861) q[0];
sx q[0];
rz(-2.200101) q[0];
rz(0.9990274) q[1];
sx q[1];
rz(-2.2143054) q[1];
sx q[1];
rz(-3.0587382) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1076528) q[0];
sx q[0];
rz(-1.4705188) q[0];
sx q[0];
rz(-0.73579181) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98495667) q[2];
sx q[2];
rz(-2.7212866) q[2];
sx q[2];
rz(-2.676385) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1513843) q[1];
sx q[1];
rz(-2.4308204) q[1];
sx q[1];
rz(0.73020331) q[1];
x q[2];
rz(0.64506809) q[3];
sx q[3];
rz(-2.0760025) q[3];
sx q[3];
rz(0.92407214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.89722172) q[2];
sx q[2];
rz(-1.3542078) q[2];
sx q[2];
rz(-1.8841891) q[2];
rz(2.2952378) q[3];
sx q[3];
rz(-2.9823163) q[3];
sx q[3];
rz(2.4634821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9427247) q[0];
sx q[0];
rz(-1.6833479) q[0];
sx q[0];
rz(-0.22920907) q[0];
rz(-0.21408679) q[1];
sx q[1];
rz(-2.2702718) q[1];
sx q[1];
rz(-0.3814989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5435181) q[0];
sx q[0];
rz(-1.2208573) q[0];
sx q[0];
rz(-2.475481) q[0];
rz(-pi) q[1];
rz(1.7842641) q[2];
sx q[2];
rz(-2.3980283) q[2];
sx q[2];
rz(1.5539757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9126622) q[1];
sx q[1];
rz(-3.0813363) q[1];
sx q[1];
rz(1.0723537) q[1];
rz(-pi) q[2];
rz(-1.5656785) q[3];
sx q[3];
rz(-1.2237142) q[3];
sx q[3];
rz(2.0818613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7233589) q[2];
sx q[2];
rz(-0.25622076) q[2];
sx q[2];
rz(-2.5733433) q[2];
rz(-1.7226284) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(1.0770146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028246183) q[0];
sx q[0];
rz(-2.009511) q[0];
sx q[0];
rz(-0.25076732) q[0];
rz(0.084550683) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(1.2006522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47287175) q[0];
sx q[0];
rz(-0.78470147) q[0];
sx q[0];
rz(-2.0913823) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24742207) q[2];
sx q[2];
rz(-1.1411925) q[2];
sx q[2];
rz(-2.1199494) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1657287) q[1];
sx q[1];
rz(-1.4714676) q[1];
sx q[1];
rz(2.4286859) q[1];
rz(-pi) q[2];
rz(-2.0968998) q[3];
sx q[3];
rz(-1.8013489) q[3];
sx q[3];
rz(-1.5073204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4667929) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(1.2910845) q[2];
rz(2.71991) q[3];
sx q[3];
rz(-2.6899874) q[3];
sx q[3];
rz(-1.5765367) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52294937) q[0];
sx q[0];
rz(-2.4186501) q[0];
sx q[0];
rz(1.500754) q[0];
rz(-0.067642637) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(-1.1281475) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65928946) q[0];
sx q[0];
rz(-1.3829675) q[0];
sx q[0];
rz(-1.6085621) q[0];
rz(-pi) q[1];
rz(2.1352) q[2];
sx q[2];
rz(-0.60205209) q[2];
sx q[2];
rz(2.4133701) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0813898) q[1];
sx q[1];
rz(-0.12559016) q[1];
sx q[1];
rz(-1.1274028) q[1];
x q[2];
rz(-2.2637783) q[3];
sx q[3];
rz(-2.7846309) q[3];
sx q[3];
rz(-1.9138543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1077914) q[2];
sx q[2];
rz(-1.1148323) q[2];
sx q[2];
rz(-0.6768674) q[2];
rz(0.90773165) q[3];
sx q[3];
rz(-1.2264484) q[3];
sx q[3];
rz(2.7270253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21861976) q[0];
sx q[0];
rz(-2.3899879) q[0];
sx q[0];
rz(0.76939097) q[0];
rz(1.9006624) q[1];
sx q[1];
rz(-2.2878094) q[1];
sx q[1];
rz(-2.9262537) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0701987) q[0];
sx q[0];
rz(-0.60766534) q[0];
sx q[0];
rz(-1.9032701) q[0];
rz(-pi) q[1];
x q[1];
rz(1.239993) q[2];
sx q[2];
rz(-1.5869405) q[2];
sx q[2];
rz(-2.6789897) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0907545) q[1];
sx q[1];
rz(-1.3267446) q[1];
sx q[1];
rz(2.4825216) q[1];
rz(-pi) q[2];
rz(1.9086544) q[3];
sx q[3];
rz(-0.90150276) q[3];
sx q[3];
rz(0.026642628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.874959) q[2];
sx q[2];
rz(-1.6338438) q[2];
sx q[2];
rz(0.61666644) q[2];
rz(-0.2002317) q[3];
sx q[3];
rz(-2.4112406) q[3];
sx q[3];
rz(1.1327789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-3.1277593) q[0];
sx q[0];
rz(-0.56202373) q[0];
sx q[0];
rz(-3.1264937) q[0];
rz(-0.12241441) q[1];
sx q[1];
rz(-1.2950803) q[1];
sx q[1];
rz(-2.1133568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5153582) q[0];
sx q[0];
rz(-1.6263824) q[0];
sx q[0];
rz(2.4856011) q[0];
rz(-pi) q[1];
rz(-1.2790658) q[2];
sx q[2];
rz(-2.5030067) q[2];
sx q[2];
rz(1.2995468) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15785698) q[1];
sx q[1];
rz(-2.7089556) q[1];
sx q[1];
rz(-3.0745201) q[1];
rz(-pi) q[2];
rz(2.8260303) q[3];
sx q[3];
rz(-1.3129819) q[3];
sx q[3];
rz(0.97168789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5369109) q[2];
sx q[2];
rz(-1.2102419) q[2];
sx q[2];
rz(-0.49986419) q[2];
rz(-2.8258421) q[3];
sx q[3];
rz(-2.4707268) q[3];
sx q[3];
rz(-2.1226814) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4278118) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(-2.1977303) q[0];
rz(1.4048514) q[1];
sx q[1];
rz(-1.2662788) q[1];
sx q[1];
rz(-0.92165438) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2093344) q[0];
sx q[0];
rz(-0.32296041) q[0];
sx q[0];
rz(-2.9684116) q[0];
rz(-pi) q[1];
rz(1.2723075) q[2];
sx q[2];
rz(-1.8137852) q[2];
sx q[2];
rz(1.9858433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5068622) q[1];
sx q[1];
rz(-2.1980739) q[1];
sx q[1];
rz(2.2192713) q[1];
x q[2];
rz(1.5632024) q[3];
sx q[3];
rz(-2.6287327) q[3];
sx q[3];
rz(-0.22554413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1450119) q[2];
sx q[2];
rz(-1.2715481) q[2];
sx q[2];
rz(-1.5550295) q[2];
rz(1.0541213) q[3];
sx q[3];
rz(-1.8127706) q[3];
sx q[3];
rz(1.740295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15040511) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(0.64252585) q[0];
rz(1.7036899) q[1];
sx q[1];
rz(-2.3846886) q[1];
sx q[1];
rz(-0.14370758) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3966345) q[0];
sx q[0];
rz(-1.8734583) q[0];
sx q[0];
rz(-0.58505262) q[0];
rz(0.92387284) q[2];
sx q[2];
rz(-1.1387741) q[2];
sx q[2];
rz(-2.1763755) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9160812) q[1];
sx q[1];
rz(-2.4298432) q[1];
sx q[1];
rz(-2.7307778) q[1];
rz(-pi) q[2];
rz(2.4293184) q[3];
sx q[3];
rz(-2.5872719) q[3];
sx q[3];
rz(-1.2620776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6045427) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(-0.26270467) q[2];
rz(-0.25775868) q[3];
sx q[3];
rz(-0.38193211) q[3];
sx q[3];
rz(-2.2885382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2951374) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(-1.2204131) q[0];
rz(-0.48183164) q[1];
sx q[1];
rz(-0.82492963) q[1];
sx q[1];
rz(-2.2442472) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8751971) q[0];
sx q[0];
rz(-2.0854073) q[0];
sx q[0];
rz(1.4304016) q[0];
rz(-pi) q[1];
rz(2.8209575) q[2];
sx q[2];
rz(-0.33679397) q[2];
sx q[2];
rz(1.3367261) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.022696115) q[1];
sx q[1];
rz(-1.5210946) q[1];
sx q[1];
rz(-1.8291874) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7190211) q[3];
sx q[3];
rz(-2.5276042) q[3];
sx q[3];
rz(2.2999291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4974978) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(-3.0268055) q[2];
rz(0.74350205) q[3];
sx q[3];
rz(-1.8758352) q[3];
sx q[3];
rz(2.6928597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262309) q[0];
sx q[0];
rz(-0.76114934) q[0];
sx q[0];
rz(2.1834955) q[0];
rz(-1.505898) q[1];
sx q[1];
rz(-1.1264569) q[1];
sx q[1];
rz(1.1503848) q[1];
rz(2.3957797) q[2];
sx q[2];
rz(-2.0796761) q[2];
sx q[2];
rz(2.833119) q[2];
rz(0.43396797) q[3];
sx q[3];
rz(-1.7229736) q[3];
sx q[3];
rz(2.4255502) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
