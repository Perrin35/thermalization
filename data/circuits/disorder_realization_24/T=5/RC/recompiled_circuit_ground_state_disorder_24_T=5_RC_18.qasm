OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.71819031) q[0];
sx q[0];
rz(-3.1388404) q[0];
sx q[0];
rz(1.0116853) q[0];
rz(-2.6818795) q[1];
sx q[1];
rz(-1.3016394) q[1];
sx q[1];
rz(-2.969892) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8780706) q[0];
sx q[0];
rz(-1.6524914) q[0];
sx q[0];
rz(0.6058713) q[0];
x q[1];
rz(2.6284211) q[2];
sx q[2];
rz(-0.22226873) q[2];
sx q[2];
rz(1.5813511) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6376892) q[1];
sx q[1];
rz(-1.1451964) q[1];
sx q[1];
rz(-0.63987672) q[1];
rz(0.054084973) q[3];
sx q[3];
rz(-1.4776609) q[3];
sx q[3];
rz(0.86291158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1823938) q[2];
sx q[2];
rz(-2.6540519) q[2];
sx q[2];
rz(1.7215151) q[2];
rz(-1.9567664) q[3];
sx q[3];
rz(-1.607837) q[3];
sx q[3];
rz(2.0737295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083953388) q[0];
sx q[0];
rz(-1.6211442) q[0];
sx q[0];
rz(2.6481096) q[0];
rz(0.77589846) q[1];
sx q[1];
rz(-2.6319365) q[1];
sx q[1];
rz(-0.50481558) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0579266) q[0];
sx q[0];
rz(-2.2424881) q[0];
sx q[0];
rz(2.2855285) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8680598) q[2];
sx q[2];
rz(-1.0631764) q[2];
sx q[2];
rz(1.0647237) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2097834) q[1];
sx q[1];
rz(-1.1714913) q[1];
sx q[1];
rz(-2.4834391) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22498954) q[3];
sx q[3];
rz(-2.8691022) q[3];
sx q[3];
rz(2.837473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33727553) q[2];
sx q[2];
rz(-1.3108459) q[2];
sx q[2];
rz(-0.47067019) q[2];
rz(-0.63052952) q[3];
sx q[3];
rz(-1.0258521) q[3];
sx q[3];
rz(1.4001018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0640963) q[0];
sx q[0];
rz(-0.12501669) q[0];
sx q[0];
rz(-0.65573829) q[0];
rz(-0.31133044) q[1];
sx q[1];
rz(-0.89404023) q[1];
sx q[1];
rz(3.0701367) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8690259) q[0];
sx q[0];
rz(-1.9211968) q[0];
sx q[0];
rz(-0.053215543) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33768968) q[2];
sx q[2];
rz(-1.9200385) q[2];
sx q[2];
rz(1.2472635) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7332257) q[1];
sx q[1];
rz(-2.0587615) q[1];
sx q[1];
rz(-3.0674175) q[1];
rz(-pi) q[2];
rz(-1.3585594) q[3];
sx q[3];
rz(-1.0490514) q[3];
sx q[3];
rz(-2.9468342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2283356) q[2];
sx q[2];
rz(-1.5847289) q[2];
sx q[2];
rz(-0.060001686) q[2];
rz(1.9289198) q[3];
sx q[3];
rz(-0.69827497) q[3];
sx q[3];
rz(0.76446271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54670984) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(2.5643964) q[0];
rz(0.82950854) q[1];
sx q[1];
rz(-2.0894158) q[1];
sx q[1];
rz(0.83782354) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4234377) q[0];
sx q[0];
rz(-2.1912247) q[0];
sx q[0];
rz(1.4761094) q[0];
rz(-0.98060645) q[2];
sx q[2];
rz(-1.5204805) q[2];
sx q[2];
rz(-0.8129763) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29463331) q[1];
sx q[1];
rz(-0.75645743) q[1];
sx q[1];
rz(-2.129617) q[1];
rz(-pi) q[2];
rz(-2.6959723) q[3];
sx q[3];
rz(-0.94404781) q[3];
sx q[3];
rz(-2.8901951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.025658) q[2];
sx q[2];
rz(-1.5237153) q[2];
sx q[2];
rz(1.9177829) q[2];
rz(2.6754248) q[3];
sx q[3];
rz(-0.40188447) q[3];
sx q[3];
rz(-0.60552067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8849628) q[0];
sx q[0];
rz(-1.7054568) q[0];
sx q[0];
rz(-0.10109854) q[0];
rz(-2.7283607) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(1.0879999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0383069) q[0];
sx q[0];
rz(-1.8977471) q[0];
sx q[0];
rz(2.0825544) q[0];
x q[1];
rz(-2.3067683) q[2];
sx q[2];
rz(-1.2918279) q[2];
sx q[2];
rz(-0.5912515) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9718379) q[1];
sx q[1];
rz(-2.7248451) q[1];
sx q[1];
rz(0.93081148) q[1];
rz(-1.5002046) q[3];
sx q[3];
rz(-0.90259457) q[3];
sx q[3];
rz(-1.6833504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92855144) q[2];
sx q[2];
rz(-0.88625675) q[2];
sx q[2];
rz(1.5555596) q[2];
rz(1.0726311) q[3];
sx q[3];
rz(-1.6173247) q[3];
sx q[3];
rz(-3.0090289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9685386) q[0];
sx q[0];
rz(-0.88468495) q[0];
sx q[0];
rz(1.435085) q[0];
rz(-0.75812078) q[1];
sx q[1];
rz(-2.4393612) q[1];
sx q[1];
rz(-3.039956) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0291189) q[0];
sx q[0];
rz(-0.85364193) q[0];
sx q[0];
rz(-2.4957335) q[0];
rz(-pi) q[1];
rz(2.0254214) q[2];
sx q[2];
rz(-1.897942) q[2];
sx q[2];
rz(-2.0899555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58858777) q[1];
sx q[1];
rz(-1.2887452) q[1];
sx q[1];
rz(-1.4600919) q[1];
rz(-pi) q[2];
rz(-0.93433617) q[3];
sx q[3];
rz(-2.015967) q[3];
sx q[3];
rz(-2.1195153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8018084) q[2];
sx q[2];
rz(-2.0399751) q[2];
sx q[2];
rz(1.8088809) q[2];
rz(-1.0821651) q[3];
sx q[3];
rz(-2.3887631) q[3];
sx q[3];
rz(1.4235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69671714) q[0];
sx q[0];
rz(-1.2514021) q[0];
sx q[0];
rz(0.74309293) q[0];
rz(0.7630868) q[1];
sx q[1];
rz(-2.5021195) q[1];
sx q[1];
rz(-0.9185763) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10906405) q[0];
sx q[0];
rz(-1.5458414) q[0];
sx q[0];
rz(-3.1335013) q[0];
x q[1];
rz(-1.0985561) q[2];
sx q[2];
rz(-1.2621321) q[2];
sx q[2];
rz(-2.204208) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1700701) q[1];
sx q[1];
rz(-2.5405209) q[1];
sx q[1];
rz(-0.55723377) q[1];
x q[2];
rz(-2.4174446) q[3];
sx q[3];
rz(-1.6586496) q[3];
sx q[3];
rz(1.4081362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.7243728) q[2];
sx q[2];
rz(-1.959274) q[2];
sx q[2];
rz(0.43583885) q[2];
rz(-3.0271652) q[3];
sx q[3];
rz(-0.2468214) q[3];
sx q[3];
rz(0.59823263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33335394) q[0];
sx q[0];
rz(-2.5840608) q[0];
sx q[0];
rz(-0.58407855) q[0];
rz(-2.7059879) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(1.6216507) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0348957) q[0];
sx q[0];
rz(-2.582552) q[0];
sx q[0];
rz(0.99026545) q[0];
rz(-pi) q[1];
rz(-1.0978537) q[2];
sx q[2];
rz(-2.6885037) q[2];
sx q[2];
rz(0.69437969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3585289) q[1];
sx q[1];
rz(-1.7768246) q[1];
sx q[1];
rz(-1.5611783) q[1];
rz(-pi) q[2];
rz(-0.55840839) q[3];
sx q[3];
rz(-2.8388151) q[3];
sx q[3];
rz(1.4608135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.074177563) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(1.1161067) q[2];
rz(-1.3793147) q[3];
sx q[3];
rz(-1.1032871) q[3];
sx q[3];
rz(-3.0232159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8922358) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(0.78999162) q[0];
rz(3.0780011) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(-0.44073179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1871619) q[0];
sx q[0];
rz(-1.2557898) q[0];
sx q[0];
rz(0.31401547) q[0];
rz(-0.49479167) q[2];
sx q[2];
rz(-1.2107009) q[2];
sx q[2];
rz(2.7419213) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72183441) q[1];
sx q[1];
rz(-1.7350619) q[1];
sx q[1];
rz(-1.0895776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7837672) q[3];
sx q[3];
rz(-0.19366385) q[3];
sx q[3];
rz(1.1676027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.26238394) q[2];
sx q[2];
rz(-0.79664207) q[2];
sx q[2];
rz(2.4764496) q[2];
rz(0.77196676) q[3];
sx q[3];
rz(-2.3699581) q[3];
sx q[3];
rz(1.6316679) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36122286) q[0];
sx q[0];
rz(-0.64798111) q[0];
sx q[0];
rz(-2.1369456) q[0];
rz(1.7338344) q[1];
sx q[1];
rz(-1.6561457) q[1];
sx q[1];
rz(-1.8515324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4125975) q[0];
sx q[0];
rz(-2.3069803) q[0];
sx q[0];
rz(0.33552977) q[0];
x q[1];
rz(0.97083135) q[2];
sx q[2];
rz(-1.7391676) q[2];
sx q[2];
rz(-1.4860934) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2659246) q[1];
sx q[1];
rz(-2.7920715) q[1];
sx q[1];
rz(-2.4825579) q[1];
x q[2];
rz(-0.1741039) q[3];
sx q[3];
rz(-2.1645567) q[3];
sx q[3];
rz(1.6792149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0066234) q[2];
sx q[2];
rz(-1.9776191) q[2];
sx q[2];
rz(0.25035614) q[2];
rz(2.1641459) q[3];
sx q[3];
rz(-1.9340632) q[3];
sx q[3];
rz(-1.4704963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95494315) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(0.74465887) q[1];
sx q[1];
rz(-0.79569334) q[1];
sx q[1];
rz(-2.4462499) q[1];
rz(0.087564115) q[2];
sx q[2];
rz(-0.14391218) q[2];
sx q[2];
rz(-1.1434855) q[2];
rz(2.0674191) q[3];
sx q[3];
rz(-1.9575098) q[3];
sx q[3];
rz(-2.248275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
