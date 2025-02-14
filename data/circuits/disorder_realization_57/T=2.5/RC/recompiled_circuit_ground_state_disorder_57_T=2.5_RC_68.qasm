OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8785414) q[0];
sx q[0];
rz(-0.6113373) q[0];
sx q[0];
rz(1.4933458) q[0];
rz(-2.2892294) q[1];
sx q[1];
rz(-1.747793) q[1];
sx q[1];
rz(1.2300904) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8297079) q[0];
sx q[0];
rz(-2.4915016) q[0];
sx q[0];
rz(2.4310396) q[0];
rz(-pi) q[1];
rz(2.4357721) q[2];
sx q[2];
rz(-0.56669826) q[2];
sx q[2];
rz(1.3659878) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38065698) q[1];
sx q[1];
rz(-1.120693) q[1];
sx q[1];
rz(0.37827079) q[1];
rz(1.854585) q[3];
sx q[3];
rz(-1.6403682) q[3];
sx q[3];
rz(2.8049935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40296951) q[2];
sx q[2];
rz(-1.3961184) q[2];
sx q[2];
rz(-1.9782861) q[2];
rz(-2.2033384) q[3];
sx q[3];
rz(-2.3353751) q[3];
sx q[3];
rz(1.7136542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6152978) q[0];
sx q[0];
rz(-1.2328923) q[0];
sx q[0];
rz(-1.1616608) q[0];
rz(-1.42234) q[1];
sx q[1];
rz(-1.6484345) q[1];
sx q[1];
rz(0.086440451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.739839) q[0];
sx q[0];
rz(-2.5883543) q[0];
sx q[0];
rz(1.2673926) q[0];
rz(-pi) q[1];
rz(-2.6506054) q[2];
sx q[2];
rz(-1.7684968) q[2];
sx q[2];
rz(-0.035015496) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.045043) q[1];
sx q[1];
rz(-3.109906) q[1];
sx q[1];
rz(-2.6331054) q[1];
x q[2];
rz(-3.0753022) q[3];
sx q[3];
rz(-1.4) q[3];
sx q[3];
rz(1.2571245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5689759) q[2];
sx q[2];
rz(-0.47875753) q[2];
sx q[2];
rz(-0.044744603) q[2];
rz(-2.4217126) q[3];
sx q[3];
rz(-1.5174815) q[3];
sx q[3];
rz(-1.450479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2429263) q[0];
sx q[0];
rz(-1.1761605) q[0];
sx q[0];
rz(-1.9157008) q[0];
rz(-2.2236845) q[1];
sx q[1];
rz(-1.8501015) q[1];
sx q[1];
rz(-1.8570522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63375134) q[0];
sx q[0];
rz(-1.0695262) q[0];
sx q[0];
rz(-3.0847163) q[0];
rz(-pi) q[1];
rz(1.578737) q[2];
sx q[2];
rz(-0.94915945) q[2];
sx q[2];
rz(-1.4162743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94938147) q[1];
sx q[1];
rz(-0.95451372) q[1];
sx q[1];
rz(0.90009113) q[1];
rz(-pi) q[2];
rz(-2.7906832) q[3];
sx q[3];
rz(-2.4882462) q[3];
sx q[3];
rz(-2.2782183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.43117493) q[2];
sx q[2];
rz(-1.147889) q[2];
sx q[2];
rz(-2.5244024) q[2];
rz(2.2271473) q[3];
sx q[3];
rz(-0.72943288) q[3];
sx q[3];
rz(-2.7845553) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0964088) q[0];
sx q[0];
rz(-3.0953188) q[0];
sx q[0];
rz(0.98840493) q[0];
rz(1.4811966) q[1];
sx q[1];
rz(-0.55440569) q[1];
sx q[1];
rz(-0.60715094) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84548002) q[0];
sx q[0];
rz(-2.76351) q[0];
sx q[0];
rz(-0.81645672) q[0];
x q[1];
rz(-1.1902624) q[2];
sx q[2];
rz(-0.96559994) q[2];
sx q[2];
rz(2.2303892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97990943) q[1];
sx q[1];
rz(-1.9075127) q[1];
sx q[1];
rz(1.0830131) q[1];
x q[2];
rz(2.061421) q[3];
sx q[3];
rz(-0.20397025) q[3];
sx q[3];
rz(-2.5734165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0561169) q[2];
sx q[2];
rz(-1.3885219) q[2];
sx q[2];
rz(2.910854) q[2];
rz(2.3750316) q[3];
sx q[3];
rz(-2.9975588) q[3];
sx q[3];
rz(1.8421596) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6799927) q[0];
sx q[0];
rz(-1.6277286) q[0];
sx q[0];
rz(-3.0738714) q[0];
rz(1.8374775) q[1];
sx q[1];
rz(-1.291052) q[1];
sx q[1];
rz(-1.4363323) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67187655) q[0];
sx q[0];
rz(-0.8133406) q[0];
sx q[0];
rz(-2.1123625) q[0];
x q[1];
rz(1.3678694) q[2];
sx q[2];
rz(-0.5022011) q[2];
sx q[2];
rz(1.3853488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1166628) q[1];
sx q[1];
rz(-0.58758333) q[1];
sx q[1];
rz(-0.95984103) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1077638) q[3];
sx q[3];
rz(-2.7269533) q[3];
sx q[3];
rz(-1.6920167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0453673) q[2];
sx q[2];
rz(-2.2787091) q[2];
sx q[2];
rz(0.53675845) q[2];
rz(1.2196352) q[3];
sx q[3];
rz(-2.2370971) q[3];
sx q[3];
rz(1.1902683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9959975) q[0];
sx q[0];
rz(-2.1639316) q[0];
sx q[0];
rz(-1.4033432) q[0];
rz(0.504269) q[1];
sx q[1];
rz(-1.206617) q[1];
sx q[1];
rz(2.9020342) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5546416) q[0];
sx q[0];
rz(-1.3901781) q[0];
sx q[0];
rz(2.9819072) q[0];
rz(1.0198258) q[2];
sx q[2];
rz(-1.9293377) q[2];
sx q[2];
rz(-1.3761686) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36167046) q[1];
sx q[1];
rz(-2.4677271) q[1];
sx q[1];
rz(-0.14072953) q[1];
x q[2];
rz(-2.2277545) q[3];
sx q[3];
rz(-2.5322134) q[3];
sx q[3];
rz(2.1003124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2704894) q[2];
sx q[2];
rz(-1.6051382) q[2];
sx q[2];
rz(-0.25320539) q[2];
rz(-2.1829055) q[3];
sx q[3];
rz(-0.91512338) q[3];
sx q[3];
rz(-1.8668713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36509982) q[0];
sx q[0];
rz(-2.4204142) q[0];
sx q[0];
rz(-1.391063) q[0];
rz(-2.5513388) q[1];
sx q[1];
rz(-1.2140467) q[1];
sx q[1];
rz(-2.6388157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25994008) q[0];
sx q[0];
rz(-1.5852196) q[0];
sx q[0];
rz(-1.676286) q[0];
x q[1];
rz(1.1994792) q[2];
sx q[2];
rz(-2.1938112) q[2];
sx q[2];
rz(-0.83500221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16044438) q[1];
sx q[1];
rz(-0.85373053) q[1];
sx q[1];
rz(-2.8129548) q[1];
rz(-pi) q[2];
rz(-0.49248691) q[3];
sx q[3];
rz(-0.60103184) q[3];
sx q[3];
rz(-0.81274569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34660029) q[2];
sx q[2];
rz(-1.184) q[2];
sx q[2];
rz(2.6896175) q[2];
rz(0.31451264) q[3];
sx q[3];
rz(-1.9032685) q[3];
sx q[3];
rz(-3.0934635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5101584) q[0];
sx q[0];
rz(-0.87961125) q[0];
sx q[0];
rz(-1.6360224) q[0];
rz(-3.0360041) q[1];
sx q[1];
rz(-1.0786723) q[1];
sx q[1];
rz(-0.67965913) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1733579) q[0];
sx q[0];
rz(-0.86167406) q[0];
sx q[0];
rz(2.3308332) q[0];
rz(0.95288743) q[2];
sx q[2];
rz(-1.2373239) q[2];
sx q[2];
rz(3.0798387) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9951156) q[1];
sx q[1];
rz(-1.5676284) q[1];
sx q[1];
rz(-1.100844) q[1];
rz(2.258753) q[3];
sx q[3];
rz(-2.853009) q[3];
sx q[3];
rz(1.2478873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6203561) q[2];
sx q[2];
rz(-0.53400365) q[2];
sx q[2];
rz(2.3109069) q[2];
rz(-1.0384809) q[3];
sx q[3];
rz(-2.4887648) q[3];
sx q[3];
rz(1.7415107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-2.0251544) q[0];
sx q[0];
rz(-1.9898299) q[0];
sx q[0];
rz(-2.108216) q[0];
rz(-2.0694464) q[1];
sx q[1];
rz(-1.9629581) q[1];
sx q[1];
rz(0.56195608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011312927) q[0];
sx q[0];
rz(-2.2743158) q[0];
sx q[0];
rz(-0.18248002) q[0];
rz(-pi) q[1];
rz(-1.1216036) q[2];
sx q[2];
rz(-0.70708924) q[2];
sx q[2];
rz(2.0769151) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.55511715) q[1];
sx q[1];
rz(-2.3080864) q[1];
sx q[1];
rz(-1.1039102) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4212332) q[3];
sx q[3];
rz(-3.0834167) q[3];
sx q[3];
rz(0.10567372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.098103913) q[2];
sx q[2];
rz(-1.6489112) q[2];
sx q[2];
rz(1.4081504) q[2];
rz(-3.0812541) q[3];
sx q[3];
rz(-1.8294168) q[3];
sx q[3];
rz(-1.2618802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14335808) q[0];
sx q[0];
rz(-2.3682605) q[0];
sx q[0];
rz(1.7711357) q[0];
rz(2.1342318) q[1];
sx q[1];
rz(-1.7528563) q[1];
sx q[1];
rz(0.14152424) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816982) q[0];
sx q[0];
rz(-2.0260677) q[0];
sx q[0];
rz(-2.0596402) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3958489) q[2];
sx q[2];
rz(-0.088938449) q[2];
sx q[2];
rz(2.3740785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.05647612) q[1];
sx q[1];
rz(-2.543138) q[1];
sx q[1];
rz(2.2920777) q[1];
x q[2];
rz(0.022079682) q[3];
sx q[3];
rz(-1.0754657) q[3];
sx q[3];
rz(-1.2091523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0274028) q[2];
sx q[2];
rz(-1.1833311) q[2];
sx q[2];
rz(1.2664504) q[2];
rz(0.13628422) q[3];
sx q[3];
rz(-2.2253021) q[3];
sx q[3];
rz(2.9436679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4895353) q[0];
sx q[0];
rz(-1.0132402) q[0];
sx q[0];
rz(-1.5860438) q[0];
rz(-0.96792211) q[1];
sx q[1];
rz(-0.5236917) q[1];
sx q[1];
rz(2.0655469) q[1];
rz(0.10836149) q[2];
sx q[2];
rz(-2.32794) q[2];
sx q[2];
rz(1.9447504) q[2];
rz(0.051911671) q[3];
sx q[3];
rz(-1.3811227) q[3];
sx q[3];
rz(1.5863252) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
