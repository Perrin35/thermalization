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
rz(-1.3574358) q[0];
sx q[0];
rz(3.7665851) q[0];
sx q[0];
rz(7.9023043) q[0];
rz(1.1433262) q[1];
sx q[1];
rz(3.7781236) q[1];
sx q[1];
rz(14.349024) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56464897) q[0];
sx q[0];
rz(-1.859382) q[0];
sx q[0];
rz(0.99228575) q[0];
rz(-pi) q[1];
rz(-0.29523103) q[2];
sx q[2];
rz(-1.2157939) q[2];
sx q[2];
rz(-0.4685185) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3431817) q[1];
sx q[1];
rz(-0.34920563) q[1];
sx q[1];
rz(-2.3073372) q[1];
x q[2];
rz(-1.3728849) q[3];
sx q[3];
rz(-2.6555914) q[3];
sx q[3];
rz(2.3305691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4437359) q[2];
sx q[2];
rz(-2.3190658) q[2];
sx q[2];
rz(-2.6653384) q[2];
rz(3.0016628) q[3];
sx q[3];
rz(-1.5014476) q[3];
sx q[3];
rz(-0.073277624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78854617) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(2.7781558) q[0];
rz(-2.4711171) q[1];
sx q[1];
rz(-1.3251708) q[1];
sx q[1];
rz(-1.0796116) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19980656) q[0];
sx q[0];
rz(-0.58871709) q[0];
sx q[0];
rz(-1.1325645) q[0];
x q[1];
rz(-0.13021664) q[2];
sx q[2];
rz(-1.0417582) q[2];
sx q[2];
rz(-2.0646273) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55077165) q[1];
sx q[1];
rz(-1.6147922) q[1];
sx q[1];
rz(1.3293056) q[1];
rz(2.7538746) q[3];
sx q[3];
rz(-2.5250375) q[3];
sx q[3];
rz(-1.316837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4662027) q[2];
sx q[2];
rz(-2.4066996) q[2];
sx q[2];
rz(-2.7030763) q[2];
rz(1.0112666) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(-1.6605759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48280516) q[0];
sx q[0];
rz(-0.52849448) q[0];
sx q[0];
rz(2.3129789) q[0];
rz(-1.8939182) q[1];
sx q[1];
rz(-1.9745461) q[1];
sx q[1];
rz(2.8819328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8369401) q[0];
sx q[0];
rz(-0.90422179) q[0];
sx q[0];
rz(0.34662342) q[0];
rz(-0.76918645) q[2];
sx q[2];
rz(-2.63317) q[2];
sx q[2];
rz(1.891013) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6344172) q[1];
sx q[1];
rz(-2.2440662) q[1];
sx q[1];
rz(0.61051621) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8494495) q[3];
sx q[3];
rz(-2.202987) q[3];
sx q[3];
rz(1.8504253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9344249) q[2];
sx q[2];
rz(-1.5639037) q[2];
sx q[2];
rz(-0.820532) q[2];
rz(-0.82646838) q[3];
sx q[3];
rz(-2.1491437) q[3];
sx q[3];
rz(1.3124527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25927037) q[0];
sx q[0];
rz(-0.81040183) q[0];
sx q[0];
rz(-0.93233863) q[0];
rz(1.8238292) q[1];
sx q[1];
rz(-1.980314) q[1];
sx q[1];
rz(-0.12277776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9230292) q[0];
sx q[0];
rz(-1.6549131) q[0];
sx q[0];
rz(0.53972678) q[0];
rz(2.3244395) q[2];
sx q[2];
rz(-1.3434402) q[2];
sx q[2];
rz(1.0447234) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9726561) q[1];
sx q[1];
rz(-2.3984809) q[1];
sx q[1];
rz(1.572322) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2776385) q[3];
sx q[3];
rz(-1.0256301) q[3];
sx q[3];
rz(2.4923862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5968898) q[2];
sx q[2];
rz(-2.2201316) q[2];
sx q[2];
rz(0.60453647) q[2];
rz(0.70945159) q[3];
sx q[3];
rz(-0.76990288) q[3];
sx q[3];
rz(-2.3179222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9619047) q[0];
sx q[0];
rz(-0.12374319) q[0];
sx q[0];
rz(-1.7748348) q[0];
rz(-0.99864117) q[1];
sx q[1];
rz(-0.81592453) q[1];
sx q[1];
rz(2.6008115) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75434573) q[0];
sx q[0];
rz(-1.7779113) q[0];
sx q[0];
rz(-1.4282754) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3983058) q[2];
sx q[2];
rz(-2.1319445) q[2];
sx q[2];
rz(2.4744792) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1652226) q[1];
sx q[1];
rz(-1.762855) q[1];
sx q[1];
rz(2.8403175) q[1];
x q[2];
rz(-2.861372) q[3];
sx q[3];
rz(-2.5570405) q[3];
sx q[3];
rz(1.1506129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6775386) q[2];
sx q[2];
rz(-0.81391922) q[2];
sx q[2];
rz(0.82826725) q[2];
rz(-1.6895435) q[3];
sx q[3];
rz(-1.6477081) q[3];
sx q[3];
rz(-2.233708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-1.117711) q[0];
sx q[0];
rz(-1.1957059) q[0];
sx q[0];
rz(-1.9697795) q[0];
rz(-0.68012971) q[1];
sx q[1];
rz(-1.9155733) q[1];
sx q[1];
rz(2.5853058) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85324641) q[0];
sx q[0];
rz(-0.46015209) q[0];
sx q[0];
rz(0.96921885) q[0];
rz(-pi) q[1];
x q[1];
rz(1.589619) q[2];
sx q[2];
rz(-0.14655098) q[2];
sx q[2];
rz(2.3702247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.43832512) q[1];
sx q[1];
rz(-0.57619737) q[1];
sx q[1];
rz(3.0238815) q[1];
x q[2];
rz(2.5056434) q[3];
sx q[3];
rz(-1.9665641) q[3];
sx q[3];
rz(1.1102939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3351626) q[2];
sx q[2];
rz(-1.8589636) q[2];
sx q[2];
rz(-0.20509091) q[2];
rz(-0.80275503) q[3];
sx q[3];
rz(-2.0358678) q[3];
sx q[3];
rz(0.55454379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71451181) q[0];
sx q[0];
rz(-2.0682122) q[0];
sx q[0];
rz(0.22739534) q[0];
rz(2.3511476) q[1];
sx q[1];
rz(-2.3643654) q[1];
sx q[1];
rz(2.3048185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28115434) q[0];
sx q[0];
rz(-1.7895164) q[0];
sx q[0];
rz(2.6118096) q[0];
x q[1];
rz(-0.58800943) q[2];
sx q[2];
rz(-1.8711897) q[2];
sx q[2];
rz(-2.755185) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6001125) q[1];
sx q[1];
rz(-1.2573693) q[1];
sx q[1];
rz(-2.7164943) q[1];
rz(-pi) q[2];
rz(-2.3804139) q[3];
sx q[3];
rz(-1.8463148) q[3];
sx q[3];
rz(2.8305284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5493912) q[2];
sx q[2];
rz(-1.4791919) q[2];
sx q[2];
rz(0.57360348) q[2];
rz(0.32522374) q[3];
sx q[3];
rz(-0.89541382) q[3];
sx q[3];
rz(-0.4044683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0322872) q[0];
sx q[0];
rz(-2.9635297) q[0];
sx q[0];
rz(-0.67894116) q[0];
rz(-3.0067054) q[1];
sx q[1];
rz(-1.0604246) q[1];
sx q[1];
rz(1.7178242) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7622691) q[0];
sx q[0];
rz(-2.0998619) q[0];
sx q[0];
rz(-2.6075415) q[0];
x q[1];
rz(-2.1982694) q[2];
sx q[2];
rz(-0.8119785) q[2];
sx q[2];
rz(1.8705778) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9805124) q[1];
sx q[1];
rz(-2.3337165) q[1];
sx q[1];
rz(-3.1001311) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0342977) q[3];
sx q[3];
rz(-2.6141659) q[3];
sx q[3];
rz(-0.81125439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5941102) q[2];
sx q[2];
rz(-2.2308733) q[2];
sx q[2];
rz(-0.2317079) q[2];
rz(2.1314651) q[3];
sx q[3];
rz(-1.908327) q[3];
sx q[3];
rz(-2.2505984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6758839) q[0];
sx q[0];
rz(-2.3864585) q[0];
sx q[0];
rz(2.6278507) q[0];
rz(-1.4328009) q[1];
sx q[1];
rz(-1.865973) q[1];
sx q[1];
rz(-1.5076216) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3604853) q[0];
sx q[0];
rz(-1.6193266) q[0];
sx q[0];
rz(-0.13412181) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7319674) q[2];
sx q[2];
rz(-1.320376) q[2];
sx q[2];
rz(-2.7572244) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3097174) q[1];
sx q[1];
rz(-0.17942218) q[1];
sx q[1];
rz(2.37876) q[1];
x q[2];
rz(-1.6674897) q[3];
sx q[3];
rz(-2.8887199) q[3];
sx q[3];
rz(-2.0677572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0214009) q[2];
sx q[2];
rz(-0.84453619) q[2];
sx q[2];
rz(2.8821442) q[2];
rz(-1.1546968) q[3];
sx q[3];
rz(-2.3754933) q[3];
sx q[3];
rz(1.8416789) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4907289) q[0];
sx q[0];
rz(-1.4800973) q[0];
sx q[0];
rz(1.4890626) q[0];
rz(-0.64000714) q[1];
sx q[1];
rz(-2.044544) q[1];
sx q[1];
rz(2.1515813) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8973577) q[0];
sx q[0];
rz(-1.4264297) q[0];
sx q[0];
rz(-1.6598527) q[0];
rz(2.0406538) q[2];
sx q[2];
rz(-0.85950101) q[2];
sx q[2];
rz(-2.2180722) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4943018) q[1];
sx q[1];
rz(-1.1539946) q[1];
sx q[1];
rz(0.17636645) q[1];
rz(-pi) q[2];
rz(1.3481989) q[3];
sx q[3];
rz(-1.5951459) q[3];
sx q[3];
rz(2.3241732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2085569) q[2];
sx q[2];
rz(-1.8427589) q[2];
sx q[2];
rz(-2.9162858) q[2];
rz(-2.2202668) q[3];
sx q[3];
rz(-0.39042979) q[3];
sx q[3];
rz(0.050617378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0113572) q[0];
sx q[0];
rz(-2.7738032) q[0];
sx q[0];
rz(-2.2875447) q[0];
rz(-1.8232952) q[1];
sx q[1];
rz(-1.5250991) q[1];
sx q[1];
rz(-0.93223882) q[1];
rz(3.0443424) q[2];
sx q[2];
rz(-1.3122059) q[2];
sx q[2];
rz(-2.9065804) q[2];
rz(1.8665707) q[3];
sx q[3];
rz(-0.56850962) q[3];
sx q[3];
rz(-3.1103984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
