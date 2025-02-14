OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.652777) q[0];
sx q[0];
rz(-2.1762159) q[0];
sx q[0];
rz(0.91685796) q[0];
rz(-1.7011473) q[1];
sx q[1];
rz(-2.1135766) q[1];
sx q[1];
rz(-2.4457959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0442565) q[0];
sx q[0];
rz(-1.5166266) q[0];
sx q[0];
rz(-1.143881) q[0];
x q[1];
rz(-1.5231499) q[2];
sx q[2];
rz(-1.5334763) q[2];
sx q[2];
rz(-0.19633987) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5196788) q[1];
sx q[1];
rz(-1.7413472) q[1];
sx q[1];
rz(-0.38003731) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0654109) q[3];
sx q[3];
rz(-1.8997703) q[3];
sx q[3];
rz(1.9291725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9838788) q[2];
sx q[2];
rz(-0.82008755) q[2];
sx q[2];
rz(-2.2411818) q[2];
rz(2.3484717) q[3];
sx q[3];
rz(-1.6499062) q[3];
sx q[3];
rz(2.4429863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124436) q[0];
sx q[0];
rz(-1.1684893) q[0];
sx q[0];
rz(-0.68552619) q[0];
rz(-2.7087052) q[1];
sx q[1];
rz(-1.8750178) q[1];
sx q[1];
rz(-1.8276851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066931574) q[0];
sx q[0];
rz(-1.0767796) q[0];
sx q[0];
rz(0.95290174) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0643343) q[2];
sx q[2];
rz(-0.80776513) q[2];
sx q[2];
rz(-1.1459717) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6095143) q[1];
sx q[1];
rz(-1.1066965) q[1];
sx q[1];
rz(1.0232693) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8223313) q[3];
sx q[3];
rz(-2.1135181) q[3];
sx q[3];
rz(-2.4595878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1757397) q[2];
sx q[2];
rz(-0.41030914) q[2];
sx q[2];
rz(-0.24743323) q[2];
rz(3.0257709) q[3];
sx q[3];
rz(-1.7146866) q[3];
sx q[3];
rz(-0.57870948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1568569) q[0];
sx q[0];
rz(-0.37608376) q[0];
sx q[0];
rz(1.0648741) q[0];
rz(0.6330601) q[1];
sx q[1];
rz(-1.9828911) q[1];
sx q[1];
rz(-0.13654581) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8234607) q[0];
sx q[0];
rz(-0.7366418) q[0];
sx q[0];
rz(0.036286906) q[0];
rz(2.0338796) q[2];
sx q[2];
rz(-1.255724) q[2];
sx q[2];
rz(-1.3797494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28440839) q[1];
sx q[1];
rz(-2.2940002) q[1];
sx q[1];
rz(-2.492849) q[1];
rz(-pi) q[2];
rz(2.7259215) q[3];
sx q[3];
rz(-1.9686254) q[3];
sx q[3];
rz(1.3597387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69918767) q[2];
sx q[2];
rz(-1.2523315) q[2];
sx q[2];
rz(-1.5415972) q[2];
rz(1.0934746) q[3];
sx q[3];
rz(-1.5019865) q[3];
sx q[3];
rz(0.34539616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.71444756) q[0];
sx q[0];
rz(-0.26145014) q[0];
sx q[0];
rz(2.5647822) q[0];
rz(-2.414074) q[1];
sx q[1];
rz(-2.0501761) q[1];
sx q[1];
rz(-2.2732546) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5142347) q[0];
sx q[0];
rz(-1.9666202) q[0];
sx q[0];
rz(-2.4577599) q[0];
x q[1];
rz(-0.98855726) q[2];
sx q[2];
rz(-1.120627) q[2];
sx q[2];
rz(2.6188713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7817745) q[1];
sx q[1];
rz(-1.9137772) q[1];
sx q[1];
rz(-2.4510259) q[1];
x q[2];
rz(2.1959325) q[3];
sx q[3];
rz(-1.1468317) q[3];
sx q[3];
rz(1.6062007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19520983) q[2];
sx q[2];
rz(-1.752744) q[2];
sx q[2];
rz(0.69551224) q[2];
rz(2.9269311) q[3];
sx q[3];
rz(-1.6404459) q[3];
sx q[3];
rz(-2.3282839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5401841) q[0];
sx q[0];
rz(-2.5267127) q[0];
sx q[0];
rz(2.0414798) q[0];
rz(-2.9310138) q[1];
sx q[1];
rz(-0.63340488) q[1];
sx q[1];
rz(-1.6483773) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8029047) q[0];
sx q[0];
rz(-2.4409823) q[0];
sx q[0];
rz(2.3466097) q[0];
rz(-1.3927685) q[2];
sx q[2];
rz(-1.653394) q[2];
sx q[2];
rz(-0.96197739) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1112342) q[1];
sx q[1];
rz(-2.6270297) q[1];
sx q[1];
rz(1.8909992) q[1];
rz(-1.8214956) q[3];
sx q[3];
rz(-0.9189824) q[3];
sx q[3];
rz(1.0221611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4139999) q[2];
sx q[2];
rz(-1.4294581) q[2];
sx q[2];
rz(0.89699888) q[2];
rz(1.8619079) q[3];
sx q[3];
rz(-1.0106267) q[3];
sx q[3];
rz(0.047209386) q[3];
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
rz(1.9996662) q[0];
sx q[0];
rz(-2.5157295) q[0];
sx q[0];
rz(2.8114787) q[0];
rz(0.64078981) q[1];
sx q[1];
rz(-1.375744) q[1];
sx q[1];
rz(-2.1660588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27811089) q[0];
sx q[0];
rz(-1.6275102) q[0];
sx q[0];
rz(-0.5910763) q[0];
rz(-0.87557809) q[2];
sx q[2];
rz(-0.40603128) q[2];
sx q[2];
rz(-0.46689597) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2745191) q[1];
sx q[1];
rz(-2.0049067) q[1];
sx q[1];
rz(-1.8810424) q[1];
rz(-pi) q[2];
rz(1.9721071) q[3];
sx q[3];
rz(-2.6875477) q[3];
sx q[3];
rz(1.3523449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2834125) q[2];
sx q[2];
rz(-2.8874669) q[2];
sx q[2];
rz(2.251909) q[2];
rz(-2.5455918) q[3];
sx q[3];
rz(-1.069205) q[3];
sx q[3];
rz(-3.0349351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0058873) q[0];
sx q[0];
rz(-1.1361253) q[0];
sx q[0];
rz(1.1329875) q[0];
rz(2.5021482) q[1];
sx q[1];
rz(-2.0638128) q[1];
sx q[1];
rz(-2.1741672) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9854162) q[0];
sx q[0];
rz(-1.420891) q[0];
sx q[0];
rz(-2.085859) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.899951) q[2];
sx q[2];
rz(-1.2094922) q[2];
sx q[2];
rz(2.6893534) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51951197) q[1];
sx q[1];
rz(-1.4556307) q[1];
sx q[1];
rz(-2.9154816) q[1];
x q[2];
rz(0.65770517) q[3];
sx q[3];
rz(-1.8175565) q[3];
sx q[3];
rz(-1.9384428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7449259) q[2];
sx q[2];
rz(-2.1325839) q[2];
sx q[2];
rz(-0.94433707) q[2];
rz(0.62263292) q[3];
sx q[3];
rz(-2.1532652) q[3];
sx q[3];
rz(-0.78701377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054976376) q[0];
sx q[0];
rz(-0.80626196) q[0];
sx q[0];
rz(0.0042313519) q[0];
rz(-1.6731693) q[1];
sx q[1];
rz(-2.3783042) q[1];
sx q[1];
rz(0.17446336) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28936548) q[0];
sx q[0];
rz(-1.4783637) q[0];
sx q[0];
rz(1.7579745) q[0];
x q[1];
rz(1.4685693) q[2];
sx q[2];
rz(-2.0172628) q[2];
sx q[2];
rz(-0.40418543) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0956438) q[1];
sx q[1];
rz(-0.87365957) q[1];
sx q[1];
rz(0.086177372) q[1];
x q[2];
rz(-3.0765303) q[3];
sx q[3];
rz(-1.8064677) q[3];
sx q[3];
rz(0.5103569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2506432) q[2];
sx q[2];
rz(-3.0477016) q[2];
sx q[2];
rz(0.54035652) q[2];
rz(-3.12449) q[3];
sx q[3];
rz(-0.98186866) q[3];
sx q[3];
rz(-2.4798149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0509725) q[0];
sx q[0];
rz(-2.3391188) q[0];
sx q[0];
rz(-2.354994) q[0];
rz(-2.0358918) q[1];
sx q[1];
rz(-1.4812171) q[1];
sx q[1];
rz(0.22122637) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.008667058) q[0];
sx q[0];
rz(-1.8485539) q[0];
sx q[0];
rz(0.8322258) q[0];
x q[1];
rz(1.7443666) q[2];
sx q[2];
rz(-1.4772869) q[2];
sx q[2];
rz(1.6261404) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3479285) q[1];
sx q[1];
rz(-1.3634487) q[1];
sx q[1];
rz(2.6479032) q[1];
rz(-1.8046326) q[3];
sx q[3];
rz(-0.89409308) q[3];
sx q[3];
rz(1.9399583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1456387) q[2];
sx q[2];
rz(-2.8990539) q[2];
sx q[2];
rz(1.185574) q[2];
rz(-2.0354185) q[3];
sx q[3];
rz(-1.0217051) q[3];
sx q[3];
rz(-2.1574028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36158654) q[0];
sx q[0];
rz(-1.3573283) q[0];
sx q[0];
rz(-0.7097882) q[0];
rz(-0.89372006) q[1];
sx q[1];
rz(-2.5137386) q[1];
sx q[1];
rz(2.6745904) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2115973) q[0];
sx q[0];
rz(-0.30147895) q[0];
sx q[0];
rz(-1.0903574) q[0];
rz(-2.8728629) q[2];
sx q[2];
rz(-2.1788283) q[2];
sx q[2];
rz(0.001645263) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4255387) q[1];
sx q[1];
rz(-1.8093093) q[1];
sx q[1];
rz(-0.038360049) q[1];
x q[2];
rz(-0.67450847) q[3];
sx q[3];
rz(-2.149636) q[3];
sx q[3];
rz(-1.9306435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.799377) q[2];
sx q[2];
rz(-0.90505427) q[2];
sx q[2];
rz(-0.13599914) q[2];
rz(2.811725) q[3];
sx q[3];
rz(-2.0743275) q[3];
sx q[3];
rz(-0.29712591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8316523) q[0];
sx q[0];
rz(-2.0257873) q[0];
sx q[0];
rz(-2.5153487) q[0];
rz(0.88824875) q[1];
sx q[1];
rz(-1.2445969) q[1];
sx q[1];
rz(0.1319763) q[1];
rz(-0.3966847) q[2];
sx q[2];
rz(-1.7783782) q[2];
sx q[2];
rz(1.5967899) q[2];
rz(1.9985401) q[3];
sx q[3];
rz(-0.88870591) q[3];
sx q[3];
rz(0.58042713) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
