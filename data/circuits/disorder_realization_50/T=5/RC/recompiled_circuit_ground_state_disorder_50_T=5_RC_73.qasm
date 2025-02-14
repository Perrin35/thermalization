OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.1640373) q[0];
sx q[0];
rz(1.2241751) q[0];
sx q[0];
rz(13.847142) q[0];
rz(-0.7723074) q[1];
sx q[1];
rz(-0.63900715) q[1];
sx q[1];
rz(-2.8746936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9123906) q[0];
sx q[0];
rz(-1.4773737) q[0];
sx q[0];
rz(1.6673198) q[0];
x q[1];
rz(-2.2592531) q[2];
sx q[2];
rz(-1.7710476) q[2];
sx q[2];
rz(2.3941674) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7952629) q[1];
sx q[1];
rz(-0.85546934) q[1];
sx q[1];
rz(0.99154179) q[1];
rz(-pi) q[2];
rz(-2.9270615) q[3];
sx q[3];
rz(-0.45154587) q[3];
sx q[3];
rz(2.8928234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.018517) q[2];
sx q[2];
rz(-2.1776431) q[2];
sx q[2];
rz(2.4289971) q[2];
rz(2.9739042) q[3];
sx q[3];
rz(-2.2919787) q[3];
sx q[3];
rz(-0.4445506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1009104) q[0];
sx q[0];
rz(-1.7829144) q[0];
sx q[0];
rz(-1.4605301) q[0];
rz(-1.679861) q[1];
sx q[1];
rz(-2.0748383) q[1];
sx q[1];
rz(1.1163968) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3202964) q[0];
sx q[0];
rz(-1.5460933) q[0];
sx q[0];
rz(-1.6139834) q[0];
rz(-pi) q[1];
rz(-1.8113715) q[2];
sx q[2];
rz(-0.5690388) q[2];
sx q[2];
rz(2.9441058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48076377) q[1];
sx q[1];
rz(-2.6763981) q[1];
sx q[1];
rz(-0.61252266) q[1];
x q[2];
rz(-0.41023572) q[3];
sx q[3];
rz(-1.6950775) q[3];
sx q[3];
rz(-1.4924442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98550335) q[2];
sx q[2];
rz(-2.6458793) q[2];
sx q[2];
rz(-0.25225857) q[2];
rz(1.0559399) q[3];
sx q[3];
rz(-0.85774937) q[3];
sx q[3];
rz(1.0004388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5591739) q[0];
sx q[0];
rz(-0.90455872) q[0];
sx q[0];
rz(-2.3140123) q[0];
rz(-2.6170392) q[1];
sx q[1];
rz(-1.2453715) q[1];
sx q[1];
rz(-0.96447271) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8480769) q[0];
sx q[0];
rz(-1.6381253) q[0];
sx q[0];
rz(0.99253722) q[0];
rz(-1.3521503) q[2];
sx q[2];
rz(-0.3836358) q[2];
sx q[2];
rz(-2.2927257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71389329) q[1];
sx q[1];
rz(-1.5212219) q[1];
sx q[1];
rz(-2.5109641) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12196864) q[3];
sx q[3];
rz(-0.44443529) q[3];
sx q[3];
rz(2.0214863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7437637) q[2];
sx q[2];
rz(-2.9167794) q[2];
sx q[2];
rz(1.942983) q[2];
rz(-1.648929) q[3];
sx q[3];
rz(-2.1677833) q[3];
sx q[3];
rz(-2.5627513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8156133) q[0];
sx q[0];
rz(-2.9878243) q[0];
sx q[0];
rz(1.9258668) q[0];
rz(-0.99023306) q[1];
sx q[1];
rz(-0.74140397) q[1];
sx q[1];
rz(0.56891099) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7051681) q[0];
sx q[0];
rz(-1.2878083) q[0];
sx q[0];
rz(2.1678336) q[0];
rz(-pi) q[1];
rz(2.1808938) q[2];
sx q[2];
rz(-0.44126661) q[2];
sx q[2];
rz(0.91741761) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4601656) q[1];
sx q[1];
rz(-0.54124628) q[1];
sx q[1];
rz(-0.98249225) q[1];
x q[2];
rz(2.3936233) q[3];
sx q[3];
rz(-1.0587721) q[3];
sx q[3];
rz(1.0486163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8625921) q[2];
sx q[2];
rz(-0.96070868) q[2];
sx q[2];
rz(-0.29435364) q[2];
rz(2.477008) q[3];
sx q[3];
rz(-0.75988257) q[3];
sx q[3];
rz(-1.4544646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6890474) q[0];
sx q[0];
rz(-1.4280467) q[0];
sx q[0];
rz(0.31387615) q[0];
rz(-0.9017871) q[1];
sx q[1];
rz(-1.7867463) q[1];
sx q[1];
rz(-1.4656167) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0756158) q[0];
sx q[0];
rz(-0.1217723) q[0];
sx q[0];
rz(-1.3568702) q[0];
rz(-pi) q[1];
rz(-1.8214297) q[2];
sx q[2];
rz(-1.4610372) q[2];
sx q[2];
rz(0.23792371) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61205381) q[1];
sx q[1];
rz(-1.1837237) q[1];
sx q[1];
rz(-1.4899292) q[1];
rz(-0.05081049) q[3];
sx q[3];
rz(-2.5103995) q[3];
sx q[3];
rz(2.7630591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9471211) q[2];
sx q[2];
rz(-2.3919545) q[2];
sx q[2];
rz(-1.0844024) q[2];
rz(-2.8179152) q[3];
sx q[3];
rz(-2.4249228) q[3];
sx q[3];
rz(-0.22423854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1502458) q[0];
sx q[0];
rz(-0.24557376) q[0];
sx q[0];
rz(2.7299951) q[0];
rz(-2.4161074) q[1];
sx q[1];
rz(-1.2440224) q[1];
sx q[1];
rz(1.4010319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65832018) q[0];
sx q[0];
rz(-1.7651157) q[0];
sx q[0];
rz(-0.65388314) q[0];
rz(-0.7404003) q[2];
sx q[2];
rz(-1.6639189) q[2];
sx q[2];
rz(1.653506) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.42565517) q[1];
sx q[1];
rz(-1.2285985) q[1];
sx q[1];
rz(-1.9075431) q[1];
x q[2];
rz(-2.6903166) q[3];
sx q[3];
rz(-1.5022931) q[3];
sx q[3];
rz(-2.9344104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97658849) q[2];
sx q[2];
rz(-0.62825957) q[2];
sx q[2];
rz(1.4264301) q[2];
rz(0.12380883) q[3];
sx q[3];
rz(-2.0721469) q[3];
sx q[3];
rz(-2.6482705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.8641149) q[0];
sx q[0];
rz(-1.9954229) q[0];
sx q[0];
rz(1.3364828) q[0];
rz(-0.9388963) q[1];
sx q[1];
rz(-2.0195596) q[1];
sx q[1];
rz(-0.28405651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1868121) q[0];
sx q[0];
rz(-2.4249888) q[0];
sx q[0];
rz(2.6096662) q[0];
rz(0.94921066) q[2];
sx q[2];
rz(-0.27602592) q[2];
sx q[2];
rz(-0.29601184) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8255684) q[1];
sx q[1];
rz(-0.033016769) q[1];
sx q[1];
rz(1.2570639) q[1];
rz(-pi) q[2];
rz(-1.4316931) q[3];
sx q[3];
rz(-1.3918607) q[3];
sx q[3];
rz(-0.034497189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17860086) q[2];
sx q[2];
rz(-2.4602349) q[2];
sx q[2];
rz(2.2006688) q[2];
rz(0.00087794463) q[3];
sx q[3];
rz(-1.9160756) q[3];
sx q[3];
rz(-3.0555449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0451999) q[0];
sx q[0];
rz(-2.2470076) q[0];
sx q[0];
rz(-0.17564242) q[0];
rz(-1.9215709) q[1];
sx q[1];
rz(-2.2717387) q[1];
sx q[1];
rz(2.2697935) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5593431) q[0];
sx q[0];
rz(-1.7212369) q[0];
sx q[0];
rz(-1.3036672) q[0];
rz(-pi) q[1];
x q[1];
rz(1.201725) q[2];
sx q[2];
rz(-2.3184274) q[2];
sx q[2];
rz(1.5368423) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51556522) q[1];
sx q[1];
rz(-2.5142418) q[1];
sx q[1];
rz(0.27317843) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4967887) q[3];
sx q[3];
rz(-1.0153061) q[3];
sx q[3];
rz(-3.1233226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.72622) q[2];
sx q[2];
rz(-1.8612334) q[2];
sx q[2];
rz(-1.0225164) q[2];
rz(0.38698777) q[3];
sx q[3];
rz(-1.756668) q[3];
sx q[3];
rz(2.3302087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33685327) q[0];
sx q[0];
rz(-1.9751208) q[0];
sx q[0];
rz(2.8785896) q[0];
rz(0.31245843) q[1];
sx q[1];
rz(-2.0461021) q[1];
sx q[1];
rz(-1.9591263) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4561037) q[0];
sx q[0];
rz(-1.0668328) q[0];
sx q[0];
rz(2.9188372) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1698007) q[2];
sx q[2];
rz(-2.3070719) q[2];
sx q[2];
rz(2.4311299) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2559465) q[1];
sx q[1];
rz(-1.3722536) q[1];
sx q[1];
rz(0.37862733) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.617917) q[3];
sx q[3];
rz(-1.0645691) q[3];
sx q[3];
rz(2.5158213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.075228127) q[2];
sx q[2];
rz(-1.790739) q[2];
sx q[2];
rz(-0.24825516) q[2];
rz(0.17812854) q[3];
sx q[3];
rz(-1.0823366) q[3];
sx q[3];
rz(-2.3586912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42661509) q[0];
sx q[0];
rz(-0.39418945) q[0];
sx q[0];
rz(2.07975) q[0];
rz(-1.3007523) q[1];
sx q[1];
rz(-1.6164833) q[1];
sx q[1];
rz(1.1311857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4507732) q[0];
sx q[0];
rz(-1.1191219) q[0];
sx q[0];
rz(2.0913238) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6386191) q[2];
sx q[2];
rz(-1.9679234) q[2];
sx q[2];
rz(-0.056588825) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50000514) q[1];
sx q[1];
rz(-1.4785071) q[1];
sx q[1];
rz(0.39938836) q[1];
x q[2];
rz(0.60148539) q[3];
sx q[3];
rz(-1.7643398) q[3];
sx q[3];
rz(-2.7880965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0608369) q[2];
sx q[2];
rz(-0.26641014) q[2];
sx q[2];
rz(0.60144919) q[2];
rz(1.0414177) q[3];
sx q[3];
rz(-1.6626549) q[3];
sx q[3];
rz(-1.3781594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1489442) q[0];
sx q[0];
rz(-2.5288378) q[0];
sx q[0];
rz(0.80768325) q[0];
rz(3.0638937) q[1];
sx q[1];
rz(-2.6014889) q[1];
sx q[1];
rz(-1.4059975) q[1];
rz(2.7355657) q[2];
sx q[2];
rz(-2.1320504) q[2];
sx q[2];
rz(-0.8036094) q[2];
rz(-2.5004456) q[3];
sx q[3];
rz(-2.2401886) q[3];
sx q[3];
rz(-0.2284579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
