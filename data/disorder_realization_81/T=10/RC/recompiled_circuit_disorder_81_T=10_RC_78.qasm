OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(1.3786432) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(0.4508957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9076951) q[0];
sx q[0];
rz(-1.6595708) q[0];
sx q[0];
rz(1.5560454) q[0];
x q[1];
rz(1.0916753) q[2];
sx q[2];
rz(-1.4716822) q[2];
sx q[2];
rz(1.5168158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9860984) q[1];
sx q[1];
rz(-1.7928147) q[1];
sx q[1];
rz(-2.6296031) q[1];
x q[2];
rz(1.4536742) q[3];
sx q[3];
rz(-0.68713596) q[3];
sx q[3];
rz(0.62108921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(-0.84428865) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(-2.8785008) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(-1.1862322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037558) q[0];
sx q[0];
rz(-3.0825966) q[0];
sx q[0];
rz(2.8179413) q[0];
x q[1];
rz(-1.6329174) q[2];
sx q[2];
rz(-1.3717692) q[2];
sx q[2];
rz(1.3028499) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8801404) q[1];
sx q[1];
rz(-1.6958691) q[1];
sx q[1];
rz(0.90420453) q[1];
rz(-pi) q[2];
rz(0.22963345) q[3];
sx q[3];
rz(-0.7080871) q[3];
sx q[3];
rz(-1.1585483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0120323) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(-1.9821232) q[2];
rz(2.7705079) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3804669) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(2.3348715) q[0];
rz(2.9280248) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(-2.321373) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8504282) q[0];
sx q[0];
rz(-0.69520742) q[0];
sx q[0];
rz(-1.3948963) q[0];
rz(1.9633031) q[2];
sx q[2];
rz(-0.5272534) q[2];
sx q[2];
rz(0.96780992) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8169176) q[1];
sx q[1];
rz(-1.1500689) q[1];
sx q[1];
rz(0.38751985) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1180325) q[3];
sx q[3];
rz(-2.0327838) q[3];
sx q[3];
rz(3.1363917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(-2.9860949) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61313066) q[0];
sx q[0];
rz(-2.4202132) q[0];
sx q[0];
rz(0.91127515) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(-1.8211676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62990084) q[0];
sx q[0];
rz(-1.1612079) q[0];
sx q[0];
rz(0.011261777) q[0];
rz(-1.3643866) q[2];
sx q[2];
rz(-2.0136535) q[2];
sx q[2];
rz(2.4368311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28543132) q[1];
sx q[1];
rz(-2.397575) q[1];
sx q[1];
rz(-0.39399778) q[1];
x q[2];
rz(-2.1774125) q[3];
sx q[3];
rz(-2.2607431) q[3];
sx q[3];
rz(-0.53897688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0358255) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(-2.7992115) q[2];
rz(-2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8300366) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(2.2763021) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(-1.8409761) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1814225) q[0];
sx q[0];
rz(-1.2167861) q[0];
sx q[0];
rz(1.5976853) q[0];
x q[1];
rz(-1.6426716) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(-0.14771151) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2758267) q[1];
sx q[1];
rz(-1.2424801) q[1];
sx q[1];
rz(0.57002108) q[1];
x q[2];
rz(-2.8793094) q[3];
sx q[3];
rz(-1.7119006) q[3];
sx q[3];
rz(-1.9067681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0098003) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(-0.47719964) q[2];
rz(2.9495083) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(-2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.79648298) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(-3.1298424) q[0];
rz(-2.5911962) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5884429) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1588622) q[0];
sx q[0];
rz(-1.9059062) q[0];
sx q[0];
rz(1.9431252) q[0];
rz(-1.9366829) q[2];
sx q[2];
rz(-1.5503746) q[2];
sx q[2];
rz(-0.4991971) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.8223871) q[1];
sx q[1];
rz(-2.2777595) q[1];
sx q[1];
rz(1.6824526) q[1];
x q[2];
rz(1.6586967) q[3];
sx q[3];
rz(-0.70471901) q[3];
sx q[3];
rz(2.8040915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6340296) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(-1.7717308) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(2.4643331) q[0];
rz(-0.15180763) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(-0.97704926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517294) q[0];
sx q[0];
rz(-1.6461194) q[0];
sx q[0];
rz(-0.59318869) q[0];
rz(-1.5708959) q[2];
sx q[2];
rz(-1.7028371) q[2];
sx q[2];
rz(-0.11167234) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2077142) q[1];
sx q[1];
rz(-1.9095451) q[1];
sx q[1];
rz(-3.0599041) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7988775) q[3];
sx q[3];
rz(-1.2667155) q[3];
sx q[3];
rz(-2.8976687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7523505) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(1.0127257) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-2.0195122) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(1.8922071) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92266868) q[0];
sx q[0];
rz(-1.7814753) q[0];
sx q[0];
rz(-1.6181437) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47301936) q[2];
sx q[2];
rz(-1.0968536) q[2];
sx q[2];
rz(2.3854286) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74832143) q[1];
sx q[1];
rz(-1.5652579) q[1];
sx q[1];
rz(-3.1194411) q[1];
rz(-pi) q[2];
rz(-1.1864248) q[3];
sx q[3];
rz(-0.66925183) q[3];
sx q[3];
rz(2.3068908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8119048) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(1.9630986) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4998528) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(1.2930124) q[0];
rz(1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(2.5440149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7062089) q[0];
sx q[0];
rz(-1.7058813) q[0];
sx q[0];
rz(0.038890966) q[0];
x q[1];
rz(2.0963247) q[2];
sx q[2];
rz(-1.0768441) q[2];
sx q[2];
rz(-1.8906821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7087047) q[1];
sx q[1];
rz(-1.8678027) q[1];
sx q[1];
rz(-1.2328641) q[1];
rz(-3.0258614) q[3];
sx q[3];
rz(-0.58905187) q[3];
sx q[3];
rz(-0.91050402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22275816) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(-0.90138609) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047886588) q[0];
sx q[0];
rz(-0.77195764) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(0.67217174) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2490847) q[0];
sx q[0];
rz(-1.5605643) q[0];
sx q[0];
rz(-0.0052878629) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1250481) q[2];
sx q[2];
rz(-0.51157727) q[2];
sx q[2];
rz(1.7553154) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.086239554) q[1];
sx q[1];
rz(-1.9085911) q[1];
sx q[1];
rz(-0.71725459) q[1];
rz(-0.43283312) q[3];
sx q[3];
rz(-1.4308617) q[3];
sx q[3];
rz(-2.7452552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6293634) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(2.771634) q[2];
rz(1.6379179) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5794012) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(0.72369408) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(-1.205668) q[2];
sx q[2];
rz(-2.7177313) q[2];
sx q[2];
rz(2.360366) q[2];
rz(-1.1643812) q[3];
sx q[3];
rz(-1.8462528) q[3];
sx q[3];
rz(0.13312199) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];