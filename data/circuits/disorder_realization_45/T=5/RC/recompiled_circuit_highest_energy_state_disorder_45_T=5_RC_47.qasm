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
rz(1.6765321) q[0];
sx q[0];
rz(-2.9002011) q[0];
sx q[0];
rz(-0.13225947) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(-0.71813923) q[1];
sx q[1];
rz(3.1225966) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26532224) q[0];
sx q[0];
rz(-0.92060584) q[0];
sx q[0];
rz(2.8702535) q[0];
x q[1];
rz(-1.7573518) q[2];
sx q[2];
rz(-1.6636408) q[2];
sx q[2];
rz(-1.873119) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2680696) q[1];
sx q[1];
rz(-1.8933081) q[1];
sx q[1];
rz(1.3871865) q[1];
rz(-pi) q[2];
rz(2.1601342) q[3];
sx q[3];
rz(-0.27772003) q[3];
sx q[3];
rz(3.115417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.12048177) q[2];
sx q[2];
rz(-2.3349473) q[2];
sx q[2];
rz(-2.3497904) q[2];
rz(1.0857438) q[3];
sx q[3];
rz(-1.2191685) q[3];
sx q[3];
rz(2.048548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2708112) q[0];
sx q[0];
rz(-0.15853156) q[0];
sx q[0];
rz(-2.812401) q[0];
rz(1.5022494) q[1];
sx q[1];
rz(-2.2396125) q[1];
sx q[1];
rz(-2.7194068) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3361583) q[0];
sx q[0];
rz(-1.0350739) q[0];
sx q[0];
rz(-1.0725007) q[0];
rz(-pi) q[1];
rz(-2.8453526) q[2];
sx q[2];
rz(-2.5736817) q[2];
sx q[2];
rz(-0.30496773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77779639) q[1];
sx q[1];
rz(-1.2470804) q[1];
sx q[1];
rz(-1.6196851) q[1];
rz(1.4677736) q[3];
sx q[3];
rz(-2.3405082) q[3];
sx q[3];
rz(-2.5689606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72055703) q[2];
sx q[2];
rz(-1.1853508) q[2];
sx q[2];
rz(-1.1400918) q[2];
rz(3.1371878) q[3];
sx q[3];
rz(-1.5604228) q[3];
sx q[3];
rz(0.13979039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36188257) q[0];
sx q[0];
rz(-0.10040586) q[0];
sx q[0];
rz(0.34565872) q[0];
rz(1.0935621) q[1];
sx q[1];
rz(-0.82229096) q[1];
sx q[1];
rz(2.0872769) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39798007) q[0];
sx q[0];
rz(-2.6184887) q[0];
sx q[0];
rz(2.4260957) q[0];
rz(-0.81297154) q[2];
sx q[2];
rz(-1.4958428) q[2];
sx q[2];
rz(-0.048150657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9011827) q[1];
sx q[1];
rz(-2.6844271) q[1];
sx q[1];
rz(2.9967876) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90462138) q[3];
sx q[3];
rz(-1.5890749) q[3];
sx q[3];
rz(-1.3410904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4652319) q[2];
sx q[2];
rz(-0.54764843) q[2];
sx q[2];
rz(2.5993627) q[2];
rz(0.17255653) q[3];
sx q[3];
rz(-1.4901284) q[3];
sx q[3];
rz(1.1057314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3365823) q[0];
sx q[0];
rz(-1.35291) q[0];
sx q[0];
rz(1.8804869) q[0];
rz(-1.1359967) q[1];
sx q[1];
rz(-2.0295862) q[1];
sx q[1];
rz(-0.13066185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4379556) q[0];
sx q[0];
rz(-1.2017631) q[0];
sx q[0];
rz(-0.28965182) q[0];
x q[1];
rz(-0.39009266) q[2];
sx q[2];
rz(-2.60139) q[2];
sx q[2];
rz(-0.15946968) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3522722) q[1];
sx q[1];
rz(-2.6323458) q[1];
sx q[1];
rz(1.6357842) q[1];
rz(-pi) q[2];
rz(-1.2011124) q[3];
sx q[3];
rz(-0.94025984) q[3];
sx q[3];
rz(-0.66555221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7589492) q[2];
sx q[2];
rz(-0.073315695) q[2];
sx q[2];
rz(0.22155133) q[2];
rz(-0.70212901) q[3];
sx q[3];
rz(-1.5789072) q[3];
sx q[3];
rz(-0.73739541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5704983) q[0];
sx q[0];
rz(-1.8966738) q[0];
sx q[0];
rz(-3.0076497) q[0];
rz(-0.36930034) q[1];
sx q[1];
rz(-1.9422653) q[1];
sx q[1];
rz(-1.7514924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5969191) q[0];
sx q[0];
rz(-1.8159165) q[0];
sx q[0];
rz(2.0952203) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8763387) q[2];
sx q[2];
rz(-1.8657639) q[2];
sx q[2];
rz(2.2504928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1410603) q[1];
sx q[1];
rz(-2.5203325) q[1];
sx q[1];
rz(-0.48807524) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69306121) q[3];
sx q[3];
rz(-0.45940889) q[3];
sx q[3];
rz(-1.0074248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.2964581) q[2];
sx q[2];
rz(-0.90201169) q[2];
sx q[2];
rz(1.585539) q[2];
rz(2.8499917) q[3];
sx q[3];
rz(-2.6652938) q[3];
sx q[3];
rz(-0.27879032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45048243) q[0];
sx q[0];
rz(-0.99523681) q[0];
sx q[0];
rz(-2.6158748) q[0];
rz(-0.37824962) q[1];
sx q[1];
rz(-1.6098166) q[1];
sx q[1];
rz(0.27413109) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4843338) q[0];
sx q[0];
rz(-2.5509301) q[0];
sx q[0];
rz(3.0953437) q[0];
rz(2.3304105) q[2];
sx q[2];
rz(-2.5497782) q[2];
sx q[2];
rz(2.9331911) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.200572) q[1];
sx q[1];
rz(-1.1686687) q[1];
sx q[1];
rz(-0.27718039) q[1];
rz(-pi) q[2];
rz(-1.0438221) q[3];
sx q[3];
rz(-0.29114215) q[3];
sx q[3];
rz(3.1028735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5452925) q[2];
sx q[2];
rz(-1.2240852) q[2];
sx q[2];
rz(-2.8833585) q[2];
rz(1.5457414) q[3];
sx q[3];
rz(-1.3629379) q[3];
sx q[3];
rz(2.7365007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3426568) q[0];
sx q[0];
rz(-0.5181784) q[0];
sx q[0];
rz(3.0361191) q[0];
rz(-0.27490973) q[1];
sx q[1];
rz(-1.7173488) q[1];
sx q[1];
rz(0.28604937) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496083) q[0];
sx q[0];
rz(-1.7546856) q[0];
sx q[0];
rz(-1.6894132) q[0];
rz(-2.8477564) q[2];
sx q[2];
rz(-2.2851737) q[2];
sx q[2];
rz(1.9686955) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0487236) q[1];
sx q[1];
rz(-1.9561617) q[1];
sx q[1];
rz(-3.0720622) q[1];
x q[2];
rz(2.4838832) q[3];
sx q[3];
rz(-2.4793787) q[3];
sx q[3];
rz(-0.42218226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1470571) q[2];
sx q[2];
rz(-0.68238443) q[2];
sx q[2];
rz(0.4500173) q[2];
rz(0.58722812) q[3];
sx q[3];
rz(-1.3415895) q[3];
sx q[3];
rz(-1.2500866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8477972) q[0];
sx q[0];
rz(-1.4601409) q[0];
sx q[0];
rz(1.1867123) q[0];
rz(-0.31934357) q[1];
sx q[1];
rz(-1.0217383) q[1];
sx q[1];
rz(-0.26225463) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7569993) q[0];
sx q[0];
rz(-2.5525064) q[0];
sx q[0];
rz(-0.37269816) q[0];
rz(-pi) q[1];
rz(-1.9840711) q[2];
sx q[2];
rz(-1.6954633) q[2];
sx q[2];
rz(-2.0365798) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1848534) q[1];
sx q[1];
rz(-2.112084) q[1];
sx q[1];
rz(3.0459896) q[1];
rz(2.6899509) q[3];
sx q[3];
rz(-1.181895) q[3];
sx q[3];
rz(-1.9967784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.78476) q[2];
sx q[2];
rz(-0.46669745) q[2];
sx q[2];
rz(-2.812815) q[2];
rz(-1.5152991) q[3];
sx q[3];
rz(-1.7833775) q[3];
sx q[3];
rz(-0.40858194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3560155) q[0];
sx q[0];
rz(-2.6308036) q[0];
sx q[0];
rz(-2.862634) q[0];
rz(-0.51697671) q[1];
sx q[1];
rz(-2.7644283) q[1];
sx q[1];
rz(1.4252211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0128929) q[0];
sx q[0];
rz(-2.1919247) q[0];
sx q[0];
rz(1.5170044) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2478156) q[2];
sx q[2];
rz(-2.2502459) q[2];
sx q[2];
rz(-0.51365863) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60630262) q[1];
sx q[1];
rz(-2.6043476) q[1];
sx q[1];
rz(0.37588889) q[1];
rz(-0.29028671) q[3];
sx q[3];
rz(-0.93433524) q[3];
sx q[3];
rz(-2.2009714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4031389) q[2];
sx q[2];
rz(-0.29610115) q[2];
sx q[2];
rz(2.9816755) q[2];
rz(-2.3219409) q[3];
sx q[3];
rz(-1.22217) q[3];
sx q[3];
rz(0.73614365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4495471) q[0];
sx q[0];
rz(-1.291438) q[0];
sx q[0];
rz(-0.29577574) q[0];
rz(-2.6462818) q[1];
sx q[1];
rz(-0.43423978) q[1];
sx q[1];
rz(-1.5089418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2343649) q[0];
sx q[0];
rz(-1.1746658) q[0];
sx q[0];
rz(-0.38568003) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0894439) q[2];
sx q[2];
rz(-1.2189157) q[2];
sx q[2];
rz(-0.032648409) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7612696) q[1];
sx q[1];
rz(-1.4810303) q[1];
sx q[1];
rz(1.5029546) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1811352) q[3];
sx q[3];
rz(-1.1077322) q[3];
sx q[3];
rz(-0.72805007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32790023) q[2];
sx q[2];
rz(-1.4699961) q[2];
sx q[2];
rz(-0.63391614) q[2];
rz(-3.0564803) q[3];
sx q[3];
rz(-1.2462933) q[3];
sx q[3];
rz(2.3688721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2278628) q[0];
sx q[0];
rz(-1.5999404) q[0];
sx q[0];
rz(3.0318442) q[0];
rz(-1.2211424) q[1];
sx q[1];
rz(-0.84538645) q[1];
sx q[1];
rz(0.56199817) q[1];
rz(2.0094677) q[2];
sx q[2];
rz(-1.7747468) q[2];
sx q[2];
rz(-2.940831) q[2];
rz(-1.802465) q[3];
sx q[3];
rz(-1.5944407) q[3];
sx q[3];
rz(-3.1138314) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
