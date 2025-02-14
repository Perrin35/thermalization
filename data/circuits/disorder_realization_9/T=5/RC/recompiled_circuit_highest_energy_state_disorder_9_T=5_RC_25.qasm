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
rz(0.63737386) q[0];
sx q[0];
rz(-2.0847991) q[0];
sx q[0];
rz(-0.771653) q[0];
rz(2.9906988) q[1];
sx q[1];
rz(-2.8291191) q[1];
sx q[1];
rz(1.6627275) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61687311) q[0];
sx q[0];
rz(-1.9025073) q[0];
sx q[0];
rz(-1.2613368) q[0];
rz(-2.0885299) q[2];
sx q[2];
rz(-1.7969171) q[2];
sx q[2];
rz(-2.256611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7141908) q[1];
sx q[1];
rz(-1.8334249) q[1];
sx q[1];
rz(2.3841969) q[1];
rz(-pi) q[2];
rz(-0.21680253) q[3];
sx q[3];
rz(-2.4507782) q[3];
sx q[3];
rz(1.3135214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5375157) q[2];
sx q[2];
rz(-2.5274369) q[2];
sx q[2];
rz(2.3094731) q[2];
rz(0.65699792) q[3];
sx q[3];
rz(-1.6712345) q[3];
sx q[3];
rz(2.7895797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81545365) q[0];
sx q[0];
rz(-2.5322545) q[0];
sx q[0];
rz(-2.4916008) q[0];
rz(-3.0187712) q[1];
sx q[1];
rz(-2.717412) q[1];
sx q[1];
rz(-2.4688683) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8433326) q[0];
sx q[0];
rz(-1.6125868) q[0];
sx q[0];
rz(3.1298547) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8287436) q[2];
sx q[2];
rz(-1.8331844) q[2];
sx q[2];
rz(-2.69115) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1993701) q[1];
sx q[1];
rz(-1.7666715) q[1];
sx q[1];
rz(-1.0088527) q[1];
rz(-pi) q[2];
rz(3.1150713) q[3];
sx q[3];
rz(-2.0904198) q[3];
sx q[3];
rz(-0.61557239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0279205) q[2];
sx q[2];
rz(-1.547926) q[2];
sx q[2];
rz(0.38635722) q[2];
rz(1.1682642) q[3];
sx q[3];
rz(-1.2315653) q[3];
sx q[3];
rz(-1.108235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.2439709) q[0];
sx q[0];
rz(-1.6672927) q[0];
sx q[0];
rz(3.0834055) q[0];
rz(-2.8367786) q[1];
sx q[1];
rz(-2.0084281) q[1];
sx q[1];
rz(2.286639) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.917934) q[0];
sx q[0];
rz(-0.068327986) q[0];
sx q[0];
rz(3.1177417) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5714855) q[2];
sx q[2];
rz(-0.7749346) q[2];
sx q[2];
rz(2.3910445) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82423174) q[1];
sx q[1];
rz(-0.79272017) q[1];
sx q[1];
rz(-2.8522245) q[1];
rz(0.26400395) q[3];
sx q[3];
rz(-1.0314634) q[3];
sx q[3];
rz(-1.7738938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.052224) q[2];
sx q[2];
rz(-1.3244018) q[2];
sx q[2];
rz(1.268187) q[2];
rz(-2.7005699) q[3];
sx q[3];
rz(-0.62097582) q[3];
sx q[3];
rz(2.2997901) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0190401) q[0];
sx q[0];
rz(-0.95128107) q[0];
sx q[0];
rz(-2.4265491) q[0];
rz(2.7252281) q[1];
sx q[1];
rz(-2.6519897) q[1];
sx q[1];
rz(0.29410902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.153115) q[0];
sx q[0];
rz(-1.6286958) q[0];
sx q[0];
rz(-2.2672976) q[0];
x q[1];
rz(-1.2671183) q[2];
sx q[2];
rz(-2.9275002) q[2];
sx q[2];
rz(-1.0372727) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70620437) q[1];
sx q[1];
rz(-1.7559373) q[1];
sx q[1];
rz(1.7035083) q[1];
x q[2];
rz(2.5542166) q[3];
sx q[3];
rz(-1.0738292) q[3];
sx q[3];
rz(2.9413829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2405582) q[2];
sx q[2];
rz(-1.3996539) q[2];
sx q[2];
rz(-2.412839) q[2];
rz(2.6567904) q[3];
sx q[3];
rz(-1.0032283) q[3];
sx q[3];
rz(0.17759594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4575551) q[0];
sx q[0];
rz(-1.8192679) q[0];
sx q[0];
rz(2.8628602) q[0];
rz(-0.067525603) q[1];
sx q[1];
rz(-2.4261256) q[1];
sx q[1];
rz(0.50594893) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8794969) q[0];
sx q[0];
rz(-0.56096062) q[0];
sx q[0];
rz(3.0039588) q[0];
rz(-pi) q[1];
rz(-1.672411) q[2];
sx q[2];
rz(-1.7676465) q[2];
sx q[2];
rz(-2.5972899) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7655395) q[1];
sx q[1];
rz(-2.3552822) q[1];
sx q[1];
rz(2.3885352) q[1];
x q[2];
rz(0.26624532) q[3];
sx q[3];
rz(-0.91495401) q[3];
sx q[3];
rz(-0.4474934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94023306) q[2];
sx q[2];
rz(-1.4505922) q[2];
sx q[2];
rz(3.039321) q[2];
rz(2.8376132) q[3];
sx q[3];
rz(-2.961048) q[3];
sx q[3];
rz(-2.6612018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893696) q[0];
sx q[0];
rz(-2.3134573) q[0];
sx q[0];
rz(2.4612259) q[0];
rz(0.96087372) q[1];
sx q[1];
rz(-1.5545132) q[1];
sx q[1];
rz(-2.5804677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8770804) q[0];
sx q[0];
rz(-1.3662369) q[0];
sx q[0];
rz(-1.2007942) q[0];
x q[1];
rz(-1.8041925) q[2];
sx q[2];
rz(-2.4037558) q[2];
sx q[2];
rz(-0.8066752) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.61455124) q[1];
sx q[1];
rz(-2.8379662) q[1];
sx q[1];
rz(2.9580081) q[1];
rz(-2.0797727) q[3];
sx q[3];
rz(-2.1474995) q[3];
sx q[3];
rz(2.1901166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23201021) q[2];
sx q[2];
rz(-0.97753111) q[2];
sx q[2];
rz(1.5634465) q[2];
rz(1.0163418) q[3];
sx q[3];
rz(-1.7137824) q[3];
sx q[3];
rz(1.1004755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72961724) q[0];
sx q[0];
rz(-0.53891861) q[0];
sx q[0];
rz(2.6475661) q[0];
rz(-1.5771075) q[1];
sx q[1];
rz(-0.99948519) q[1];
sx q[1];
rz(0.090242537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7800538) q[0];
sx q[0];
rz(-1.1403196) q[0];
sx q[0];
rz(-1.3132877) q[0];
x q[1];
rz(-2.266062) q[2];
sx q[2];
rz(-1.9771431) q[2];
sx q[2];
rz(-1.4515431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4524432) q[1];
sx q[1];
rz(-1.9724558) q[1];
sx q[1];
rz(2.8287877) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8537515) q[3];
sx q[3];
rz(-1.3850901) q[3];
sx q[3];
rz(2.9880092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51284853) q[2];
sx q[2];
rz(-2.5675842) q[2];
sx q[2];
rz(2.4981456) q[2];
rz(-2.511034) q[3];
sx q[3];
rz(-0.75391155) q[3];
sx q[3];
rz(-3.0923617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1132722) q[0];
sx q[0];
rz(-1.4406818) q[0];
sx q[0];
rz(1.2233618) q[0];
rz(-0.89448294) q[1];
sx q[1];
rz(-2.5136785) q[1];
sx q[1];
rz(-1.9953413) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2558738) q[0];
sx q[0];
rz(-0.71122716) q[0];
sx q[0];
rz(0.92899404) q[0];
rz(-pi) q[1];
rz(-1.1432392) q[2];
sx q[2];
rz(-1.9263066) q[2];
sx q[2];
rz(1.8380788) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1652897) q[1];
sx q[1];
rz(-2.7600651) q[1];
sx q[1];
rz(0.69451992) q[1];
rz(2.8798929) q[3];
sx q[3];
rz(-1.0013805) q[3];
sx q[3];
rz(-0.70139317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.869183) q[2];
sx q[2];
rz(-0.86239186) q[2];
sx q[2];
rz(-1.4746846) q[2];
rz(-2.9278751) q[3];
sx q[3];
rz(-1.4054047) q[3];
sx q[3];
rz(2.2739482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2992582) q[0];
sx q[0];
rz(-2.5304351) q[0];
sx q[0];
rz(0.70511955) q[0];
rz(0.57535386) q[1];
sx q[1];
rz(-1.4082785) q[1];
sx q[1];
rz(0.56952482) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3045159) q[0];
sx q[0];
rz(-1.4936556) q[0];
sx q[0];
rz(1.4894372) q[0];
x q[1];
rz(-0.51409431) q[2];
sx q[2];
rz(-0.55214685) q[2];
sx q[2];
rz(0.29468003) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7895487) q[1];
sx q[1];
rz(-2.1716699) q[1];
sx q[1];
rz(2.4070441) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6105004) q[3];
sx q[3];
rz(-0.55289024) q[3];
sx q[3];
rz(-0.031888511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9465955) q[2];
sx q[2];
rz(-2.1985998) q[2];
sx q[2];
rz(2.1627964) q[2];
rz(-0.43092522) q[3];
sx q[3];
rz(-0.48723358) q[3];
sx q[3];
rz(2.0144958) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8579213) q[0];
sx q[0];
rz(-0.019275276) q[0];
sx q[0];
rz(1.6790144) q[0];
rz(-2.1025533) q[1];
sx q[1];
rz(-1.3835399) q[1];
sx q[1];
rz(1.9336112) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6039654) q[0];
sx q[0];
rz(-0.4930636) q[0];
sx q[0];
rz(-2.2490671) q[0];
x q[1];
rz(-1.6287491) q[2];
sx q[2];
rz(-0.52554146) q[2];
sx q[2];
rz(-1.4144104) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4962311) q[1];
sx q[1];
rz(-1.4252311) q[1];
sx q[1];
rz(2.1880989) q[1];
x q[2];
rz(1.1753756) q[3];
sx q[3];
rz(-1.6334157) q[3];
sx q[3];
rz(2.5015802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.155948) q[2];
sx q[2];
rz(-2.1375103) q[2];
sx q[2];
rz(-1.8211625) q[2];
rz(2.7909347) q[3];
sx q[3];
rz(-2.3242293) q[3];
sx q[3];
rz(-0.58026522) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88631267) q[0];
sx q[0];
rz(-1.9552312) q[0];
sx q[0];
rz(2.2807688) q[0];
rz(-1.6571922) q[1];
sx q[1];
rz(-1.3927554) q[1];
sx q[1];
rz(1.4295255) q[1];
rz(1.6419353) q[2];
sx q[2];
rz(-1.9491458) q[2];
sx q[2];
rz(-0.083935621) q[2];
rz(1.2091533) q[3];
sx q[3];
rz(-0.71930887) q[3];
sx q[3];
rz(0.96439509) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
