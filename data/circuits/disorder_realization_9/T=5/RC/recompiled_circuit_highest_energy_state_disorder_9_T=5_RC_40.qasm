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
rz(-0.15089384) q[1];
sx q[1];
rz(-0.3124736) q[1];
sx q[1];
rz(-1.6627275) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61687311) q[0];
sx q[0];
rz(-1.2390854) q[0];
sx q[0];
rz(1.8802558) q[0];
rz(-1.1356699) q[2];
sx q[2];
rz(-0.56081334) q[2];
sx q[2];
rz(2.8307874) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2661765) q[1];
sx q[1];
rz(-2.3485942) q[1];
sx q[1];
rz(-2.7686053) q[1];
rz(-pi) q[2];
rz(2.462384) q[3];
sx q[3];
rz(-1.4333042) q[3];
sx q[3];
rz(0.42543558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5375157) q[2];
sx q[2];
rz(-0.61415577) q[2];
sx q[2];
rz(2.3094731) q[2];
rz(-0.65699792) q[3];
sx q[3];
rz(-1.4703581) q[3];
sx q[3];
rz(-0.35201296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81545365) q[0];
sx q[0];
rz(-2.5322545) q[0];
sx q[0];
rz(-2.4916008) q[0];
rz(-0.12282148) q[1];
sx q[1];
rz(-2.717412) q[1];
sx q[1];
rz(-0.67272433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024351748) q[0];
sx q[0];
rz(-0.043406758) q[0];
sx q[0];
rz(-1.2971334) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31284903) q[2];
sx q[2];
rz(-1.3084083) q[2];
sx q[2];
rz(-2.69115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6482268) q[1];
sx q[1];
rz(-2.1207444) q[1];
sx q[1];
rz(-2.9112766) q[1];
rz(-pi) q[2];
rz(-3.1150713) q[3];
sx q[3];
rz(-1.0511729) q[3];
sx q[3];
rz(-0.61557239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0279205) q[2];
sx q[2];
rz(-1.5936667) q[2];
sx q[2];
rz(-0.38635722) q[2];
rz(1.1682642) q[3];
sx q[3];
rz(-1.2315653) q[3];
sx q[3];
rz(-1.108235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2439709) q[0];
sx q[0];
rz(-1.4742999) q[0];
sx q[0];
rz(3.0834055) q[0];
rz(0.30481401) q[1];
sx q[1];
rz(-2.0084281) q[1];
sx q[1];
rz(2.286639) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236587) q[0];
sx q[0];
rz(-0.068327986) q[0];
sx q[0];
rz(-0.023850993) q[0];
rz(-pi) q[1];
rz(0.68944745) q[2];
sx q[2];
rz(-1.9580286) q[2];
sx q[2];
rz(2.7508798) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42306976) q[1];
sx q[1];
rz(-0.81941019) q[1];
sx q[1];
rz(-1.8526444) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1257954) q[3];
sx q[3];
rz(-1.7966174) q[3];
sx q[3];
rz(-2.8005536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.052224) q[2];
sx q[2];
rz(-1.3244018) q[2];
sx q[2];
rz(1.8734056) q[2];
rz(-0.44102272) q[3];
sx q[3];
rz(-2.5206168) q[3];
sx q[3];
rz(-0.84180251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0190401) q[0];
sx q[0];
rz(-2.1903116) q[0];
sx q[0];
rz(2.4265491) q[0];
rz(2.7252281) q[1];
sx q[1];
rz(-0.48960296) q[1];
sx q[1];
rz(2.8474836) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9884777) q[0];
sx q[0];
rz(-1.5128969) q[0];
sx q[0];
rz(2.2672976) q[0];
x q[1];
rz(-1.7753698) q[2];
sx q[2];
rz(-1.507221) q[2];
sx q[2];
rz(0.83067375) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3336368) q[1];
sx q[1];
rz(-0.2273493) q[1];
sx q[1];
rz(-0.61509404) q[1];
x q[2];
rz(0.58737602) q[3];
sx q[3];
rz(-1.0738292) q[3];
sx q[3];
rz(0.20020974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2405582) q[2];
sx q[2];
rz(-1.7419387) q[2];
sx q[2];
rz(-2.412839) q[2];
rz(-2.6567904) q[3];
sx q[3];
rz(-2.1383643) q[3];
sx q[3];
rz(-2.9639967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6840376) q[0];
sx q[0];
rz(-1.8192679) q[0];
sx q[0];
rz(-2.8628602) q[0];
rz(-0.067525603) q[1];
sx q[1];
rz(-0.7154671) q[1];
sx q[1];
rz(-0.50594893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26209575) q[0];
sx q[0];
rz(-0.56096062) q[0];
sx q[0];
rz(3.0039588) q[0];
rz(-pi) q[1];
rz(2.9437482) q[2];
sx q[2];
rz(-1.6704419) q[2];
sx q[2];
rz(-2.0951592) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3619574) q[1];
sx q[1];
rz(-1.0655626) q[1];
sx q[1];
rz(-2.510406) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2413195) q[3];
sx q[3];
rz(-2.4412324) q[3];
sx q[3];
rz(0.86802378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94023306) q[2];
sx q[2];
rz(-1.6910005) q[2];
sx q[2];
rz(0.10227164) q[2];
rz(-2.8376132) q[3];
sx q[3];
rz(-0.18054466) q[3];
sx q[3];
rz(-2.6612018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893696) q[0];
sx q[0];
rz(-0.82813534) q[0];
sx q[0];
rz(-2.4612259) q[0];
rz(2.1807189) q[1];
sx q[1];
rz(-1.5870794) q[1];
sx q[1];
rz(0.56112498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17616464) q[0];
sx q[0];
rz(-2.7211004) q[0];
sx q[0];
rz(-1.0499369) q[0];
x q[1];
rz(1.8041925) q[2];
sx q[2];
rz(-0.73783685) q[2];
sx q[2];
rz(2.3349175) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42237709) q[1];
sx q[1];
rz(-1.2724306) q[1];
sx q[1];
rz(1.5136615) q[1];
x q[2];
rz(0.64294358) q[3];
sx q[3];
rz(-2.3920632) q[3];
sx q[3];
rz(-1.7484322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9095824) q[2];
sx q[2];
rz(-0.97753111) q[2];
sx q[2];
rz(1.5781461) q[2];
rz(1.0163418) q[3];
sx q[3];
rz(-1.7137824) q[3];
sx q[3];
rz(-2.0411172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72961724) q[0];
sx q[0];
rz(-0.53891861) q[0];
sx q[0];
rz(-2.6475661) q[0];
rz(-1.5771075) q[1];
sx q[1];
rz(-0.99948519) q[1];
sx q[1];
rz(-3.0513501) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099797332) q[0];
sx q[0];
rz(-1.3372412) q[0];
sx q[0];
rz(0.44333027) q[0];
x q[1];
rz(-2.162287) q[2];
sx q[2];
rz(-2.3537209) q[2];
sx q[2];
rz(0.56172127) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68914946) q[1];
sx q[1];
rz(-1.1691368) q[1];
sx q[1];
rz(-2.8287877) q[1];
rz(0.19321038) q[3];
sx q[3];
rz(-1.8487559) q[3];
sx q[3];
rz(-1.4708468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6287441) q[2];
sx q[2];
rz(-2.5675842) q[2];
sx q[2];
rz(2.4981456) q[2];
rz(2.511034) q[3];
sx q[3];
rz(-2.3876811) q[3];
sx q[3];
rz(-3.0923617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028320463) q[0];
sx q[0];
rz(-1.7009108) q[0];
sx q[0];
rz(-1.2233618) q[0];
rz(0.89448294) q[1];
sx q[1];
rz(-2.5136785) q[1];
sx q[1];
rz(1.9953413) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47726163) q[0];
sx q[0];
rz(-2.1210175) q[0];
sx q[0];
rz(-0.47622891) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9983534) q[2];
sx q[2];
rz(-1.9263066) q[2];
sx q[2];
rz(-1.3035139) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.24487296) q[1];
sx q[1];
rz(-1.280652) q[1];
sx q[1];
rz(1.822132) q[1];
x q[2];
rz(-0.2616997) q[3];
sx q[3];
rz(-2.1402121) q[3];
sx q[3];
rz(-2.4401995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2724096) q[2];
sx q[2];
rz(-2.2792008) q[2];
sx q[2];
rz(1.4746846) q[2];
rz(2.9278751) q[3];
sx q[3];
rz(-1.4054047) q[3];
sx q[3];
rz(0.86764446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2992582) q[0];
sx q[0];
rz(-0.61115757) q[0];
sx q[0];
rz(2.4364731) q[0];
rz(-2.5662388) q[1];
sx q[1];
rz(-1.4082785) q[1];
sx q[1];
rz(-2.5720678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8370767) q[0];
sx q[0];
rz(-1.6479371) q[0];
sx q[0];
rz(-1.6521554) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6492278) q[2];
sx q[2];
rz(-1.8316744) q[2];
sx q[2];
rz(0.82792574) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.25334206) q[1];
sx q[1];
rz(-2.1565795) q[1];
sx q[1];
rz(0.82510494) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6105004) q[3];
sx q[3];
rz(-2.5887024) q[3];
sx q[3];
rz(-3.1097041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9465955) q[2];
sx q[2];
rz(-2.1985998) q[2];
sx q[2];
rz(-0.97879624) q[2];
rz(0.43092522) q[3];
sx q[3];
rz(-0.48723358) q[3];
sx q[3];
rz(-2.0144958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2836714) q[0];
sx q[0];
rz(-3.1223174) q[0];
sx q[0];
rz(-1.6790144) q[0];
rz(-1.0390394) q[1];
sx q[1];
rz(-1.7580527) q[1];
sx q[1];
rz(-1.2079814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20328612) q[0];
sx q[0];
rz(-1.948256) q[0];
sx q[0];
rz(-2.8164144) q[0];
rz(-1.6287491) q[2];
sx q[2];
rz(-0.52554146) q[2];
sx q[2];
rz(1.7271822) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9644248) q[1];
sx q[1];
rz(-0.96098122) q[1];
sx q[1];
rz(-2.9637106) q[1];
rz(-pi) q[2];
x q[2];
rz(1.966217) q[3];
sx q[3];
rz(-1.6334157) q[3];
sx q[3];
rz(-2.5015802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.98564467) q[2];
sx q[2];
rz(-2.1375103) q[2];
sx q[2];
rz(1.8211625) q[2];
rz(0.350658) q[3];
sx q[3];
rz(-0.81736332) q[3];
sx q[3];
rz(-0.58026522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.25528) q[0];
sx q[0];
rz(-1.1863615) q[0];
sx q[0];
rz(-0.86082389) q[0];
rz(1.6571922) q[1];
sx q[1];
rz(-1.7488372) q[1];
sx q[1];
rz(-1.7120672) q[1];
rz(-0.37921956) q[2];
sx q[2];
rz(-1.6368964) q[2];
sx q[2];
rz(1.5131769) q[2];
rz(1.9324393) q[3];
sx q[3];
rz(-2.4222838) q[3];
sx q[3];
rz(-2.1771976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
