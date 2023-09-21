OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(7.6927778) q[0];
sx q[0];
rz(11.132244) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(-2.5399962) q[1];
sx q[1];
rz(-0.4184202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2375862) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(2.0531274) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2566968) q[2];
sx q[2];
rz(-0.17527097) q[2];
sx q[2];
rz(-1.0980609) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0093706) q[1];
sx q[1];
rz(-1.8144061) q[1];
sx q[1];
rz(-1.1019215) q[1];
rz(-pi) q[2];
rz(2.0099785) q[3];
sx q[3];
rz(-1.7712799) q[3];
sx q[3];
rz(0.85103121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.819954) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(2.3036172) q[2];
rz(2.6485802) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943966) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(-1.863742) q[0];
rz(2.9648119) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(-2.7094254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85053274) q[0];
sx q[0];
rz(-1.2948372) q[0];
sx q[0];
rz(0.54131298) q[0];
rz(-pi) q[1];
rz(0.400153) q[2];
sx q[2];
rz(-1.1740985) q[2];
sx q[2];
rz(-0.19043365) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1951616) q[1];
sx q[1];
rz(-0.45743194) q[1];
sx q[1];
rz(-1.3034348) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1167332) q[3];
sx q[3];
rz(-2.1356574) q[3];
sx q[3];
rz(2.3853175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(0.48669997) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040722672) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(-1.4886645) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(0.506385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15784141) q[0];
sx q[0];
rz(-0.39641532) q[0];
sx q[0];
rz(2.0435964) q[0];
x q[1];
rz(0.26489139) q[2];
sx q[2];
rz(-2.3893917) q[2];
sx q[2];
rz(2.0843992) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80458927) q[1];
sx q[1];
rz(-2.5924006) q[1];
sx q[1];
rz(0.6188436) q[1];
rz(-2.0154325) q[3];
sx q[3];
rz(-2.7362842) q[3];
sx q[3];
rz(1.0682378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4425519) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(-2.5555723) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716924) q[0];
sx q[0];
rz(-0.62830347) q[0];
sx q[0];
rz(-0.83918321) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(-1.5930088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4913113) q[0];
sx q[0];
rz(-0.45931739) q[0];
sx q[0];
rz(-1.3643054) q[0];
rz(-0.71009212) q[2];
sx q[2];
rz(-2.2033764) q[2];
sx q[2];
rz(1.3131504) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3241987) q[1];
sx q[1];
rz(-0.84016582) q[1];
sx q[1];
rz(2.6170931) q[1];
x q[2];
rz(-1.2876835) q[3];
sx q[3];
rz(-1.1554171) q[3];
sx q[3];
rz(-2.0302041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8362391) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-0.099686064) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(-1.3249741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754159) q[0];
sx q[0];
rz(-1.7594936) q[0];
sx q[0];
rz(2.8856522) q[0];
rz(-0.4610962) q[1];
sx q[1];
rz(-2.0979116) q[1];
sx q[1];
rz(-0.76006132) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81239031) q[0];
sx q[0];
rz(-2.985552) q[0];
sx q[0];
rz(0.72326707) q[0];
rz(-3.0779482) q[2];
sx q[2];
rz(-1.9846989) q[2];
sx q[2];
rz(2.6457583) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8885986) q[1];
sx q[1];
rz(-2.9610486) q[1];
sx q[1];
rz(-2.351159) q[1];
rz(-pi) q[2];
rz(1.101196) q[3];
sx q[3];
rz(-2.4393775) q[3];
sx q[3];
rz(1.5576253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(-2.467353) q[2];
rz(-2.9267866) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(0.017344346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53133416) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.1791139) q[0];
rz(2.9367661) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(1.0669605) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93766312) q[0];
sx q[0];
rz(-1.5852889) q[0];
sx q[0];
rz(0.020676215) q[0];
rz(-pi) q[1];
rz(2.3382171) q[2];
sx q[2];
rz(-2.5640045) q[2];
sx q[2];
rz(-2.9327649) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8384335) q[1];
sx q[1];
rz(-1.618297) q[1];
sx q[1];
rz(0.82364239) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7339891) q[3];
sx q[3];
rz(-2.3605151) q[3];
sx q[3];
rz(1.2490602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.431488) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(-0.7263178) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85754919) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(-1.42111) q[0];
rz(0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(0.85817671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9439529) q[0];
sx q[0];
rz(-1.7174935) q[0];
sx q[0];
rz(-0.12978817) q[0];
rz(1.1602976) q[2];
sx q[2];
rz(-1.5294642) q[2];
sx q[2];
rz(1.2683887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3837636) q[1];
sx q[1];
rz(-1.541829) q[1];
sx q[1];
rz(1.5860735) q[1];
rz(-pi) q[2];
rz(-1.9167561) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(0.35017761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.28785607) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(-2.251513) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(-0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778075) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(-2.7668787) q[0];
rz(-2.162714) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(-1.7920866) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6222854) q[0];
sx q[0];
rz(-1.1847727) q[0];
sx q[0];
rz(0.65667721) q[0];
rz(-2.4235054) q[2];
sx q[2];
rz(-1.5880843) q[2];
sx q[2];
rz(-1.7503439) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6264682) q[1];
sx q[1];
rz(-1.6203468) q[1];
sx q[1];
rz(-1.7527761) q[1];
x q[2];
rz(-1.3091062) q[3];
sx q[3];
rz(-2.1013386) q[3];
sx q[3];
rz(-1.7390651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.65138856) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(0.46869579) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(1.3635925) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63012183) q[0];
sx q[0];
rz(-2.2352495) q[0];
sx q[0];
rz(1.2619031) q[0];
rz(-2.966554) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(-1.5375686) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984734) q[0];
sx q[0];
rz(-2.7612918) q[0];
sx q[0];
rz(0.93912504) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41861694) q[2];
sx q[2];
rz(-1.4756243) q[2];
sx q[2];
rz(-1.0375432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6919842) q[1];
sx q[1];
rz(-0.87032986) q[1];
sx q[1];
rz(1.7193754) q[1];
rz(-pi) q[2];
rz(-2.780464) q[3];
sx q[3];
rz(-1.1361406) q[3];
sx q[3];
rz(1.1772732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1016772) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(2.3802479) q[2];
rz(0.90041655) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(0.049023978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3392357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(0.21690579) q[0];
rz(-0.63198173) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(2.1868618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.117884) q[0];
sx q[0];
rz(-2.3291991) q[0];
sx q[0];
rz(2.6931767) q[0];
rz(-1.9913313) q[2];
sx q[2];
rz(-2.0805801) q[2];
sx q[2];
rz(-1.2812986) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2659457) q[1];
sx q[1];
rz(-1.6469643) q[1];
sx q[1];
rz(2.6850558) q[1];
rz(-0.44317742) q[3];
sx q[3];
rz(-2.157353) q[3];
sx q[3];
rz(-3.0122258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1251936) q[2];
sx q[2];
rz(-1.23896) q[2];
sx q[2];
rz(0.94474244) q[2];
rz(2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(-0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5007297) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-0.8846994) q[1];
sx q[1];
rz(-2.2365166) q[1];
sx q[1];
rz(2.8832163) q[1];
rz(2.2434071) q[2];
sx q[2];
rz(-0.85815103) q[2];
sx q[2];
rz(1.382538) q[2];
rz(-2.9363019) q[3];
sx q[3];
rz(-2.0060354) q[3];
sx q[3];
rz(2.5397186) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];