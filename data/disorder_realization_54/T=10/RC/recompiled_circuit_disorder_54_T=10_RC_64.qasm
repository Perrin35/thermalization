OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(2.8110992) q[0];
rz(-2.8117872) q[1];
sx q[1];
rz(-2.2916315) q[1];
sx q[1];
rz(-0.70911521) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5753484) q[0];
sx q[0];
rz(-2.4465804) q[0];
sx q[0];
rz(2.038875) q[0];
rz(3.1094413) q[2];
sx q[2];
rz(-1.2574147) q[2];
sx q[2];
rz(-1.1653479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46804894) q[1];
sx q[1];
rz(-1.1357726) q[1];
sx q[1];
rz(0.97462868) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0508918) q[3];
sx q[3];
rz(-1.2399779) q[3];
sx q[3];
rz(-2.1045121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7314529) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(1.6072134) q[2];
rz(-2.2065227) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(-2.5193135) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(-0.91631779) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4692357) q[0];
sx q[0];
rz(-1.0807481) q[0];
sx q[0];
rz(2.310918) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3071438) q[2];
sx q[2];
rz(-1.4730028) q[2];
sx q[2];
rz(-2.9546839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2434477) q[1];
sx q[1];
rz(-1.2444278) q[1];
sx q[1];
rz(1.9279724) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.968859) q[3];
sx q[3];
rz(-1.5916087) q[3];
sx q[3];
rz(1.2956937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(0.82143482) q[2];
rz(0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538552) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(2.1333372) q[0];
rz(0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(-0.52454138) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7099972) q[0];
sx q[0];
rz(-1.734373) q[0];
sx q[0];
rz(3.1113935) q[0];
rz(-pi) q[1];
rz(-1.3069986) q[2];
sx q[2];
rz(-2.2321777) q[2];
sx q[2];
rz(-1.7922572) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9959065) q[1];
sx q[1];
rz(-1.6325163) q[1];
sx q[1];
rz(1.5476336) q[1];
rz(-1.1509622) q[3];
sx q[3];
rz(-1.435278) q[3];
sx q[3];
rz(-2.6672222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(-1.4952205) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33048531) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(0.088949732) q[0];
rz(-0.51070172) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(-2.451992) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423858) q[0];
sx q[0];
rz(-0.33873522) q[0];
sx q[0];
rz(0.50868209) q[0];
rz(-1.7961411) q[2];
sx q[2];
rz(-2.1225404) q[2];
sx q[2];
rz(1.4622886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1414813) q[1];
sx q[1];
rz(-0.82427102) q[1];
sx q[1];
rz(0.4517171) q[1];
x q[2];
rz(-0.48456405) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66578635) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(0.62292567) q[2];
rz(-2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2899807) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(0.583453) q[0];
rz(-1.9955697) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(-1.4978283) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9153584) q[0];
sx q[0];
rz(-0.81904531) q[0];
sx q[0];
rz(-0.35924964) q[0];
rz(-pi) q[1];
rz(-3.0450902) q[2];
sx q[2];
rz(-0.8000024) q[2];
sx q[2];
rz(-0.9347136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6795265) q[1];
sx q[1];
rz(-0.70409173) q[1];
sx q[1];
rz(-0.042298869) q[1];
rz(-pi) q[2];
rz(2.6260914) q[3];
sx q[3];
rz(-1.4294129) q[3];
sx q[3];
rz(0.055671234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.65537611) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(-1.5779457) q[2];
rz(-2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435796) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-0.046982732) q[0];
rz(2.9934096) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(-1.7061589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3008266) q[0];
sx q[0];
rz(-0.78221417) q[0];
sx q[0];
rz(-1.7813111) q[0];
rz(-3.0827423) q[2];
sx q[2];
rz(-2.119679) q[2];
sx q[2];
rz(0.6097874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6378577) q[1];
sx q[1];
rz(-1.6956455) q[1];
sx q[1];
rz(1.1605074) q[1];
rz(-pi) q[2];
rz(1.8085603) q[3];
sx q[3];
rz(-1.5974345) q[3];
sx q[3];
rz(0.8386855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(-3.0701239) q[2];
rz(-1.4525157) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(-0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(0.57762161) q[0];
rz(-1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(-2.1320027) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6972835) q[0];
sx q[0];
rz(-2.5464006) q[0];
sx q[0];
rz(2.5213581) q[0];
rz(-2.6206714) q[2];
sx q[2];
rz(-0.43653566) q[2];
sx q[2];
rz(1.6233363) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2958888) q[1];
sx q[1];
rz(-1.0755952) q[1];
sx q[1];
rz(1.1843029) q[1];
x q[2];
rz(2.1647251) q[3];
sx q[3];
rz(-0.90126029) q[3];
sx q[3];
rz(0.23564786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0916831) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(0.47232929) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(1.002731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.49884477) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(-2.3587976) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(-0.4447287) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3702104) q[0];
sx q[0];
rz(-1.0078197) q[0];
sx q[0];
rz(0.51008205) q[0];
rz(-pi) q[1];
rz(2.2112446) q[2];
sx q[2];
rz(-2.1313416) q[2];
sx q[2];
rz(2.408574) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6726482) q[1];
sx q[1];
rz(-2.4441507) q[1];
sx q[1];
rz(-1.7698606) q[1];
x q[2];
rz(-2.7304857) q[3];
sx q[3];
rz(-1.6931705) q[3];
sx q[3];
rz(1.7468332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-0.81531173) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.8243272) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8274882) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-2.8544193) q[0];
rz(-2.9526967) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(2.8093991) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4591601) q[0];
sx q[0];
rz(-0.81572616) q[0];
sx q[0];
rz(2.0072323) q[0];
x q[1];
rz(-0.046499649) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(-2.7547714) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6576782) q[1];
sx q[1];
rz(-2.6473443) q[1];
sx q[1];
rz(2.5792522) q[1];
rz(1.166415) q[3];
sx q[3];
rz(-2.0669193) q[3];
sx q[3];
rz(-2.8465084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(1.8509289) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(0.65657842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0861417) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(2.0822051) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(0.25451452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33728889) q[0];
sx q[0];
rz(-2.8992607) q[0];
sx q[0];
rz(-1.2605577) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6539831) q[2];
sx q[2];
rz(-1.3541823) q[2];
sx q[2];
rz(-2.285514) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3855615) q[1];
sx q[1];
rz(-0.33547151) q[1];
sx q[1];
rz(-0.24753333) q[1];
rz(-pi) q[2];
rz(-1.47154) q[3];
sx q[3];
rz(-1.1343065) q[3];
sx q[3];
rz(1.0222767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.09482) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(2.0937031) q[2];
rz(-1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(2.5361983) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142674) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(-2.7813773) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(0.24047492) q[2];
sx q[2];
rz(-2.0839811) q[2];
sx q[2];
rz(2.6032084) q[2];
rz(-1.8923106) q[3];
sx q[3];
rz(-0.63890639) q[3];
sx q[3];
rz(-1.3756868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
