OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11256448) q[0];
sx q[0];
rz(-1.636314) q[0];
sx q[0];
rz(0.78328744) q[0];
rz(1.0822436) q[1];
sx q[1];
rz(-1.487027) q[1];
sx q[1];
rz(-2.3496871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5675674) q[0];
sx q[0];
rz(-1.8301395) q[0];
sx q[0];
rz(1.7631084) q[0];
rz(-pi) q[1];
rz(1.7246805) q[2];
sx q[2];
rz(-1.6072011) q[2];
sx q[2];
rz(-0.77637451) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.344393) q[1];
sx q[1];
rz(-2.5973592) q[1];
sx q[1];
rz(3.0809156) q[1];
x q[2];
rz(-0.79444076) q[3];
sx q[3];
rz(-2.0047965) q[3];
sx q[3];
rz(1.3484914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59114328) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(0.099451065) q[2];
rz(-0.24672306) q[3];
sx q[3];
rz(-1.5244502) q[3];
sx q[3];
rz(2.7439086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1219015) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(0.9285399) q[0];
rz(-1.6909201) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(-2.4931152) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4546616) q[0];
sx q[0];
rz(-1.0333916) q[0];
sx q[0];
rz(-0.65673687) q[0];
rz(0.5669539) q[2];
sx q[2];
rz(-1.840072) q[2];
sx q[2];
rz(-0.93967162) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9916223) q[1];
sx q[1];
rz(-2.9987965) q[1];
sx q[1];
rz(-3.0317329) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2821372) q[3];
sx q[3];
rz(-1.4573754) q[3];
sx q[3];
rz(2.3577549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6917981) q[2];
sx q[2];
rz(-0.75642502) q[2];
sx q[2];
rz(2.2890384) q[2];
rz(-0.84613386) q[3];
sx q[3];
rz(-1.6328014) q[3];
sx q[3];
rz(-1.030863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5071252) q[0];
sx q[0];
rz(-1.9159303) q[0];
sx q[0];
rz(-1.0852098) q[0];
rz(0.94404864) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(1.5012213) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1000993) q[0];
sx q[0];
rz(-1.3736808) q[0];
sx q[0];
rz(-1.5310578) q[0];
rz(-0.88884647) q[2];
sx q[2];
rz(-2.8434128) q[2];
sx q[2];
rz(-2.2928638) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2861917) q[1];
sx q[1];
rz(-0.58127379) q[1];
sx q[1];
rz(-1.3894677) q[1];
rz(0.15764938) q[3];
sx q[3];
rz(-1.7846037) q[3];
sx q[3];
rz(-2.714954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40470716) q[2];
sx q[2];
rz(-2.6077304) q[2];
sx q[2];
rz(-0.85477465) q[2];
rz(0.9084304) q[3];
sx q[3];
rz(-1.5769438) q[3];
sx q[3];
rz(2.0250208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.214355) q[0];
sx q[0];
rz(-0.4011918) q[0];
sx q[0];
rz(-1.1965363) q[0];
rz(1.475097) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(1.9570785) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9036983) q[0];
sx q[0];
rz(-1.8082976) q[0];
sx q[0];
rz(0.71126513) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75293137) q[2];
sx q[2];
rz(-1.4309466) q[2];
sx q[2];
rz(-1.1876196) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2503916) q[1];
sx q[1];
rz(-1.7970835) q[1];
sx q[1];
rz(1.0550642) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58301894) q[3];
sx q[3];
rz(-0.32430092) q[3];
sx q[3];
rz(0.7079269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.92675942) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(1.8801749) q[2];
rz(2.7091806) q[3];
sx q[3];
rz(-2.0093446) q[3];
sx q[3];
rz(-2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686907) q[0];
sx q[0];
rz(-1.2304767) q[0];
sx q[0];
rz(-0.029408971) q[0];
rz(-2.3853761) q[1];
sx q[1];
rz(-2.5543946) q[1];
sx q[1];
rz(2.3888033) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22827029) q[0];
sx q[0];
rz(-3.1359657) q[0];
sx q[0];
rz(1.6371284) q[0];
rz(-pi) q[1];
rz(1.2529066) q[2];
sx q[2];
rz(-0.28159062) q[2];
sx q[2];
rz(1.3791305) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.237287) q[1];
sx q[1];
rz(-2.5616966) q[1];
sx q[1];
rz(-2.94983) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4642508) q[3];
sx q[3];
rz(-0.65653446) q[3];
sx q[3];
rz(-1.7650106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0011562) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(-2.746554) q[2];
rz(0.47518528) q[3];
sx q[3];
rz(-1.3522215) q[3];
sx q[3];
rz(0.37676677) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659371) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(0.96250594) q[0];
rz(2.6267701) q[1];
sx q[1];
rz(-1.6074901) q[1];
sx q[1];
rz(0.33822507) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048934919) q[0];
sx q[0];
rz(-2.3065164) q[0];
sx q[0];
rz(0.2661163) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8184616) q[2];
sx q[2];
rz(-1.3405352) q[2];
sx q[2];
rz(-2.3375653) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.6392614) q[1];
sx q[1];
rz(-0.0031454589) q[1];
sx q[1];
rz(-1.5793213) q[1];
rz(1.148726) q[3];
sx q[3];
rz(-1.633612) q[3];
sx q[3];
rz(-0.40336762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6414791) q[2];
sx q[2];
rz(-1.9794455) q[2];
sx q[2];
rz(-0.55434736) q[2];
rz(2.2902299) q[3];
sx q[3];
rz(-0.35836372) q[3];
sx q[3];
rz(0.62801492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8165269) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(-2.391173) q[0];
rz(0.57506192) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(-2.4837928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039666273) q[0];
sx q[0];
rz(-2.3453379) q[0];
sx q[0];
rz(0.58923652) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70722981) q[2];
sx q[2];
rz(-1.6042738) q[2];
sx q[2];
rz(1.0302803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9896587) q[1];
sx q[1];
rz(-2.8550145) q[1];
sx q[1];
rz(1.6295939) q[1];
rz(-pi) q[2];
rz(-2.9255387) q[3];
sx q[3];
rz(-1.9441868) q[3];
sx q[3];
rz(1.7301529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49360069) q[2];
sx q[2];
rz(-1.0085663) q[2];
sx q[2];
rz(1.4823401) q[2];
rz(-3.0209387) q[3];
sx q[3];
rz(-1.3841265) q[3];
sx q[3];
rz(-0.83546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7928612) q[0];
sx q[0];
rz(-2.272235) q[0];
sx q[0];
rz(0.29801512) q[0];
rz(-1.0385849) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-0.15377741) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2033635) q[0];
sx q[0];
rz(-1.3171609) q[0];
sx q[0];
rz(1.6481736) q[0];
rz(0.91570274) q[2];
sx q[2];
rz(-0.88353744) q[2];
sx q[2];
rz(-2.5926431) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.01659) q[1];
sx q[1];
rz(-1.9167797) q[1];
sx q[1];
rz(-1.6660652) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3078717) q[3];
sx q[3];
rz(-2.701303) q[3];
sx q[3];
rz(1.057098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.095857233) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(-1.9528961) q[2];
rz(-1.3304322) q[3];
sx q[3];
rz(-1.9537787) q[3];
sx q[3];
rz(2.4082898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34340149) q[0];
sx q[0];
rz(-2.5372086) q[0];
sx q[0];
rz(-1.4991624) q[0];
rz(-2.5921953) q[1];
sx q[1];
rz(-1.5645942) q[1];
sx q[1];
rz(0.25064358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3651705) q[0];
sx q[0];
rz(-0.93342268) q[0];
sx q[0];
rz(-2.754162) q[0];
rz(2.8806503) q[2];
sx q[2];
rz(-1.5847209) q[2];
sx q[2];
rz(-1.5858638) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.033015164) q[1];
sx q[1];
rz(-2.3653154) q[1];
sx q[1];
rz(-3.0423711) q[1];
rz(-1.4446684) q[3];
sx q[3];
rz(-1.3775004) q[3];
sx q[3];
rz(0.85053683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9833019) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(-2.1251202) q[2];
rz(3.0554092) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(-2.3763954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8330399) q[0];
sx q[0];
rz(-0.85423952) q[0];
sx q[0];
rz(-2.7225851) q[0];
rz(2.6047756) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(0.73582617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2013071) q[0];
sx q[0];
rz(-2.2129411) q[0];
sx q[0];
rz(1.8117732) q[0];
rz(-pi) q[1];
rz(0.67086733) q[2];
sx q[2];
rz(-0.76422526) q[2];
sx q[2];
rz(0.66167458) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23259737) q[1];
sx q[1];
rz(-2.4687662) q[1];
sx q[1];
rz(0.35079591) q[1];
x q[2];
rz(2.5352468) q[3];
sx q[3];
rz(-2.5744573) q[3];
sx q[3];
rz(2.4727269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95911038) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(-2.518892) q[2];
rz(0.73838082) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(-0.27899376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57711346) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(2.4979757) q[1];
sx q[1];
rz(-1.2871965) q[1];
sx q[1];
rz(1.0866477) q[1];
rz(-2.9109091) q[2];
sx q[2];
rz(-1.2103969) q[2];
sx q[2];
rz(-2.8878676) q[2];
rz(-2.6162061) q[3];
sx q[3];
rz(-2.864884) q[3];
sx q[3];
rz(-1.3340193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
